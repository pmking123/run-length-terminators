#!/usr/bin/env python3
"""
per_genome_paired_test.py

The pseudoreplication control for Results 3.1, on the FROZEN run table.

Pooling millions of runs treats non-independent runs within a genome as independent
observations, so pooled p-values / z-scores are anti-conservative by orders of
magnitude. This script makes the GENOME the unit of inference:

  * for each genome, compute E_near = P(lambda >= k | d < NEAR) and the same quantity
    under two matched nulls:
        - uniform null : boundary positions resampled uniformly within the genome
        - rotation null: the boundary track circularly rotated by a random offset
          (preserves boundary count AND spacing AND the genome's compositional
           landscape; only the registration between runs and boundaries is broken)
  * compare observed vs null ACROSS genomes with a paired Wilcoxon signed-rank test
  * report the fraction of genomes with observed > null, and the fraction with
    per-genome odds ratio > 1.5

Also repeats the whole thing excluding boundary-overlapping (distance-0) runs.

Input : run_table.csv from boundary_enrichment.py (the frozen one)
Output: <outdir>/
          per_genome_paired_summary.csv     <- headline paired-test numbers
          per_genome_enrichment.csv         <- one row per genome (obs, nulls, OR)

Needs: pandas, numpy, scipy.  No network.

Example
-------
python per_genome_paired_test.py \
    --run-table analysis_final/boundary_enrichment/run_table.csv \
    --outdir analysis_final/per_genome_paired \
    --k 5 --near 25 --n-null 200 --seed 12345
"""
from __future__ import annotations
import argparse, re
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

BASES = ["A", "C", "G", "T"]
KIND_POS_RE = re.compile(r":(?:start|end):(\d+)")


def load(run_table: Path):
    usecols = ["accession", "sequence_length", "base", "run_length",
               "run_start_1based", "run_end_1based",
               "boundary_hits", "nearest_boundary"]
    accs, rl, s, e, is_at = [], [], [], [], []
    boundaries = defaultdict(set)
    seqlen = {}
    for chunk in pd.read_csv(run_table, usecols=usecols, dtype=str, chunksize=400_000):
        chunk = chunk[chunk["base"].isin(BASES)]
        a = chunk["accession"].to_numpy()
        accs.append(a)
        rl.append(pd.to_numeric(chunk["run_length"]).to_numpy())
        s.append(pd.to_numeric(chunk["run_start_1based"]).to_numpy())
        e.append(pd.to_numeric(chunk["run_end_1based"]).to_numpy())
        is_at.append(chunk["base"].isin(["A", "T"]).to_numpy())
        for acc, L in zip(a, pd.to_numeric(chunk["sequence_length"]).to_numpy()):
            seqlen.setdefault(acc, int(L))
        for acc, bh, nb in zip(a, chunk["boundary_hits"].fillna(""),
                               chunk["nearest_boundary"].fillna("")):
            for p in KIND_POS_RE.findall(bh):
                boundaries[acc].add(int(p))
            for p in KIND_POS_RE.findall(nb):
                boundaries[acc].add(int(p))
    df = pd.DataFrame({
        "acc": np.concatenate(accs),
        "rl": np.concatenate(rl).astype(np.int32),
        "s": np.concatenate(s).astype(np.int64),
        "e": np.concatenate(e).astype(np.int64),
        "at": np.concatenate(is_at),
    })
    return df, {a: np.array(sorted(b), dtype=np.int64) for a, b in boundaries.items()}, seqlen


def nearest_dist(starts, ends, b):
    if len(b) == 0:
        return np.full(len(starts), np.inf)
    lo = np.searchsorted(b, starts, "left")
    hi = np.searchsorted(b, ends, "right")
    contains = hi > lo
    below = np.where(lo > 0, starts - b[np.clip(lo - 1, 0, len(b) - 1)], np.inf)
    above = np.where(hi < len(b), b[np.clip(hi, 0, len(b) - 1)] - ends, np.inf)
    d = np.minimum(below, above)
    d[contains] = 0.0
    return d


def rate(ge, mask):
    return float(ge[mask].mean()) if mask.any() else np.nan


def odds_ratio(ge, near):
    a1 = int((ge & near).sum()); a0 = int((~ge & near).sum())
    c1 = int((ge & ~near).sum()); c0 = int((~ge & ~near).sum())
    if a0 == 0 or c1 == 0:
        return np.nan
    return (a1 * c0) / (a0 * c1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-table", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--near", type=int, default=25)
    ap.add_argument("--n-null", type=int, default=200,
                    help="null replicates per genome (uniform and rotation each)")
    ap.add_argument("--at-only", action="store_true",
                    help="restrict to A/T runs (the bases carrying the signal)")
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    print("loading run table + reconstructing boundaries ...", flush=True)
    df, B, seqlen = load(args.run_table)
    if args.at_only:
        df = df[df["at"]].copy()
    print(f"  {len(df):,} runs, {df['acc'].nunique()} genomes", flush=True)

    rows = []
    for acc, sub in df.groupby("acc", sort=False):
        b = B.get(acc, np.array([], dtype=np.int64))
        L = int(seqlen.get(acc, 0))
        if len(b) < 2 or L <= 0:
            continue
        s = sub["s"].to_numpy(); e = sub["e"].to_numpy()
        ge = (sub["rl"].to_numpy() >= args.k)
        if ge.sum() == 0:
            continue

        d = nearest_dist(s, e, b)
        near = d < args.near
        nz = d > 0                              # non-overlapping (width control)
        obs = rate(ge, near)
        obs_no = rate(ge[nz], near[nz]) if nz.any() else np.nan
        or_obs = odds_ratio(ge, near)

        un, rt, un_no, rt_no = [], [], [], []
        allpos = np.arange(1, L + 1)
        nb = min(len(b), L)
        for _ in range(args.n_null):
            bu = np.sort(rng.choice(allpos, size=nb, replace=False))
            du = nearest_dist(s, e, bu); nu = du < args.near
            un.append(rate(ge, nu))
            m = du > 0
            un_no.append(rate(ge[m], nu[m]) if m.any() else np.nan)

            off = int(rng.integers(1, L))
            br = np.sort(((b - 1 + off) % L) + 1)
            dr = nearest_dist(s, e, br); nr = dr < args.near
            rt.append(rate(ge, nr))
            m = dr > 0
            rt_no.append(rate(ge[m], nr[m]) if m.any() else np.nan)

        rows.append(dict(accession=acc, n_runs=len(sub), n_long=int(ge.sum()),
                         n_boundaries=len(b), seq_len=L,
                         E_near_obs=obs, E_near_uniform=np.nanmean(un),
                         E_near_rotation=np.nanmean(rt), OR_obs=or_obs,
                         E_near_obs_no_overlap=obs_no,
                         E_near_uniform_no_overlap=np.nanmean(un_no),
                         E_near_rotation_no_overlap=np.nanmean(rt_no)))
    per = pd.DataFrame(rows)
    per.to_csv(args.outdir / "per_genome_enrichment.csv", index=False)

    def paired(obs_col, null_col, label):
        m = per[[obs_col, null_col]].dropna()
        o = m[obs_col].to_numpy(); n = m[null_col].to_numpy()
        keep = (o != n)                       # Wilcoxon drops zero-differences anyway
        try:
            stat, p = wilcoxon(o[keep], n[keep])
        except Exception:
            stat, p = np.nan, np.nan
        return dict(comparison=label, n_genomes=int(len(m)),
                    median_obs=float(np.median(o)), median_null=float(np.median(n)),
                    median_diff=float(np.median(o - n)),
                    frac_obs_gt_null=float((o > n).mean()),
                    wilcoxon_p=float(p))

    summ = pd.DataFrame([
        paired("E_near_obs", "E_near_uniform",  "observed vs UNIFORM null"),
        paired("E_near_obs", "E_near_rotation", "observed vs ROTATION null"),
        paired("E_near_obs_no_overlap", "E_near_uniform_no_overlap",
               "observed vs UNIFORM null (overlap-excluded)"),
        paired("E_near_obs_no_overlap", "E_near_rotation_no_overlap",
               "observed vs ROTATION null (overlap-excluded)"),
    ])
    frac_or = float((per["OR_obs"] > 1.5).mean())
    summ["frac_genomes_OR_gt_1.5"] = frac_or
    summ.to_csv(args.outdir / "per_genome_paired_summary.csv", index=False)

    print("\n=== Per-genome paired tests (genome = unit of inference) ===")
    print(summ.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print(f"\nfraction of genomes with per-genome OR > 1.5: {frac_or:.3f}"
          f"  (n = {per['OR_obs'].notna().sum()} genomes)")
    print(f"\nWrote {args.outdir}/per_genome_paired_summary.csv and per_genome_enrichment.csv")


if __name__ == "__main__":
    main()
