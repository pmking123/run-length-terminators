#!/usr/bin/env python3
"""bootstrap_enrichment.py

Genus-clustered bootstrap + leave-one-genus-out for gene-end tract enrichment.

Reads per_candidate_by_k.csv (accession, genus, family, k, tract_len, gene_end,
null_rate) as emitted by clade_enrichment_table.py. The matched null is held fixed at
its published per-k value; only the candidate tracts are resampled, clustered by GENUS,
which answers the reviewer question "is the enrichment driven by dense sampling within
one genus?" without needing per-genome null counts (which the summary does not expose).

For each requested k it prints the point enrichment, the genus-clustered 95% CI, the
max single-genus candidate share, and the full leave-one-genus-out table, and it writes
a durable CSV (one row per k) so the intervals are regenerable rather than living only
in terminal scrollback.

Usage
-----
  # all thresholds present for the family, in one run:
  python bootstrap_enrichment.py --per-candidate .../per_candidate_by_k.csv --family Rhabdoviridae
  # or a specific k / subset:
  python bootstrap_enrichment.py --per-candidate .../per_candidate_by_k.csv --family Rhabdoviridae --k 6,7

Outputs
-------
  <outdir or per-candidate dir>/<family>_enrichment_ci.csv
      family, k, n_candidates, n_gene_end, n_genera, max_genus_share,
      enrichment, ci_lo, ci_hi, logo_min, logo_max, B, seed
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def enrichment(df, null_rate):
    n = len(df)
    if n == 0:
        return np.nan
    return df["gene_end"].mean() / null_rate      # frac / null, matching the frozen table


def genus_cluster_bootstrap(df, null_rate, B=4000, seed=0):
    rng = np.random.default_rng(seed)
    genera = df["genus"].unique()
    by_genus = {g: df[df.genus == g] for g in genera}     # pre-split once, not per draw
    point = enrichment(df, null_rate)
    draws = np.empty(B)
    for i in range(B):
        pick = rng.choice(genera, size=len(genera), replace=True)
        boot = pd.concat([by_genus[g] for g in pick], ignore_index=True)
        draws[i] = enrichment(boot, null_rate)
    lo, hi = np.nanpercentile(draws, [2.5, 97.5])
    return point, lo, hi, len(genera)


def leave_one_genus_out(df, null_rate):
    genera = df["genus"].unique()
    out = [(g, int((df.genus == g).sum()), enrichment(df[df.genus != g], null_rate))
           for g in genera]
    return pd.DataFrame(out, columns=["dropped_genus", "n_candidates_dropped",
                                      "enrichment_without"]).sort_values("enrichment_without")


def run_one_k(df_k, k, family, B, seed):
    null_rate = df_k["null_rate"].iloc[0]
    assert df_k["null_rate"].nunique() == 1, f"null_rate not constant within k={k}"
    point, lo, hi, ng = genus_cluster_bootstrap(df_k, null_rate, B=B, seed=seed)
    logo = leave_one_genus_out(df_k, null_rate)
    max_share = int(df_k.groupby("genus").size().max())

    print(f"{family}  k>={k}   ({len(df_k)} candidates, "
          f"{int(df_k.gene_end.sum())} gene-end, {ng} genera)")
    print(f"  enrichment          : {point:.1f}x")
    print(f"  genus-cluster 95% CI: [{lo:.1f}x, {hi:.1f}x]   (B={B})")
    print(f"  max single-genus share: {max_share}/{len(df_k)} candidates")
    print(f"  leave-one-genus-out band: [{logo.enrichment_without.min():.1f}x, "
          f"{logo.enrichment_without.max():.1f}x]")
    print("  leave-one-genus-out (sorted; top = removal lowers most):")
    print(logo.to_string(index=False, float_format=lambda v: f"{v:.1f}"))
    print()

    return dict(family=family, k=k, n_candidates=len(df_k),
                n_gene_end=int(df_k.gene_end.sum()), n_genera=ng,
                max_genus_share=max_share, enrichment=round(point, 3),
                ci_lo=round(lo, 3), ci_hi=round(hi, 3),
                logo_min=round(logo.enrichment_without.min(), 3),
                logo_max=round(logo.enrichment_without.max(), 3), B=B, seed=seed)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-candidate", required=True, type=Path)
    ap.add_argument("--family", default="Rhabdoviridae")
    ap.add_argument("--k", default=None,
                    help="comma-separated k values (default: all present for the family)")
    ap.add_argument("--B", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--outdir", default=None, type=Path,
                    help="where to write <family>_enrichment_ci.csv "
                         "(default: alongside --per-candidate)")
    args = ap.parse_args()

    df = pd.read_csv(args.per_candidate)
    df = df[df.family == args.family].copy()
    if df.empty:
        raise SystemExit(f"no candidates for family={args.family}")

    if args.k:
        ks = [int(x) for x in args.k.split(",")]
    else:
        ks = sorted(int(k) for k in df.k.unique())

    records = []
    for k in ks:
        df_k = df[df.k == k]
        if df_k.empty:
            print(f"WARNING: no candidates for {args.family} k>={k}; skipping\n")
            continue
        records.append(run_one_k(df_k, k, args.family, args.B, args.seed))

    outdir = args.outdir or args.per_candidate.parent
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"{args.family}_enrichment_ci.csv"
    pd.DataFrame(records).to_csv(out, index=False)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
