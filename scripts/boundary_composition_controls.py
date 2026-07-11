#!/usr/bin/env python3
"""
boundary_composition_controls.py

Closes the composition-control gap in Results 3.1 / Table 1, on the FROZEN run table.

For long runs (default lambda >= 5) near vs far from CDS boundaries it reports, per
base and for all bases pooled:

  * crude near/far odds ratio
  * Mantel-Haenszel odds ratio adjusted for LOCAL base composition
        (stratified by decile of local A+T fraction in a +/-W window around each run,
         EXCLUDING the run's own bases; W is swept, default 15 and 50 bp)
  * the same, EXCLUDING runs that overlap a boundary (distance 0)  -> width control
  * local A+T fraction near vs far (the compositional shift itself)

and, separately, the same odds ratios split by STRAND-CORRECTED boundary kind, i.e.
biological 5' vs 3' gene ends (a coordinate "start"/"end" is a biological 5'/3' end
only after accounting for strand).

Everything is recomputed from the run table alone: boundary positions, kinds and
strands are parsed out of the `boundary_hits` / `nearest_boundary` text columns, so no
GenBank access or network is needed. Memory is kept flat by streaming the CSV.

Input : run_table.csv from boundary_enrichment.py (the frozen one)
Output: <outdir>/
          table1_boundary_enrichment.csv      <- Table 1
          table1_boundary_enrichment.md       <- same, markdown, paste-ready
          local_composition_near_far.csv      <- the A+T shift

Needs: pandas, numpy.  No network.

Example
-------
python boundary_composition_controls.py \
    --run-table analysis_final/boundary_enrichment/run_table.csv \
    --outdir analysis_final/boundary_composition_controls \
    --k 5 --near 25 --windows 15,50
"""
from __future__ import annotations
import argparse, re
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

BASES = ["A", "C", "G", "T"]
CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
# boundary_hits / nearest_boundary look like:  ...:start:1234;   ...:end:5678:feature=..:strand=+
KIND_POS_RE = re.compile(r":(start|end):(\d+)")
STRAND_RE = re.compile(r":(?:start|end):(\d+):feature=\d+-\d+:strand=([+\-])")
NDEC = 10


# --------------------------------------------------------------------- loading
def load(run_table: Path):
    """Stream the run table; return run rows + per-accession boundary kind/strand."""
    usecols = ["accession", "sequence_length", "base", "run_length",
               "run_start_1based", "run_end_1based",
               "boundary_hits", "nearest_boundary"]
    accs, code, rl, s, e = [], [], [], [], []
    kinds = defaultdict(lambda: defaultdict(set))   # acc -> pos -> {'start','end'}
    strands = defaultdict(dict)                     # acc -> pos -> '+'/'-'
    seqlen = {}
    for chunk in pd.read_csv(run_table, usecols=usecols, dtype=str, chunksize=400_000):
        chunk = chunk[chunk["base"].isin(BASES)]
        a = chunk["accession"].to_numpy()
        accs.append(a)
        code.append(chunk["base"].map(CODE).to_numpy())
        rl.append(pd.to_numeric(chunk["run_length"]).to_numpy())
        s.append(pd.to_numeric(chunk["run_start_1based"]).to_numpy())
        e.append(pd.to_numeric(chunk["run_end_1based"]).to_numpy())
        for acc, L in zip(a, pd.to_numeric(chunk["sequence_length"]).to_numpy()):
            seqlen.setdefault(acc, int(L))
        for acc, bh, nb in zip(a, chunk["boundary_hits"].fillna(""),
                               chunk["nearest_boundary"].fillna("")):
            for k, p in KIND_POS_RE.findall(bh):
                kinds[acc][int(p)].add(k)
            for k, p in KIND_POS_RE.findall(nb):
                kinds[acc][int(p)].add(k)
            for p, st in STRAND_RE.findall(nb):
                strands[acc][int(p)] = st
    df = pd.DataFrame({
        "acc": np.concatenate(accs),
        "code": np.concatenate(code).astype(np.int8),
        "rl": np.concatenate(rl).astype(np.int32),
        "s": np.concatenate(s).astype(np.int64),
        "e": np.concatenate(e).astype(np.int64),
    })
    return df, kinds, strands, seqlen


# ----------------------------------------------------------------- annotation
def annotate(df, kinds, strands, seqlen, windows):
    """Add distance, nearest-boundary biological end (5'/3'), and local A+T per window."""
    df = df.sort_values(["acc", "s"]).reset_index(drop=True)
    n = len(df)
    dist = np.full(n, np.inf)
    biol = np.empty(n, dtype=object)
    at_cols = {W: np.full(n, np.nan) for W in windows}
    code = df["code"].to_numpy(); s = df["s"].to_numpy(); e = df["e"].to_numpy()

    for acc, idx in df.groupby("acc", sort=False).indices.items():
        idx = np.asarray(idx)
        bd = kinds.get(acc, {})
        if len(bd) < 2:
            continue
        bpos = np.array(sorted(bd.keys()), dtype=np.int64)
        ss, ee, cc = s[idx], e[idx], code[idx]
        L = int(max(seqlen.get(acc, 0), ee.max()))

        # ---- nearest-boundary distance (vectorised) ----
        lo = np.searchsorted(bpos, ss, "left")
        hi = np.searchsorted(bpos, ee, "right")
        contains = hi > lo
        below = np.where(lo > 0, ss - bpos[np.clip(lo - 1, 0, len(bpos) - 1)], np.inf)
        above = np.where(hi < len(bpos), bpos[np.clip(hi, 0, len(bpos) - 1)] - ee, np.inf)
        d = np.minimum(below, above); d[contains] = 0.0
        dist[idx] = d

        # nearest boundary position -> its kind + strand -> biological end
        cand_lo = np.clip(lo - 1, 0, len(bpos) - 1)
        cand_hi = np.clip(hi, 0, len(bpos) - 1)
        near_pos = np.where(contains, bpos[np.clip(lo, 0, len(bpos) - 1)],
                            np.where(above < below, bpos[cand_hi], bpos[cand_lo]))
        for j, p in zip(idx, near_pos):
            ks = bd.get(int(p), set())
            k = "both" if ks == {"start", "end"} else (next(iter(ks)) if ks else "")
            st = strands.get(acc, {}).get(int(p), "")
            if k == "start" and st == "+":   biol[j] = "5prime"
            elif k == "end" and st == "+":   biol[j] = "3prime"
            elif k == "start" and st == "-": biol[j] = "3prime"
            elif k == "end" and st == "-":   biol[j] = "5prime"
            else:                            biol[j] = "ambiguous"

        # ---- local A+T fraction in original coordinates, excluding the run itself ----
        at_mask = (cc == 0) | (cc == 3)                 # A or T
        diff = np.zeros(L + 2)
        np.add.at(diff, ss[at_mask], 1.0)
        np.add.at(diff, ee[at_mask] + 1, -1.0)
        at = np.cumsum(diff)[:L + 1]                    # at[p]=1 if coord p is A/T
        csum = np.concatenate([[0.0], np.cumsum(at[1:L + 1])])
        rls = ee - ss + 1
        run_at = np.where(at_mask, rls, 0)
        for W in windows:
            a1 = np.clip(ss - W, 1, L); b1 = np.clip(ee + W, 1, L)
            win_at = csum[b1] - csum[a1 - 1]
            denom = (b1 - a1 + 1) - rls
            at_cols[W][idx] = np.where(denom > 0, (win_at - run_at) / denom, np.nan)

    df["dist"] = dist
    df["biol"] = biol
    for W in windows:
        df[f"localAT_{W}"] = at_cols[W]
    return df


# ------------------------------------------------------------------ statistics
def crude_or(ge, near):
    a1 = int((ge & near).sum()); a0 = int((~ge & near).sum())
    c1 = int((ge & ~near).sum()); c0 = int((~ge & ~near).sum())
    if a0 == 0 or c1 == 0:
        return np.nan, a1, a0, c1, c0
    return (a1 * c0) / (a0 * c1), a1, a0, c1, c0


def mh_or(ge, near, strat):
    """Mantel-Haenszel odds ratio pooled over strata."""
    num = den = 0.0
    ok = np.isfinite(strat) if strat.dtype.kind == "f" else np.ones(len(strat), bool)
    try:
        q = pd.qcut(pd.Series(strat[ok]), NDEC, duplicates="drop", labels=False).to_numpy()
    except Exception:
        return np.nan
    g = ge[ok]; nr = near[ok]
    for st in np.unique(q):
        m = q == st
        T = int(m.sum())
        if T == 0:
            continue
        a1 = int((g[m] & nr[m]).sum());  a0 = int((~g[m] & nr[m]).sum())
        c1 = int((g[m] & ~nr[m]).sum()); c0 = int((~g[m] & ~nr[m]).sum())
        num += a1 * c0 / T
        den += a0 * c1 / T
    return num / den if den > 0 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-table", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--k", type=int, default=5, help="long-run threshold (lambda >= k)")
    ap.add_argument("--near", type=int, default=25, help="near-boundary distance (bp)")
    ap.add_argument("--windows", default="15,50",
                    help="comma-separated +/- windows (bp) for local composition")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    windows = [int(x) for x in args.windows.split(",")]

    print("loading run table + reconstructing boundaries ...", flush=True)
    df, kinds, strands, seqlen = load(args.run_table)
    print(f"  {len(df):,} runs, {df['acc'].nunique()} genomes", flush=True)
    print("annotating distance / biological end / local composition ...", flush=True)
    df = annotate(df, kinds, strands, seqlen, windows)
    df = df[np.isfinite(df["dist"])].copy()

    ge = (df["rl"] >= args.k).to_numpy()
    near = (df["dist"] < args.near).to_numpy()
    nonoverlap = (df["dist"] > 0).to_numpy()
    codes = df["code"].to_numpy()

    rows = []
    subsets = [("all", np.ones(len(df), bool))] + [(b, codes == CODE[b]) for b in BASES]
    for name, sel in subsets:
        cr, a1, a0, c1, c0 = crude_or(ge[sel], near[sel])
        row = dict(subset=name, crude_OR=cr,
                   n_near=int(near[sel].sum()), n_far=int((~near[sel]).sum()),
                   E_near=float(ge[sel & near].mean()) if (sel & near).any() else np.nan,
                   E_far=float(ge[sel & ~near].mean()) if (sel & ~near).any() else np.nan)
        for W in windows:
            strat = df[f"localAT_{W}"].to_numpy()
            row[f"adj_OR_W{W}"] = mh_or(ge[sel], near[sel], strat[sel])
            m = sel & nonoverlap
            row[f"adj_OR_W{W}_no_overlap"] = mh_or(ge[m], near[m], strat[m])
        rows.append(row)
    tab = pd.DataFrame(rows)

    # strand-corrected biological 5' vs 3' ends (A/T runs, the ones carrying the signal)
    at_sel = (codes == CODE["A"]) | (codes == CODE["T"])
    for end in ("5prime", "3prime"):
        m = at_sel & (df["biol"].to_numpy() == end)
        # near/far within runs whose NEAREST boundary is of this biological kind
        cr, *_ = crude_or(ge[m], near[m])
        row = dict(subset=f"A/T @ biological {end}", crude_OR=cr,
                   n_near=int(near[m].sum()), n_far=int((~near[m]).sum()),
                   E_near=float(ge[m & near].mean()) if (m & near).any() else np.nan,
                   E_far=float(ge[m & ~near].mean()) if (m & ~near).any() else np.nan)
        for W in windows:
            strat = df[f"localAT_{W}"].to_numpy()
            row[f"adj_OR_W{W}"] = mh_or(ge[m], near[m], strat[m])
            mm = m & nonoverlap
            row[f"adj_OR_W{W}_no_overlap"] = mh_or(ge[mm], near[mm], strat[mm])
        rows.append(row)
    tab = pd.DataFrame(rows)

    tab.to_csv(args.outdir / "table1_boundary_enrichment.csv", index=False)
    (args.outdir / "table1_boundary_enrichment.md").write_text(
        tab.to_markdown(index=False, floatfmt=".3f"))

    # the compositional shift itself (near vs far), estimated from SHORT runs only
    short = df["rl"] < args.k
    comp_rows = []
    for W in windows:
        v = df[f"localAT_{W}"].to_numpy()
        comp_rows.append(dict(window=W,
                              localAT_near=float(np.nanmean(v[short & near])),
                              localAT_far=float(np.nanmean(v[short & ~near]))))
    comp = pd.DataFrame(comp_rows)
    comp.to_csv(args.outdir / "local_composition_near_far.csv", index=False)

    print("\n=== Table 1: boundary enrichment, crude vs composition-adjusted ===")
    print(tab.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print("\n=== Local A+T fraction (from runs < k, to avoid circularity) ===")
    print(comp.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print(f"\nWrote {args.outdir}/table1_boundary_enrichment.{{csv,md}} and "
          f"local_composition_near_far.csv")


if __name__ == "__main__":
    main()
