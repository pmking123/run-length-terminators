#!/usr/bin/env python3
"""rnie_score_sweep.py

Re-threshold an existing Infernal cmsearch --tblout across a range of bit-score
cutoffs, and report how many hits survive at each and across how many genomes.

Reuses predictor_to_bed.parse_cmsearch so the counts are apples-to-apples with the
frozen benchmark (same column indexing, same strand handling) -- the only thing that
changes across the sweep is the threshold, which is the variable under test.

Intended use: after a permissive re-run (cmsearch --max, no -T), to establish whether
RNIE registers these genomes at all below its recommended -T 14 cutoff, and if so, at
what score. Prints the raw score distribution so you can see how far below 14 anything
sits, and (optionally) whether surviving hits fall within a gene-end window.

Usage
-----
  python rnie_score_sweep.py --input analysis_final/rnie_mono_permissive.tbl \
      --out analysis_final/bench_mono_rnie/rnie_score_sweep.csv

  # optional gene-end overlap check: pass a 4-col reference-style BED of CDS 3' ends
  # (accession, start, end, strand) and a window; reports how many swept hits land within.
  #   --gene-ends analysis_final/.../gene_end_windows.bed --window 50
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import pandas as pd

# reuse the exact parser the frozen benchmark used
sys.path.insert(0, str(Path(__file__).resolve().parent))
from predictor_to_bed import parse_cmsearch


def raw_scores(path: Path):
    """Read the bit-score column (14) directly, no threshold, for the distribution."""
    scores = []
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) >= 16:
            try:
                scores.append(float(f[14]))
            except ValueError:
                pass
    return pd.Series(scores, dtype=float)


def load_gene_end_windows(path: Path):
    """4-col reference BED: accession, start, end, strand (1-based inclusive)."""
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["accession", "start", "end", "strand"])
    return df


def overlaps_gene_end(acc, s, e, windows, pad):
    w = windows[windows.accession == acc]
    if w.empty:
        return False
    lo, hi = min(s, e) - pad, max(s, e) + pad
    return bool(((w.start <= hi) & (w.end >= lo)).any())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--cutoffs", default="14,12,10,8,6,4,2,0,-5,-10")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--gene-ends", type=Path, default=None,
                    help="optional 4-col BED of CDS 3' ends for an overlap check")
    ap.add_argument("--window", type=int, default=50)
    args = ap.parse_args()

    cutoffs = [float(x) for x in args.cutoffs.split(",")]
    windows = load_gene_end_windows(args.gene_ends) if args.gene_ends else None

    rows = []
    for c in cutoffs:
        r = parse_cmsearch(args.input, min_score=c)          # [(acc, start, end, strand), ...]
        accs = {x[0] for x in r}
        ge = None
        if windows is not None:
            ge = sum(overlaps_gene_end(a, s, e, windows, args.window) for (a, s, e, _) in r)
        rows.append(dict(min_score=c, n_hits=len(r), n_genomes=len(accs),
                         n_gene_end=("" if ge is None else ge)))
        msg = f"  min_score >= {c:>6}:  {len(r):4d} hits across {len(accs):2d} genomes"
        if ge is not None:
            msg += f";  {ge} within {args.window} bp of a gene-end"
        print(msg)

    s = raw_scores(args.input)
    print()
    if len(s):
        print(f"raw hits: {len(s)} total; score range [{s.min():.1f}, {s.max():.1f}], "
              f"median {s.median():.1f}")
        print(f"  hits with score >= 14 (RNIE's recommended -T): {(s >= 14).sum()}")
        print(f"  hits with score in [0, 14):                    {((s >= 0) & (s < 14)).sum()}")
    else:
        print("No raw hits in the table at any score.")

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(args.out, index=False)
        print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
