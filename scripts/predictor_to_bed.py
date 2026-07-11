#!/usr/bin/env python3
"""
predictor_to_bed.py

Convert a terminator predictor's output into the 4-column reference file that
terminator_benchmark.py reads via --reference-bed:

    accession  start  end  strand

Supports:
  --format rnie   RNIE GFF output (default): seqid, source, feature, start, end,
                  score, strand, ... ; interval = cols 4/5, strand = col 7.
  --format gff    generic GFF3 (same columns as rnie).
  --format bed    BED6 (chrom, start, end, name, score, strand); e.g. BATTER.
  --format cmsearch   Infernal 1.1 cmsearch --tblout table (whitespace-delimited):
                  target(0) ... seq_from(7) seq_to(8) strand(9) ... score(14) Evalue(15).

Notes
-----
* seqid/chrom must equal the GenBank record.id the benchmark keys on. If you ran
  the predictor on the genomes.fasta produced by gb_to_transterm.py (headers are
  bare accessions), this is automatic.
* --min-score filters by the score column (GFF col 6 / BED col 5). RNIE genome
  mode already applies an internal bit-score cut; use --min-score only if you
  want a stricter set.
* BED is 0-based half-open; we emit 1-based inclusive to match GFF/GenBank. The
  benchmark takes min/max of start/end, so small off-by-ones don't affect overlap
  within the default tolerance, but the conversion is done correctly regardless.

Example
-------
python predictor_to_bed.py --input rnie_mono.gff --format rnie \
    --output reference_rnie_mono.bed
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd


def parse_gff(path, min_score):
    rows = []
    for line in Path(path).read_text(errors="ignore").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 7:
            f = line.split()
            if len(f) < 7:
                continue
        seqid, start, end, score, strand = f[0], f[3], f[4], f[5], f[6]
        try:
            s, e = int(start), int(end)
        except ValueError:
            continue
        if min_score is not None:
            try:
                if float(score) < min_score:
                    continue
            except ValueError:
                pass
        if strand not in ("+", "-"):
            strand = "+"
        rows.append((seqid, min(s, e), max(s, e), strand))
    return rows


def parse_bed(path, min_score):
    rows = []
    for line in Path(path).read_text(errors="ignore").splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split()
        if len(f) < 3:
            continue
        chrom = f[0]
        s0, e0 = int(f[1]), int(f[2])            # 0-based half-open
        score = f[4] if len(f) > 4 else "."
        strand = f[5] if len(f) > 5 else "+"
        if min_score is not None:
            try:
                if float(score) < min_score:
                    continue
            except ValueError:
                pass
        if strand not in ("+", "-"):
            strand = "+"
        rows.append((chrom, s0 + 1, e0, strand))   # -> 1-based inclusive
    return rows


def parse_cmsearch(path, min_score):
    """Infernal 1.1 cmsearch --tblout: whitespace-delimited, fixed columns.
    target(0) accession(1) query(2) accession(3) mdl(4) mdl_from(5) mdl_to(6)
    seq_from(7) seq_to(8) strand(9) trunc(10) pass(11) gc(12) bias(13)
    score(14) E-value(15) inc(16) description(17..). We use score as the filter
    (bit score), since RNIE models are best thresholded by score, not E-value."""
    rows = []
    for line in Path(path).read_text(errors="ignore").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 16:
            continue
        seqid = f[0]
        try:
            s, e = int(f[7]), int(f[8]); score = float(f[14])
        except ValueError:
            continue
        strand = f[9] if f[9] in ("+", "-") else "+"
        if min_score is not None and score < min_score:
            continue
        rows.append((seqid, min(s, e), max(s, e), strand))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--format", choices=["rnie", "gff", "bed", "cmsearch"], default="rnie")
    ap.add_argument("--min-score", type=float, default=None)
    args = ap.parse_args()

    if args.format in ("rnie", "gff"):
        rows = parse_gff(args.input, args.min_score)
    elif args.format == "cmsearch":
        rows = parse_cmsearch(args.input, args.min_score)
    else:
        rows = parse_bed(args.input, args.min_score)

    df = pd.DataFrame(rows, columns=["accession", "start", "end", "strand"])
    df.to_csv(args.output, sep="\t", header=False, index=False)
    print(f"parsed {len(df)} terminator intervals from {args.input} ({args.format})")
    if len(df):
        print(f"  {df['accession'].nunique()} distinct accessions")
        print(f"  wrote {args.output}  (feed to terminator_benchmark.py --reference-bed)")
    else:
        print("  WARNING: no intervals parsed - check --format and that the file has content")


if __name__ == "__main__":
    main()
