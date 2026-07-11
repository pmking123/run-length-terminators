#!/usr/bin/env python3
"""
subset_gb_by_taxon.py

Copy the GenBank files whose taxonomy matches any of the given terms into a new
directory, so you can run the terminator benchmark on a specific clade (e.g. the
Mononegavirales, where hairpin-based predictors like TransTermHP are blind).

Matching is case-insensitive against the record's taxonomy lineage AND organism
string, so --taxa Mononegavirales catches Rhabdoviridae, Paramyxoviridae,
Filoviridae, Pneumoviridae, etc. in one go.

Example
-------
python subset_gb_by_taxon.py --gb-cache gb_cache/ --outdir gb_cache_mono/ \
    --taxa Mononegavirales

Outputs: the matching *.gb files copied into --outdir, plus
subset_manifest.csv (accession, organism, matched_term).

Needs: biopython.
"""
from __future__ import annotations
import argparse, shutil
from pathlib import Path
import pandas as pd
from Bio import SeqIO


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gb-cache", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--taxa", required=True,
                    help="comma-separated taxonomy terms to match (any-of)")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    terms = [t.strip().lower() for t in args.taxa.split(",") if t.strip()]

    gbs = sorted(p for p in args.gb_cache.glob("*.gb") if not p.name.endswith(".tmp.gb"))
    if not gbs:
        raise SystemExit(f"No .gb files in {args.gb_cache}")

    rows = []
    copied = 0
    for gb in gbs:
        matched_term = None
        organism = ""
        for record in SeqIO.parse(str(gb), "genbank"):
            lineage = " ; ".join(record.annotations.get("taxonomy", []))
            organism = record.annotations.get("organism", "")
            hay = (lineage + " ; " + organism).lower()
            for t in terms:
                if t in hay:
                    matched_term = t
                    break
            if matched_term:
                break
        if matched_term:
            shutil.copy2(gb, args.outdir / gb.name)
            copied += 1
            rows.append(dict(accession=gb.stem, organism=organism,
                             matched_term=matched_term))

    man = pd.DataFrame(rows)
    man.to_csv(args.outdir / "subset_manifest.csv", index=False)
    print(f"scanned {len(gbs)} genomes; copied {copied} matching {terms}")
    if copied:
        print(man["organism"].value_counts().to_string())
    print(f"\nsubset written to {args.outdir}/  (+ subset_manifest.csv)")


if __name__ == "__main__":
    main()
