#!/usr/bin/env python3
"""
gb_to_transterm.py

Convert a GenBank cache into the two inputs TransTermHP needs, with FASTA headers
and .crd chrom_ids guaranteed to equal the GenBank record.id -- which is the same
id terminator_benchmark.py keys on. This removes the silent ID-mismatch pitfall.

Emits, by default, one combined FASTA and one combined .crd so you can run:

    transterm -p expterm.dat genomes.fasta genomes.crd > transterm.tt

.crd line format is:  gene_name  start  end  chrom_id
Strand is encoded the TransTermHP way: forward genes have start < end, reverse
genes have start > end (coordinates swapped). TransTermHP needs >= 2 genes per
sequence; records with fewer real CDS get two fake 1-bp flanking genes (the
annotation-free trick from the TransTermHP manual, section 10) so they are still
processed -- these records are listed in the manifest so you can treat them
separately if you wish.

Outputs (in --outdir):
    genomes.fasta
    genomes.crd
    gb_to_transterm_manifest.csv   (accession, length, n_cds, used_fake_genes)

Needs: biopython, pandas.
"""
from __future__ import annotations
import argparse, re
from pathlib import Path
import pandas as pd
from Bio import SeqIO


def safe_token(text, fallback):
    """A single whitespace-free token for the gene_name column."""
    t = re.sub(r"\s+", "_", str(text)).strip("_")
    t = re.sub(r"[^A-Za-z0-9_.:-]", "", t)
    return t or fallback


def gene_name(feat, idx):
    for key in ("locus_tag", "gene", "protein_id", "product"):
        if key in feat.qualifiers and feat.qualifiers[key]:
            return safe_token(feat.qualifiers[key][0], f"CDS_{idx}")
    return f"CDS_{idx}"


def wrap(seq, width):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gb-cache", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--fasta-name", default="genomes.fasta")
    ap.add_argument("--crd-name", default="genomes.crd")
    ap.add_argument("--wrap", type=int, default=70)
    ap.add_argument("--feature-types", default="CDS",
                    help="comma-separated feature types to emit as genes")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    ftypes = set(t.strip() for t in args.feature_types.split(","))

    gbs = sorted(p for p in args.gb_cache.glob("*.gb") if not p.name.endswith(".tmp.gb"))
    if not gbs:
        raise SystemExit(f"No .gb files in {args.gb_cache}")

    fasta_lines, crd_lines, manifest = [], [], []
    seen_ids = set()
    for gb in gbs:
        for record in SeqIO.parse(str(gb), "genbank"):
            acc = record.id
            if acc in seen_ids:
                print(f"WARNING: duplicate record id {acc}; skipping later copy")
                continue
            seen_ids.add(acc)
            seq = str(record.seq).upper()
            L = len(seq)
            if L == 0:
                continue
            fasta_lines.append(f">{acc}")
            fasta_lines.append(wrap(seq, args.wrap))

            genes = []
            idx = 0
            for feat in record.features:
                if feat.type not in ftypes:
                    continue
                strand = feat.location.strand
                if strand not in (1, -1):
                    continue
                try:
                    l = int(feat.location.start) + 1   # 1-based left
                    r = int(feat.location.end)          # 1-based right
                except Exception:
                    continue
                if r <= l:
                    continue
                idx += 1
                nm = gene_name(feat, idx)
                # forward: start<end ; reverse: start>end (swap)
                if strand == 1:
                    genes.append((nm, l, r))
                else:
                    genes.append((nm, r, l))

            used_fake = False
            if len(genes) < 2:
                # TransTermHP needs >=2 genes: add fake tail-to-tail flankers
                genes = [("fakegene1", 1, 2), ("fakegene2", L - 1, L)]
                used_fake = True
            for nm, s, e in genes:
                crd_lines.append(f"{nm}\t{s}\t{e}\t{acc}")

            manifest.append(dict(accession=acc, length=L,
                                 n_cds=idx, used_fake_genes=used_fake))

    (args.outdir / args.fasta_name).write_text("\n".join(fasta_lines) + "\n")
    (args.outdir / args.crd_name).write_text("\n".join(crd_lines) + "\n")
    man = pd.DataFrame(manifest)
    man.to_csv(args.outdir / "gb_to_transterm_manifest.csv", index=False)

    n_fake = int(man["used_fake_genes"].sum()) if len(man) else 0
    print(f"records: {len(man)}   genes written: {len(crd_lines)}")
    print(f"records using fake flanking genes (<2 CDS): {n_fake}")
    print(f"\nWrote:\n  {args.outdir/args.fasta_name}\n  {args.outdir/args.crd_name}"
          f"\n  {args.outdir/'gb_to_transterm_manifest.csv'}")
    print("\nRun TransTermHP with:")
    print(f"  transterm -p expterm.dat {args.outdir/args.fasta_name} "
          f"{args.outdir/args.crd_name} > transterm.tt")


if __name__ == "__main__":
    main()
