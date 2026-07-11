#!/usr/bin/env python3
"""
clade_enrichment_table.py

Per-clade gene-end tract enrichment across tract-length thresholds - the data behind
Figure 2a and the k-sweep table in Results 3.4.

Family assignment is taken from the GenBank TAXONOMY lineage (authoritative), NOT from
organism-name pattern matching. Many rhabdoviruses are named for host or locality
("Shanxi Arboretum virus", "Daphne virus 1") and carry no taxonomic signal in the
organism string, so name-based rules silently misassign them.

For each terminator_benchmark output directory (one per k), this reads the candidate
table and that run's own matched-null rate, assigns each candidate's genome to a
family from its GenBank record, and reports per family:

    n           total run-length candidates
    n_gene_end  candidates within the gene-end window
    frac        n_gene_end / n
    enrichment  frac / (that run's matched-null gene-end rate)

Using each run's own null keeps the ratio self-consistent (the null depends on k).

Any genome whose lineage contains no family is reported as UNASSIGNED with a loud
warning; it is never silently pooled.

Inputs
------
--bench-dirs  comma-separated terminator_benchmark output dirs (one per k)
--gb-cache    the GenBank cache the benchmark was run on (for taxonomy lineages)

Outputs
-------
<outdir>/clade_enrichment_by_k.csv     long: family, k, n, n_gene_end, frac, enrichment
<outdir>/clade_enrichment_by_k.md      pivot table, paste-ready
<outdir>/genome_family_assignment.csv  accession -> organism, family, lineage (audit)

Needs: pandas, numpy, biopython.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO


def family_from_record(record):
    """Return (family, lineage_string): the lineage element ending in 'viridae'."""
    lineage = record.annotations.get("taxonomy", [])
    lin_str = "; ".join(lineage)
    for tax in lineage:
        t = tax.strip()
        if t.endswith("viridae"):
            return t, lin_str
    return "UNASSIGNED", lin_str


def load_taxonomy(gb_cache: Path) -> pd.DataFrame:
    rows = []
    for gb in sorted(p for p in gb_cache.glob("*.gb") if not p.name.endswith(".tmp.gb")):
        for record in SeqIO.parse(str(gb), "genbank"):
            fam, lin = family_from_record(record)
            rows.append(dict(accession=record.id,
                             organism=record.annotations.get("organism", ""),
                             family=fam, lineage=lin))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench-dirs", required=True)
    ap.add_argument("--gb-cache", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    tax = load_taxonomy(args.gb_cache)
    if len(tax) == 0:
        raise SystemExit(
            f"ERROR: no GenBank records parsed from {args.gb_cache}\n"
            f"  Expected a DIRECTORY containing *.gb files, e.g. data/gb_cache_mono\n"
            f"  Check: ls {args.gb_cache}/*.gb | head")
    tax.to_csv(args.outdir / "genome_family_assignment.csv", index=False)
    fam_of = dict(zip(tax["accession"].astype(str), tax["family"]))
    print(f"taxonomy: {len(tax)} genomes")
    print(tax["family"].value_counts().to_string())
    un = tax[tax["family"] == "UNASSIGNED"]
    if len(un):
        print("\n*** WARNING: genomes with no family in lineage ***")
        print(un[["accession", "organism"]].to_string(index=False))

    rows = []
    for d in [Path(x.strip()) for x in args.bench_dirs.split(",") if x.strip()]:
        cf, sf = d / "terminator_candidates.csv", d / "terminator_benchmark_summary.csv"
        if not cf.exists() or not sf.exists():
            print(f"WARNING: skipping {d}")
            continue
        C = pd.read_csv(cf); S = pd.read_csv(sf)
        null = float(S["null_gene_end_proximal_rate"].iloc[0])
        k = int(C["tract_len"].min()) if len(C) else np.nan
        C["gene_end"] = C["gene_end_proximal"].astype(bool)
        C["family"] = C["accession"].astype(str).map(fam_of).fillna("UNASSIGNED")
        for fam, sub in C.groupby("family"):
            n = len(sub); ge = int(sub["gene_end"].sum())
            frac = ge / n if n else np.nan
            rows.append(dict(family=fam, k=k, n=n, n_gene_end=ge, frac=frac,
                             null_rate=null, enrichment=(frac / null) if null else np.nan))

    out = pd.DataFrame(rows).sort_values(["family", "k"])
    out.to_csv(args.outdir / "clade_enrichment_by_k.csv", index=False)

    piv = out.pivot(index="family", columns="k", values="enrichment")
    cnt = out.pivot(index="family", columns="k", values="n_gene_end")
    tot = out.pivot(index="family", columns="k", values="n")
    md = ["# Per-clade gene-end enrichment (vs matched null)", "",
          "Enrichment, with (gene-end tracts / total candidates):", "",
          "| family | " + " | ".join(f"k>={int(k)}" for k in piv.columns) + " |",
          "|" + "---|" * (len(piv.columns) + 1)]
    for fam in piv.index:
        cells = []
        for k in piv.columns:
            e, c, t = piv.loc[fam, k], cnt.loc[fam, k], tot.loc[fam, k]
            cells.append("-" if pd.isna(e) else f"{e:.1f}x ({int(c)}/{int(t)})")
        md.append(f"| {fam} | " + " | ".join(cells) + " |")
    (args.outdir / "clade_enrichment_by_k.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    print("\n=== Per-clade gene-end enrichment by threshold ===")
    print(out.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    n_un = int(out.loc[out["family"] == "UNASSIGNED", "n"].sum()) if len(out) else 0
    if n_un:
        print(f"\n*** {n_un} candidates UNASSIGNED - fix before using this table ***")
    print(f"\nWrote {args.outdir}/clade_enrichment_by_k.{{csv,md}}")


if __name__ == "__main__":
    main()
