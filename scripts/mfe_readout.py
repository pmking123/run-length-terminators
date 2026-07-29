#!/usr/bin/env python3
"""mfe_readout.py

Direct structural evidence for the "hairpin-poor" claim in Results 3.2, replacing the
G+C proxy with the folding energy the caller already computed per candidate.

Three groups, same folding pipeline, no taxonomy:

  1. reference-matched hairpin terminators  (matches_reference == True, full 200-set)
     -- tracts a hairpin-based tool (TransTermHP/RNIE) actually recognised: the
        operational definition of "structured", in the benchmark's own terms.
  2. rhabdovirus gene-end tracts            (gene_end_proximal == True, mono set,
        family == Rhabdoviridae via the genus/family join)
     -- the class we claim is hairpin-poor.
  3. matched null                            -- optional, if a null candidate table with
        the same columns is supplied; otherwise reported as N/A.

For each group: n, fraction passing the caller's own hairpin call (has_hairpin, i.e.
MFE <= -3 kcal/mol), and the MFE distribution (median, mean, IQR). The claim is a
DISTRIBUTION SHIFT, not a clean separation -- some rhabdovirus gene-end tracts do fold
(seen in the raw rows) -- so the readout reports rates and quantiles, and the wording it
supports is "hairpin-poor", never "hairpin-free".

Inputs
------
--full-candidates  full 200-set candidate table (e.g. bench_all_transterm/terminator_candidates.csv)
--mono-candidates  Mononegavirales candidate table at the headline k (e.g. bench_mono_k6/...)
--family-map       genome_family_assignment.csv (accession -> family) from clade_enrichment
--null-candidates  OPTIONAL matched-null candidate table with the same columns
--out              CSV summary path

Columns expected in candidate tables:
  accession, tract_len, hairpin_mfe, has_hairpin, hairpin_gc, gene_end_proximal,
  matches_reference, corroborated
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def as_bool(s):
    if s.dtype == bool:
        return s
    return s.astype(str).str.strip().str.lower().isin(("true", "1", "yes"))


def summarise(name, g):
    if len(g) == 0:
        return dict(group=name, n=0, has_hairpin_rate=np.nan, mfe_median=np.nan,
                    mfe_mean=np.nan, mfe_q25=np.nan, mfe_q75=np.nan)
    mfe = pd.to_numeric(g["hairpin_mfe"], errors="coerce")
    return dict(group=name, n=len(g),
                has_hairpin_rate=round(as_bool(g["has_hairpin"]).mean(), 3),
                mfe_median=round(mfe.median(), 2), mfe_mean=round(mfe.mean(), 2),
                mfe_q25=round(mfe.quantile(0.25), 2), mfe_q75=round(mfe.quantile(0.75), 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full-candidates", required=True, type=Path)
    ap.add_argument("--mono-candidates", required=True, type=Path)
    ap.add_argument("--family-map", required=True, type=Path)
    ap.add_argument("--null-candidates", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    full = pd.read_csv(args.full_candidates)
    mono = pd.read_csv(args.mono_candidates)
    fam = pd.read_csv(args.family_map)[["accession", "family"]]
    fam_of = dict(zip(fam["accession"].astype(str), fam["family"]))

    full["matches_reference"] = as_bool(full["matches_reference"])
    full["gene_end_proximal"] = as_bool(full["gene_end_proximal"])
    mono["gene_end_proximal"] = as_bool(mono["gene_end_proximal"])
    mono["family"] = mono["accession"].astype(str).map(fam_of)

    # group 1: reference-matched hairpin terminators (the "structured" comparator)
    ref_hairpin = full[full["matches_reference"]]
    # group 2: rhabdovirus gene-end tracts (the class claimed hairpin-poor)
    rhabdo_ge = mono[(mono["family"] == "Rhabdoviridae") & (mono["gene_end_proximal"])]

    records = [summarise("reference-matched hairpin terminators", ref_hairpin),
               summarise("rhabdovirus gene-end tracts", rhabdo_ge)]

    if args.null_candidates and args.null_candidates.exists():
        nul = pd.read_csv(args.null_candidates)
        records.append(summarise("matched null", nul))
    else:
        records.append(dict(group="matched null", n=0, has_hairpin_rate=np.nan,
                            mfe_median=np.nan, mfe_mean=np.nan, mfe_q25=np.nan, mfe_q75=np.nan))

    out = pd.DataFrame(records)
    pd.set_option("display.width", 140)
    print(out.to_string(index=False))

    # the one-line contrast the paper needs
    r = {d["group"]: d for d in records}
    ref = r["reference-matched hairpin terminators"]
    rh = r["rhabdovirus gene-end tracts"]
    if ref["n"] and rh["n"]:
        print(f"\nCONTRAST: reference-matched hairpin terminators pass the hairpin call at "
              f"{ref['has_hairpin_rate']:.0%} (median MFE {ref['mfe_median']} kcal/mol); "
              f"rhabdovirus gene-end tracts at {rh['has_hairpin_rate']:.0%} "
              f"(median MFE {rh['mfe_median']} kcal/mol).")
        print("NB: distribution shift, not clean separation -- some rhabdovirus tracts do "
              "fold. Wording: 'hairpin-poor', NOT 'hairpin-free'.")

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.out, index=False)
        print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
