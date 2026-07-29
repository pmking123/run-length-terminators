# Annotation-free run-length spectra detect rhabdovirus gene-end signals missed by intrinsic-terminator predictors

Code and data for the paper *"Annotation-free run-length spectra detect rhabdovirus
gene-end signals missed by intrinsic-terminator predictors."*

**Headline result.** An annotation-free homopolymer run-length caller flags
rhabdovirus gene-end (polyadenylation) signals at up to 17.7× the matched-null rate
(95% CI 10.9–23.1, genus-clustered bootstrap; monotonic across thresholds,
2.8× → 5.9× → 17.7× at k ≥ 5, 6, 7), while two structurally distinct terminator
predictors — TransTermHP (pattern-based) and RNIE (covariance-model-based) — detect
3 and 0 such features respectively across the same 36 *Mononegavirales* genomes,
despite both finding hundreds of terminators on the full 200-genome set. These
hairpin-based intrinsic-terminator predictors are blind to this signal class because
it lacks the stem they require — shown directly by folding energy (gene-end tracts
fold to a median −2.65 kcal/mol versus −11.4 for reference-matched terminators) and by
a maximum-sensitivity RNIE re-run that recovers no callable hit at any score.

---

## Reproducing the paper

Every number, table and figure in the paper is regenerated from the frozen inputs by:

```bash
export EXPTERM=/path/to/expterm.dat            # ships with TransTermHP
export RNIE_CM=/path/to/rnie_genomic_1p1.cm    # see "RNIE models" below
bash run_all.sh analysis_final
```

This is deterministic and runs entirely offline. Expect ~1–2 h single-core (the
`--n-null` resampling and ViennaRNA folding dominate).

The confidence intervals, the RNIE threshold-independence check, and the folding-energy
readout are produced by a small set of post-processing steps that consume `run_all.sh`
outputs — see **[Reproducing the supplementary analyses](#reproducing-the-supplementary-analyses)** below.

### Frozen inputs

The analysis is pinned to a fixed set of records, **not** to a live NCBI query:

| file | contents |
|---|---|
| `data/gb_cache/` | 200 GenBank records — the analysed set |
| `data/gb_cache_mono/` | 36 *Mononegavirales* records (a subset of the above) |
| `data/accessions_all.txt` | the 200 accessions |
| `data/accessions_mono.txt` | the 36 *Mononegavirales* accessions |
| `data/features.csv` | per-genome feature table driving the pipeline |

These were retrieved from NCBI nuccore on 11 June 2026 with the Entrez query

```
txid10239[Organism] AND complete genome[Title] AND RefSeq[filter]
```

taking the **first 200 records under NCBI's default ordering** (`esearch`/`efetch`,
`db=nuccore`, no `sort` specified) with **no dereplication**. The largest genome is
466,767 bp. **Re-running that query today will return a different set** — RefSeq grows,
`[Title]` is a free-text match, and default-order retrieval is not stable over time —
which is precisely why the GenBank cache is archived here as the reproducible reference
rather than a live query. Nothing in the pipeline contacts the network.

Family composition (25 families) is dominated by bacteriophages (Drexlerviridae 47,
Steitzviridae 25), with a 36-genome *Mononegavirales* bloc (Rhabdoviridae 29,
Paramyxoviridae 5, Pneumoviridae 2); the full breakdown is Supplementary Table S1 in
the paper. Because RefSeq contains many near-identical strains, the analyses handle the
resulting phylogenetic non-independence explicitly (per-genome, genus-clustered, and
clade-aware) rather than by pre-filtering — see **Statistical design notes**.

### Dependencies

Python (see `environment.yml`): pandas, numpy, scipy, biopython, matplotlib, ViennaRNA.

External tools:

* **TransTermHP** v2.09 — provides `transterm` and `expterm.dat`.
* **Infernal** 1.1.x — provides `cmsearch` and `cmconvert`.
* **ViennaRNA** — install the *Python bindings* (`pip install ViennaRNA`). Note the
  scripts fall back to the `RNAfold` binary, and then to a crude inverted-repeat
  heuristic, if the bindings are unavailable; only the bindings/binary paths give
  meaningful hairpin statistics.

### RNIE models

RNIE (Gardner et al. 2011) was built against Infernal 1.0 and does not run against
Infernal 1.1. We use its trained covariance models directly:

```bash
cmconvert RNIE/models/genomic.cm > rnie_genomic_1p1.cm
cmsearch --tblout out.tbl -T 14 rnie_genomic_1p1.cm genomes.fasta
```

Because the converted models are not calibrated for 1.1, we threshold on **bit score**
(`-T 14`, RNIE's recommended genomic-mode cut), not E-value. The paper additionally
reports a maximum-sensitivity re-run that removes this threshold entirely to show the
blindness is threshold-independent — see the supplementary analyses below.

---

## Reproducing the supplementary analyses

These three steps run **after** `bash run_all.sh analysis_final` and consume its
outputs. Each writes a small CSV and prints the numbers used in the paper. They are
separated from `run_all.sh` because two of them (the bootstrap, the RNIE re-run) are
post-hoc robustness analyses rather than part of the primary pipeline.

### 1. Genus-clustered bootstrap CIs (Table 3, abstract)

`clade_enrichment_table.py` (run inside `run_all.sh`) emits
`analysis_final/clade_enrichment/per_candidate_by_k.csv`, one row per candidate tract
with its genus. The bootstrap resamples genera with replacement to put confidence
intervals on the per-clade enrichment, and reports leave-one-genus-out to show no
single genus drives the effect:

```bash
python scripts/bootstrap_enrichment.py \
    --per-candidate analysis_final/clade_enrichment/per_candidate_by_k.csv \
    --family Rhabdoviridae
# -> analysis_final/clade_enrichment/Rhabdoviridae_enrichment_ci.csv
#    (rows for k = 5, 6, 7; enrichment, 95% CI, leave-one-genus-out band, n_genera)
```

The null is held fixed at its published per-k value; only the candidate tracts are
resampled, clustered by genus. This answers "is the enrichment driven by dense sampling
within one genus?" without needing per-genome null counts (which the benchmark summary
does not expose).

### 2. RNIE threshold-independence (maximum-sensitivity re-run + sweep)

The frozen benchmark applies RNIE's recommended `-T 14` cut. To show the zero count is
not an artefact of that threshold, re-run `cmsearch` at maximum sensitivity — all
heuristic filters off, no score threshold — then walk the bit-score cutoff down over
the resulting hits:

```bash
cmsearch --tblout analysis_final/rnie_mono_permissive.tbl --max --cpu 8 \
    "$RNIE_CM" analysis_final/transterm_inputs_mono/genomes.fasta

python scripts/rnie_score_sweep.py \
    --input analysis_final/rnie_mono_permissive.tbl \
    --out analysis_final/bench_mono_rnie/rnie_score_sweep.csv
```

The re-run returns only four hits across the 36 genomes, all below `-T 14` (scores
11.4–13.1) and all structurally degenerate poly-A/U matches rather than hairpins; the
count is flat from score −10 up to 14, so no cutoff recovers a callable terminator.
`rnie_score_sweep.py` reuses `predictor_to_bed.parse_cmsearch`, so its counts are
directly comparable to the frozen benchmark.

### 3. Folding-energy readout (Results 3.2)

The "hairpin-poor" claim is shown directly from the per-candidate `hairpin_mfe` /
`has_hairpin` columns, contrasting rhabdovirus gene-end tracts against the tracts a
hairpin tool actually recognised (`matches_reference`), with no taxonomy required:

```bash
python scripts/mfe_readout.py \
    --full-candidates analysis_final/bench_all_transterm/terminator_candidates.csv \
    --mono-candidates analysis_final/bench_mono_k6/terminator_candidates.csv \
    --family-map analysis_final/clade_enrichment/genome_family_assignment.csv \
    --out analysis_final/bench_mono_rnie/mfe_readout.csv
```

Reference-matched terminators fold to a median −11.4 kcal/mol (98% pass the hairpin
call); rhabdovirus gene-end tracts fold to −2.65 (46% clear the permissive −3 screen,
but only 6.8% reach −8 and none reach −11), i.e. marginal structure well short of a
terminator stem.

---

## Layout

```
scripts/                        analysis code (see below)
data/                           frozen inputs (GenBank caches, accession lists)
analysis_final/                 all derived outputs (regenerated by run_all.sh)
run_all.sh                      one-command reproduction
environment.yml                 conda environment
```

This repository contains the code and frozen data only; the manuscript is not included
here (see **Citation**).

### Scripts

| script | role |
|---|---|
| `ccdna.py` | core run-length decomposition and per-genome features |
| `boundary_enrichment.py` | builds `run_table.csv`: every run + distance to nearest CDS boundary |
| `distance_profile.py` | E_k(d) profiles by distance bin |
| `random_boundary_control.py` | uniform/rotation boundary nulls |
| `boundary_composition_controls.py` | **Table 1** — crude vs Mantel–Haenszel composition-adjusted ORs; strand-corrected 5′/3′ split |
| `per_genome_paired_test.py` | pseudoreplication control — genome as unit, paired Wilcoxon vs matched nulls |
| `threeprime_motifs.py` | **Figure 1** — 3′ positional profiles and gene-end motif enrichment |
| `gb_to_transterm.py` | GenBank cache → TransTermHP FASTA + `.crd` (IDs guaranteed to match) |
| `predictor_to_bed.py` | TransTermHP / RNIE / cmsearch / BED → common reference format |
| `subset_gb_by_taxon.py` | taxonomy-based cache subsetting (e.g. *Mononegavirales*) |
| `terminator_benchmark.py` | **Table 2** — run-length caller vs reference predictors, with matched null |
| `clade_enrichment_table.py` | **Figure 2a, Table 3** — per-clade gene-end enrichment across thresholds; also emits `per_candidate_by_k.csv` (genus-tagged, for the bootstrap) |
| `bootstrap_enrichment.py` | **Table 3 CIs** — genus-clustered bootstrap + leave-one-genus-out on the per-clade enrichment |
| `rnie_score_sweep.py` | **RNIE threshold-independence** — re-thresholds a cmsearch table across bit-score cutoffs |
| `mfe_readout.py` | **Results 3.2** — folding-energy contrast, rhabdovirus gene-end tracts vs reference-matched terminators |

---

## Statistical design notes

Three controls do the load-bearing work, and are worth understanding before
reinterpreting any output:

1. **Pseudoreplication.** Runs within a genome are not independent. Pooling ~3.2 M runs
   yields wildly anti-conservative p-values. All inference in the paper uses the
   **genome** as the unit (`per_genome_paired_test.py`): per-genome rates compared
   against matched nulls by paired Wilcoxon across 199 genomes (one genome, NC_139179.1,
   has no run of length ≥ 5 and so no defined long-run rate).

2. **Composition matching.** Boundaries sit in locally A+T-rich sequence, and run-length
   tails scale steeply with base frequency. `boundary_composition_controls.py` therefore
   reports Mantel–Haenszel odds ratios stratified by local A+T decile, at two window
   sizes, with and without boundary-overlapping runs.

3. **Matched nulls in the benchmark.** A "candidate the predictor missed" means nothing
   without knowing how often the corroboration criteria fire by chance.
   `terminator_benchmark.py --null-per-genome` evaluates the identical hairpin and
   gene-end tests at random positions. Report the *enrichment over that null*, never the
   raw corroboration count.

**Phylogenetic non-independence.** RefSeq contains dense clusters of near-identical
strains, so naive pooling or random cross-validation is severely inflated by
phylogenetic leakage (close relatives fall in both train and test). The per-clade
enrichment is therefore reported with genus-clustered bootstrap CIs and
leave-one-genus-out (`bootstrap_enrichment.py`), and genus is assigned from the GenBank
lineage rank-aware (last single-token `-virus` element), not by organism-name matching.

**The strand model matters.** `--gene-end-mode same` encodes intrinsic-terminator
geometry (U-tract on the gene strand); `polyA` encodes *Mononegavirales* gene-end
geometry (A-tract 3′ of a same-sense gene). Applying `same` to *Mononegavirales*
suppresses the signal ~17-fold — an artifact, not a negative result. The paper reports
`any` (the conservative choice); `polyA` gives a higher enrichment.

---

## Citation

If you use this code, please cite the paper and the archived dataset. The manuscript is
maintained separately and is not included in this repository; its full reference will be
added here on publication. The dataset is archived at Zenodo under concept DOI
10.5281/zenodo.21311759 (which always resolves to the latest version); `CITATION.cff`
carries the machine-readable citation metadata.

## Licence

Code: MIT (see `LICENSE`). GenBank records in `data/` are public-domain NCBI records,
redistributed here to pin the analysed set.
