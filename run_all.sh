#!/usr/bin/env bash
#
# run_all.sh — regenerate every number in the paper from the frozen GenBank cache.
#
# Usage:
#   bash run_all.sh [OUTDIR]
#
# Prerequisites (see environment.yml / README):
#   * python 3.10+ with pandas, numpy, scipy, biopython, matplotlib, ViennaRNA
#   * transterm (TransTermHP v2.09) on PATH, and its expterm.dat
#   * cmsearch  (Infernal 1.1.x) on PATH, and the converted RNIE model
#
# Inputs (must exist; these are the FROZEN inputs archived with this repo):
#   data/gb_cache/        200 GenBank records (the analysed set)
#   data/gb_cache_mono/   34 Mononegavirales records (subset of the above)
#   data/accessions_all.txt, data/accessions_mono.txt
#
# Everything below is deterministic given those inputs. No network access is used.
#
set -euo pipefail

OUT="${1:-analysis_final}"
GB=data/gb_cache
GB_MONO=data/gb_cache_mono

# --- external tool locations: EDIT THESE or export before running -------------
EXPTERM="${EXPTERM:?set EXPTERM=/path/to/expterm.dat}"
RNIE_CM="${RNIE_CM:?set RNIE_CM=/path/to/rnie_genomic_1p1.cm}"

mkdir -p "$OUT"
echo "== outputs -> $OUT"

# -----------------------------------------------------------------------------
# Stage 1. Core pipeline: run table + boundary enrichment
# -----------------------------------------------------------------------------
echo "== [1/5] boundary enrichment (run table)"
python scripts/boundary_enrichment.py \
  --features data/features.csv \
  --email placeholder@example.org \
  --outdir "$OUT/boundary_enrichment" \
  --cache-dir "$GB" \
  --boundary-window 250 \
  --resume
# NOTE: no --force-refetch. With a complete cache this runs fully offline.

RT="$OUT/boundary_enrichment/run_table.csv"

echo "== [2/5] distance profile + random-boundary control"
python scripts/distance_profile.py \
  --run-table "$RT" --outdir "$OUT/distance_profile" \
  --thresholds 5,10,15 --bins 0,25,50,100,200,400,800,1600,inf

python scripts/random_boundary_control.py \
  --run-table "$RT" --outdir "$OUT/random_boundary_control" \
  --thresholds 5,10,15 --n-random 100 --save-null-profiles

# -----------------------------------------------------------------------------
# Stage 2. Composition + pseudoreplication controls  (Table 1, Section 3.1)
# -----------------------------------------------------------------------------
echo "== [3/5] composition-matched and per-genome controls (Table 1)"
python scripts/boundary_composition_controls.py \
  --run-table "$RT" --outdir "$OUT/boundary_composition_controls" \
  --k 5 --near 25 --windows 15,50

python scripts/per_genome_paired_test.py \
  --run-table "$RT" --outdir "$OUT/per_genome_paired" \
  --k 5 --near 25 --n-null 200 --seed 12345

# -----------------------------------------------------------------------------
# Stage 3. 3' positional / motif analysis  (Figure 1)
# -----------------------------------------------------------------------------
echo "== [4/5] 3-prime positional and motif analysis (Figure 1)"
python scripts/threeprime_motifs.py \
  --gb-cache "$GB" --outdir "$OUT/threeprime_motifs" \
  --window 50 --min-run 5

# -----------------------------------------------------------------------------
# Stage 4. Terminator benchmark  (Table 2, Figure 2)
# -----------------------------------------------------------------------------
echo "== [5/5] terminator benchmark"

# 4a. predictor inputs
python scripts/gb_to_transterm.py --gb-cache "$GB"      --outdir "$OUT/transterm_inputs"
python scripts/gb_to_transterm.py --gb-cache "$GB_MONO" --outdir "$OUT/transterm_inputs_mono"

# 4b. TransTermHP
transterm -p "$EXPTERM" "$OUT/transterm_inputs/genomes.fasta"  "$OUT/transterm_inputs/genomes.crd"  > "$OUT/transterm_all.tt"
transterm -p "$EXPTERM" "$OUT/transterm_inputs_mono/genomes.fasta" "$OUT/transterm_inputs_mono/genomes.crd" > "$OUT/transterm_mono.tt"

# 4c. RNIE (Infernal covariance-model search)
cmsearch --tblout "$OUT/rnie_all.tbl"  -T 14 "$RNIE_CM" "$OUT/transterm_inputs/genomes.fasta"
cmsearch --tblout "$OUT/rnie_mono.tbl" -T 14 "$RNIE_CM" "$OUT/transterm_inputs_mono/genomes.fasta"
python scripts/predictor_to_bed.py --input "$OUT/rnie_all.tbl"  --format cmsearch --output "$OUT/reference_rnie_all.bed"  --min-score 14
python scripts/predictor_to_bed.py --input "$OUT/rnie_mono.tbl" --format cmsearch --output "$OUT/reference_rnie_mono.bed" --min-score 14

# 4d. the four benchmark runs (Table 2)
python scripts/terminator_benchmark.py --gb-cache "$GB" \
  --transterm "$OUT/transterm_all.tt" --outdir "$OUT/bench_all_transterm" \
  --k-tract 6 --sweep 4,5,6,7,8,10 --null-per-genome 300 --gene-end-mode same

python scripts/terminator_benchmark.py --gb-cache "$GB" \
  --reference-bed "$OUT/reference_rnie_all.bed" --outdir "$OUT/bench_all_rnie" \
  --k-tract 6 --sweep 4,5,6,7,8,10 --null-per-genome 300 --gene-end-mode same

python scripts/terminator_benchmark.py --gb-cache "$GB_MONO" \
  --transterm "$OUT/transterm_mono.tt" --outdir "$OUT/bench_mono_transterm" \
  --k-tract 6 --sweep 4,5,6,7,8,10 --null-per-genome 300 --gene-end-mode any

python scripts/terminator_benchmark.py --gb-cache "$GB_MONO" \
  --reference-bed "$OUT/reference_rnie_mono.bed" --outdir "$OUT/bench_mono_rnie" \
  --k-tract 6 --sweep 4,5,6,7,8,10 --null-per-genome 300 --gene-end-mode any

# 4e. strand-model cross-check and k-stability (Section 3.4)
python scripts/terminator_benchmark.py --gb-cache "$GB_MONO" \
  --transterm "$OUT/transterm_mono.tt" --outdir "$OUT/bench_mono_polyA" \
  --k-tract 6 --sweep 4,5,6,7,8,10 --null-per-genome 300 --gene-end-mode polyA

for K in 5 6 7; do
  python scripts/terminator_benchmark.py --gb-cache "$GB_MONO" \
    --transterm "$OUT/transterm_mono.tt" --outdir "$OUT/bench_mono_k${K}" \
    --k-tract "${K}" --sweep "${K}" --null-per-genome 300 --gene-end-mode any
done

# 4f. per-clade gene-end enrichment (Figure 2a)
python scripts/clade_enrichment_table.py \
  --bench-dirs "$OUT/bench_mono_k5,$OUT/bench_mono_k6,$OUT/bench_mono_k7" \
  --gb-cache "$GB_MONO" \
  --outdir "$OUT/clade_enrichment"

echo
echo "== done. Key outputs:"
echo "   Table 1 : $OUT/boundary_composition_controls/table1_boundary_enrichment.md"
echo "             $OUT/per_genome_paired/per_genome_paired_summary.csv"
echo "   Table 2 : $OUT/bench_*/terminator_benchmark_summary.csv"
echo "   Fig 1   : $OUT/threeprime_motifs/"
echo "   Fig 2   : $OUT/clade_enrichment/"
