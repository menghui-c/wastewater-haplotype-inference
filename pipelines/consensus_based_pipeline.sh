#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# consensus_based_pipeline.sh
#
# A one‐step, parameterized MAFFT‐based pipeline:
#   0) Copy over the read‐based variant table from read_based_results_files
#   1) Align sample genomes to a reference
#   2) Summarize per‐site allele counts
#   3) Enumerate & normalize all haplotypes
#
# Usage:
#   chmod +x consensus_based_pipeline.sh
#   ./consensus_based_pipeline.sh \
#     -r reference.fasta \
#     -s samples.fasta \
#     --readdir path/to/read_based_output \
#     [-t threads] \
#     [-o outdir]
# ==============================================================================

# — parse arguments —
while [[ $# -gt 0 ]]; do
  case $1 in
    -r|--reference) REF="$2"; shift 2;;
    -s|--samples)   SAMPLES="$2"; shift 2;;
    -t|--threads)   THREADS="$2"; shift 2;;
    -o|--outdir)    OUTDIR="$2"; shift 2;;
    --readdir)      READDIR="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done

: "${REF:?Need -r reference.fasta}"
: "${SAMPLES:?Need -s samples.fasta}"
: "${READDIR:?Need --readdir to point to read_based_results directory}"
THREADS="${THREADS:-4}"
OUTDIR="${OUTDIR:-consensus_based_results}"
mkdir -p "$OUTDIR"

echo ">>> Copying read‐based variant table from $READDIR/"
cp "${READDIR}/variants_final_table.tsv" "${OUTDIR}/variants_final_table.tsv"

echo ">>> Aligning samples to reference with MAFFT (threads=$THREADS)"
mafft --thread "$THREADS" --auto \
      --keeplength --addfragments "$SAMPLES" \
      "$REF" > "${OUTDIR}/aligned.fasta"

echo ">>> Summarizing variant alleles"
./gisaid_summarize_variant_alleles.py \
  --variants-table "${OUTDIR}/variants_final_table.tsv" \
  --aligned-fasta "${OUTDIR}/aligned.fasta" \
  --reference-fasta "$REF" \
  --output "${OUTDIR}/consensus_variant_allele_summary.tsv"

echo ">>> Enumerating & normalizing haplotypes"
./gisaid_haplotype_freq.py \
  --summary "${OUTDIR}/consensus_variant_allele_summary.tsv" \
  --output "${OUTDIR}/consensus_all_haplotypes.tsv"

echo ">>> Done! Results in $OUTDIR/"
echo "    • ${OUTDIR}/aligned.fasta"
echo "    • ${OUTDIR}/consensus_variant_allele_summary.tsv"
echo "    • ${OUTDIR}/consensus_all_haplotypes.tsv"


