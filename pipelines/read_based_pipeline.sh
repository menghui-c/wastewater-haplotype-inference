#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# read_based_pipeline.sh  (paired-end + optional iVar primer trimming)
#
# Steps:
#   1) Index reference (bwa, samtools faidx)
#   2) Map reads (bwa mem) → sort & index BAM
#   3) [optional] iVar primer trimming (requires --primer-bed)
#   4) Indel-quality recalibration (lofreq indelqual --dindel)
#   5) SNV calling (lofreq call)
#   6) Secondary VCF filtering (your python script)
#   7) Build variant table (your python script)
#   8) Extract SNV list (your python script)
#   9) Encode reads → pattern counts (your python script)
#  10) Enumerate & normalize haplotypes (your python script)
#
# Usage examples:
#   # pair-ended reads（recommended）：
#   ./read_based_pipeline.sh \
#     -r NC_045512.2.fasta \
#     -1 xxxx_R1.fastq \
#     -2 xxxx_R2.fastq \
#     --primer-bed primers.bed \
#     -t 8 --min-cov 30 --min-alt-freq 0.10 --min-ref-freq 0.10 \
#     --outdir results/xxxx
#
#   # Single-file interleaving（interleaved）：
#   ./read_based_pipeline.sh \
#     -r NC_045512.2.fasta \
#     -q SRR30606112.fastq \
#     --interleaved \
#     --primer-bed primers.bed \
#     -t 8 --outdir results/xxxx
#
#   # single-ended：
#   ./read_based_pipeline.sh -r ref.fa -q reads.fq -t 8 --outdir results/run1
# -----------------------------------------------------------------------------

# ---------- Parse args ----------
REF="" ; READS="" ; R1="" ; R2=""
THREADS=8
MIN_COV=30
MIN_ALT=0.10
MIN_REF=0.10
OUTDIR="read_based_results_files"
PRIMER_BED=""
INTERLEAVED="false"
IVAR_Q=20          
IVAR_MINLEN=30    
IVAR_EXCLUDE="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--reference)     REF="$2"; shift 2;;
    -q|--reads)         READS="$2"; shift 2;;
    -1)                 R1="$2"; shift 2;;
    -2)                 R2="$2"; shift 2;;
    --interleaved)      INTERLEAVED="true"; shift 1;;
    --primer-bed)       PRIMER_BED="$2"; shift 2;;
    --ivar-q)           IVAR_Q="$2"; shift 2;;
    --ivar-minlen)      IVAR_MINLEN="$2"; shift 2;;
    --ivar-exclude)     IVAR_EXCLUDE="true"; shift 1;;
    -t|--threads)       THREADS="$2"; shift 2;;
    --min-cov)          MIN_COV="$2"; shift 2;;
    --min-alt-freq)     MIN_ALT="$2"; shift 2;;
    --min-ref-freq)     MIN_REF="$2"; shift 2;;
    --outdir)           OUTDIR="$2"; shift 2;;
    -h|--help)
      sed -n '1,120p' "$0"; exit 0;;
    *) echo "[ERROR] Unknown arg: $1"; exit 1;;
  esac
done

# ---------- Checks ----------
[[ -z "$REF" ]] && { echo "[ERROR] Need -r/--reference"; exit 1; }

if [[ -n "$R1" || -n "$R2" ]]; then
  [[ -z "$R1" || -z "$R2" ]] && { echo "[ERROR] Need both -1 and -2 for paired-end."; exit 1; }
elif [[ -n "$READS" ]]; then
  : 
else
  echo "[ERROR] Provide either (-1 R1 -2 R2) or (-q READS)."; exit 1
fi

# Tools check
for x in bwa samtools lofreq; do
  command -v "$x" >/dev/null 2>&1 || { echo "[ERROR] $x not found in PATH"; exit 1; }
done
if [[ -n "$PRIMER_BED" ]]; then
  command -v ivar >/dev/null 2>&1 || { echo "[ERROR] iVar not found but --primer-bed specified"; exit 1; }
fi

mkdir -p "$OUTDIR"

if [[ -n "$R1" && -n "$R2" ]]; then
  bn="$(basename "$R1")"
  bn="${bn%%.fastq.gz}"; bn="${bn%%.fq.gz}"; bn="${bn%%.fastq}"; bn="${bn%%.fq}"
  PREFIX="${bn%_R1*}"
elif [[ -n "$READS" ]]; then
  bn="$(basename "$READS")"
  bn="${bn%%.fastq.gz}"; bn="${bn%%.fq.gz}"; bn="${bn%%.fastq}"; bn="${bn%%.fq}"
  PREFIX="$bn"
fi

LOG="$OUTDIR/${PREFIX}.pipeline.log"
exec > >(tee -i "$LOG") 2>&1

echo "========== PIPELINE START =========="
date
echo "[INFO] Reference         : $REF"
echo "[INFO] Threads           : $THREADS"
if [[ -n "$R1" && -n "$R2" ]]; then
  echo "[INFO] Mode              : Paired-end (R1/R2)"
  echo "[INFO] R1                : $R1"
  echo "[INFO] R2                : $R2"
elif [[ "$INTERLEAVED" == "true" ]]; then
  echo "[INFO] Mode              : Paired-end (single interleaved file)"
  echo "[INFO] Reads             : $READS"
else
  echo "[INFO] Mode              : Single-end"
  echo "[INFO] Reads             : $READS"
fi
if [[ -n "$PRIMER_BED" ]]; then
  echo "[INFO] Primer BED        : $PRIMER_BED"
  echo "[INFO] iVar trim params  : -q $IVAR_Q -m $IVAR_MINLEN $([[ $IVAR_EXCLUDE == "true" ]] && echo "-e" || echo "")"
else
  echo "[INFO] Primer trimming   : SKIP (no --primer-bed)"
fi
echo "[INFO] LoFreq thresholds : --min-cov $MIN_COV --min-alt-freq $MIN_ALT --min-ref-freq $MIN_REF"
echo "[INFO] Output dir        : $OUTDIR"
echo "===================================="

# ---------- 1) Index reference ----------
echo ">>> [1/9] Indexing reference"
[[ -f "${REF}.bwt" ]] || bwa index "$REF"
[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

# ---------- 2) Mapping ----------
echo ">>> [2/9] Mapping with bwa mem (threads=$THREADS)"
SAM="$OUTDIR/${PREFIX}.sam"
BAM_SORT="$OUTDIR/${PREFIX}.sorted.bam"

if [[ -n "$R1" && -n "$R2" ]]; then
  bwa mem -t "$THREADS" "$REF" "$R1" "$R2" > "$SAM"
elif [[ "$INTERLEAVED" == "true" ]]; then
  bwa mem -t "$THREADS" -p "$REF" "$READS" > "$SAM"
else
  bwa mem -t "$THREADS" "$REF" "$READS" > "$SAM"
fi

samtools view -bS "$SAM" | samtools sort -@ "$THREADS" -o "$BAM_SORT"
samtools index "$BAM_SORT"

echo ">>> Mapping summary:"
samtools flagstat "$BAM_SORT" | tee "$OUTDIR/${PREFIX}.mapping_stats.txt"

TOTAL=$(samtools flagstat "$BAM_SORT" | head -n 1 | awk '{print $1}')
MAPPED=$(samtools flagstat "$BAM_SORT" | awk '/ mapped \(/ {print $1; exit}')
PCT=$(awk "BEGIN{ if($TOTAL>0) printf \"%.2f\", 100*$MAPPED/$TOTAL; else print \"0.00\" }")
echo ">>> Total reads         : $TOTAL"
echo ">>> Successfully mapped : $MAPPED ($PCT%)"

# ---------- 3) [optional] iVar primer trimming ----------
WORK_BAM="$BAM_SORT"
if [[ -n "$PRIMER_BED" ]]; then
  echo ">>> [3/9] iVar primer trimming"
  BAM_NSORT="$OUTDIR/${PREFIX}.namesort.bam"
  samtools sort -n -@ "$THREADS" -o "$BAM_NSORT" "$BAM_SORT"

  IVAR_OUT_PREFIX="$OUTDIR/${PREFIX}.trimmed"
  IVAR_CMD=(ivar trim -i "$BAM_NSORT" -b "$PRIMER_BED" -p "$IVAR_OUT_PREFIX" -q "$IVAR_Q" -m "$IVAR_MINLEN")
  [[ "$IVAR_EXCLUDE" == "true" ]] && IVAR_CMD+=(-e)

  echo ">>> Running: ${IVAR_CMD[*]}"
  "${IVAR_CMD[@]}"

  BAM_TRIM_RAW="${IVAR_OUT_PREFIX}.bam"
  BAM_TRIM_SORT="$OUTDIR/${PREFIX}.trimmed.csorted.bam"
  samtools sort -@ "$THREADS" -o "$BAM_TRIM_SORT" "$BAM_TRIM_RAW"
  samtools index "$BAM_TRIM_SORT"

  WORK_BAM="$BAM_TRIM_SORT"
  echo ">>> Primer-trimmed BAM  : $WORK_BAM"
else
  echo ">>> [3/9] iVar primer trimming : SKIP"
fi

# ---------- 4) Indel-quality recalibration ----------
echo ">>> [4/9] LoFreq indelqual --dindel"
BAM_INDEL="$OUTDIR/${PREFIX}.indelqual.bam"
rm -f "$BAM_INDEL" "${BAM_INDEL}.bai"
lofreq indelqual --dindel -f "$REF" -o "$BAM_INDEL" "$WORK_BAM"
samtools index "$BAM_INDEL"

# ---------- 5) SNV calling ----------
echo ">>> [5/9] LoFreq SNV calling"
VCF_RAW="$OUTDIR/${PREFIX}.lofreq.raw.vcf"
rm -f "$VCF_RAW"
lofreq call --min-cov "$MIN_COV" --sig 1e-6 -f "$REF" -o "$VCF_RAW" "$BAM_INDEL"

# ---------- 6) VCF secondary filtering ----------
echo ">>> [6/9] Secondary VCF filtering (Python)"
VCF_FILT="$OUTDIR/${PREFIX}.lofreq.filtered.vcf"
./filter_lofreq_vcf.py \
  -i "$VCF_RAW" \
  -o "$VCF_FILT" \
  --min-cov "$MIN_COV" \
  --min-alt-freq "$MIN_ALT" \
  --min-ref-freq "$MIN_REF"

# ---------- 7) Variant table ----------
echo ">>> [7/9] Build variant table (Python)"
VAR_TSV="$OUTDIR/variants_final_table.tsv"
./extract_snv_table.py -i "$VCF_FILT" -o "$VAR_TSV"

# ---------- 8) SNV list ----------
echo ">>> [8/9] Extract SNV list (Python)"
SNV_JSON="$OUTDIR/snv_list.json"
./get_snvs.py -i "$VCF_FILT" -o "$SNV_JSON"

# ---------- 9) Encode reads → pattern counts ----------
echo ">>> [9/9] Encode reads → pattern counts (Python)"
ENC_COUNTS="$OUTDIR/encoding_counts.tsv"
./read_encoding.py -s "$SNV_JSON" -b "$BAM_INDEL" -o "$ENC_COUNTS"

# ---------- 10) Enumerate & normalize haplotypes ----------
echo ">>> [+] Enumerate & normalize haplotypes (Python)"
HAP_FREQS="$OUTDIR/all_hap_allele_freqs.tsv"
./compute_haplotype_freqs_allele.py --counts "$ENC_COUNTS" -s "$SNV_JSON" -o "$HAP_FREQS"

# ---------- Wrap up ----------
echo ">>> Outputs:"
ls -lh "$OUTDIR" | sed 's/^/    /'

echo "========== PIPELINE DONE =========="
date