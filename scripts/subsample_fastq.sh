#!/usr/bin/env bash

# -----------------------------------------------------
# Subsample FASTQ using BBTools reformat.sh
# Author: Kasonde Chewe
# -----------------------------------------------------

set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: bash subsample_reformat.sh <fraction> <outdir> <R1> [R2]"
    exit 1
fi

FRACTION=$1
OUTDIR=$2
R1=$3
R2=${4:-""}
SAMPLE=$5

mkdir -p "$OUTDIR"

# ---- Create log file with timestamp ----------------------------------
CURRENT_DATETIME=$(date +"%Y-%m-%d_%H-%M-%S")

R1_BASE=$(basename "$R1")
R1_BASE="${R1_BASE%.fastq.gz}"
R1_BASE="${R1_BASE%.fq.gz}"
R1_BASE="${R1_BASE%.fastq}"
R1_BASE="${R1_BASE%.fq}"

LOG="$OUTDIR/subsample_${R1_BASE}_${FRACTION}_${CURRENT_DATETIME}.log"

# ---- Initialize Conda ------------------------------------
source /home/kchewe/anaconda3/etc/profile.d/conda.sh
conda activate bissnp_benchmark

# # ---- Clean output names ----------------------------------------------
# clean_name() {
#     local f=$(basename "$1")
#     f="${f%.fastq.gz}"
#     f="${f%.fq.gz}"
#     f="${f%.fastq}"
#     f="${f%.fq}"
#     echo "${f}_sub.fastq.gz"
# }

# R1_OUT="$OUTDIR/$(clean_name "$R1")"
# [[ -n "$R2" ]] && R2_OUT="$OUTDIR/$(clean_name "$R2")"

# Override output names using SAMPLE
R1_OUT="$OUTDIR/${SAMPLE}_R1_sub.fastq.gz"
R2_OUT="$OUTDIR/${SAMPLE}_R2_sub.fastq.gz"

# ---- Header ----------------------------------------------------------
{
echo "-----------------------------------------------------"
echo " Subsampling FASTQ using reformat.sh"
echo " Timestamp:   $CURRENT_DATETIME"
echo " Fraction:    $FRACTION"
echo " R1:          $R1"
[[ -n "$R2" ]] && echo " R2:          $R2"
echo " Output R1:   $R1_OUT"
[[ -n "$R2" ]] && echo " Output R2:   $R2_OUT"
echo " Log file:    $LOG"
echo "-----------------------------------------------------"
} | tee "$LOG"

# ---- Run Single-end or Paired-end ----------------------------------
if [[ -z "$R2" ]]; then
    echo "[Mode] Single-end" | tee -a "$LOG"

    reformat.sh \
        in="$R1" \
        out="$R1_OUT" \
        samplerate="$FRACTION" \
        overwrite=t 2>&1 | tee -a "$LOG"

else
    echo "[Mode] Paired-end" | tee -a "$LOG"

    reformat.sh \
        in="$R1" in2="$R2" \
        out="$R1_OUT" out2="$R2_OUT" \
        samplerate="$FRACTION" \
        overwrite=t 2>&1 | tee -a "$LOG"
fi

echo "Done!" | tee -a "$LOG"
echo "-----------------------------------------------------" | tee -a "$LOG"
