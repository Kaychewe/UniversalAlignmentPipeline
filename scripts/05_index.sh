#!/usr/bin/env bash

# ------------------------------------------------------------
# 05_index.sh
# Index a sorted BAM using samtools
#
# Usage:
#   bash 05_index.sh <OUTDIR> <LOGDIR> <SAMPLE> <THREADS> <INPUT_BAM>
#
# Output:
#   <OUTDIR>/<SAMPLE>/sorted/<SAMPLE>.sorted.bam.bai
# Logs:
#   <LOGDIR>/<SAMPLE>/index.log
# ------------------------------------------------------------

set -euo pipefail

# ---------------------- USAGE CHECK -------------------------
if [[ $# -lt 5 ]]; then
    echo "Usage: bash 05_index.sh <OUTDIR> <LOGDIR> <SAMPLE> <THREADS> <INPUT_BAM>"
    exit 1
fi

OUTDIR=$1
LOGDIR=$2
SAMPLE=$3
THREADS=$4
IN_BAM=$5

SORT_DIR="$OUTDIR/$SAMPLE/sorted"
SAMPLE_LOGDIR="$LOGDIR/$SAMPLE"
LOGFILE="$SAMPLE_LOGDIR/index.log"

OUT_INDEX="${IN_BAM}.bai"

mkdir -p "$SORT_DIR"
mkdir -p "$SAMPLE_LOGDIR"

# ---------------------- Conda Init ---------------------------
source /home/kchewe/anaconda3/etc/profile.d/conda.sh
conda activate bissnp_benchmark

# ---------------------- Tool Check ---------------------------
if ! command -v samtools &> /dev/null; then
    echo "[ERROR] samtools not found." | tee -a "$LOGFILE"
    exit 1
fi

# ---------------------- Log Header ---------------------------
{
echo "------------------------------------------------------------"
echo " [INDEX] Indexing BAM for sample: $SAMPLE"
echo " Timestamp: $(date)"
echo " Input BAM:  $IN_BAM"
echo " Output BAI: $OUT_INDEX"
echo " Threads:    $THREADS"
echo " Log:        $LOGFILE"
echo "------------------------------------------------------------"
} | tee "$LOGFILE"

# ---------------------- Indexing ------------------------------
samtools index -@ "$THREADS" "$IN_BAM" 2>&1 | tee -a "$LOGFILE"

# ---------------------- Finish -------------------------------
{
echo "[INDEX] Completed BAM indexing for $SAMPLE"
echo " Output index: $OUT_INDEX"
echo "------------------------------------------------------------"
} | tee -a "$LOGFILE"
