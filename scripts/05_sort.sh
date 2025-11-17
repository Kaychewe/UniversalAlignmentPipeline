#!/usr/bin/env bash

# ------------------------------------------------------------
# 04_sort.sh
# Sort a BAM file with samtools
#
# Usage:
#   bash 04_sort.sh <OUTDIR> <LOGDIR> <SAMPLE> <THREADS> <INPUT_BAM>
#
# Output:
#   <OUTDIR>/<SAMPLE>/sorted/<SAMPLE>.sorted.bam
# Logs:
#   <LOGDIR>/<SAMPLE>/sort.log
# ------------------------------------------------------------

set -euo pipefail

# ---------------------- USAGE CHECK -------------------------
if [[ $# -lt 5 ]]; then
    echo "Usage: bash 04_sort.sh <OUTDIR> <LOGDIR> <SAMPLE> <THREADS> <INPUT_BAM>"
    exit 1
fi

OUTDIR=$1
LOGDIR=$2
SAMPLE=$3
THREADS=$4
IN_BAM=$5

SORT_DIR="$OUTDIR/$SAMPLE/sorted"
SAMPLE_LOGDIR="$LOGDIR/$SAMPLE"
LOGFILE="$SAMPLE_LOGDIR/sort.log"

OUT_SORTED="$SORT_DIR/${SAMPLE}.sorted.bam"

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
echo " [SORT] Sorting BAM for sample: $SAMPLE"
echo " Timestamp: $(date)"
echo " Input BAM:  $IN_BAM"
echo " Output BAM: $OUT_SORTED"
echo " Threads:    $THREADS"
echo " Log:        $LOGFILE"
echo "------------------------------------------------------------"
} | tee "$LOGFILE"

# ---------------------- Sorting ------------------------------
samtools sort \
    -@ "$THREADS" \
    -o "$OUT_SORTED" \
    "$IN_BAM" \
    2>&1 | tee -a "$LOGFILE"

# ---------------------- Finish -------------------------------
{
echo "[SORT] Completed BAM sorting for $SAMPLE"
echo " Output: $OUT_SORTED"
echo "------------------------------------------------------------"
} | tee -a "$LOGFILE"
