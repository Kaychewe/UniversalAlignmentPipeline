#!/usr/bin/env bash

# ------------------------------------------------------------
# 06_markdup.sh
# Mark duplicate reads using Picard MarkDuplicates
#
# Usage:
#   bash 06_markdup.sh <OUTDIR> <LOGDIR> <SAMPLE> <THREADS> <INPUT_BAM>
#
# Output:
#   <OUTDIR>/<SAMPLE>/dedup/<SAMPLE>.dedup.bam
# Logs:
#   <LOGDIR>/<SAMPLE>/markdup.log
#   <LOGDIR>/<SAMPLE>/markdup.metrics.txt
# ------------------------------------------------------------

set -euo pipefail

# ---------------------- USAGE CHECK -------------------------
if [[ $# -lt 5 ]]; then
    echo "Usage: bash 06_markdup.sh <OUTDIR> <LOGDIR> <SAMPLE> <THREADS> <INPUT_BAM>"
    exit 1
fi

OUTDIR=$1
LOGDIR=$2
SAMPLE=$3
THREADS=$4
IN_BAM=$5

DEDUP_DIR="$OUTDIR/$SAMPLE/dedup"
SAMPLE_LOGDIR="$LOGDIR/$SAMPLE"

LOGFILE="$SAMPLE_LOGDIR/markdup.log"
METRICS="$SAMPLE_LOGDIR/markdup.metrics.txt"

OUT_BAM="$DEDUP_DIR/${SAMPLE}.dedup.bam"

mkdir -p "$DEDUP_DIR"
mkdir -p "$SAMPLE_LOGDIR"

# ---------------------- Conda Init ---------------------------
source /home/kchewe/anaconda3/etc/profile.d/conda.sh
conda activate bissnp_benchmark

# ---------------------- Tool Check ---------------------------
if ! command -v picard &> /dev/null; then
    echo "[ERROR] picard command not found. Ensure it's in the environment." | tee -a "$LOGFILE"
    exit 1
fi

# ---------------------- Log Header ---------------------------
{
echo "------------------------------------------------------------"
echo " [MARKDUP] Marking duplicates for sample: $SAMPLE"
echo " Timestamp: $(date)"
echo " Input BAM:  $IN_BAM"
echo " Output BAM: $OUT_BAM"
echo " Metrics:    $METRICS"
echo " Threads:    $THREADS (Picard is single-threaded for MarkDuplicates)"
echo " Log:        $LOGFILE"
echo "------------------------------------------------------------"
} | tee "$LOGFILE"

# ---------------------- Run Picard ---------------------------
picard MarkDuplicates \
    I="$IN_BAM" \
    O="$OUT_BAM" \
    M="$METRICS" \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=false \
    VALIDATION_STRINGENCY=SILENT \
    2>&1 | tee -a "$LOGFILE"

# ---------------------- Finish -------------------------------
{
echo "[MARKDUP] Completed duplicate marking for $SAMPLE"
echo " Output BAM: $OUT_BAM"
echo " Metrics:    $METRICS"
echo "------------------------------------------------------------"
} | tee -a "$LOGFILE"
