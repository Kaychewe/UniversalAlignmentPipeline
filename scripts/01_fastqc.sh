#!/usr/bin/env bash

# -------------------------------------------------------------
# 01_fastqc.sh
# Run FastQC for a sample (SE or PE)
#
# Usage:
#   bash 01_fastqc.sh <OUTDIR> <LOGDIR> <SAMPLE> <THREADS> <FASTQ_R1> [FASTQ_R2]
#
# Output:
#   <OUTDIR>/<SAMPLE>/fastqc/
# Logs:
#   <LOGDIR>/<SAMPLE>/fastqc.log
# -------------------------------------------------------------

set -euo pipefail

# ---------------------- USAGE CHECK --------------------------
if [[ $# -lt 5 ]]; then
    echo "Usage: bash 01_fastqc.sh <OUTDIR> <LOGDIR> <SAMPLE> <THREADS> <FASTQ_R1> [FASTQ_R2]"
    exit 1
fi

OUTDIR=$1
LOGDIR=$2
SAMPLE=$3
THREADS=$4
R1=$5
R2=${6:-}

FASTQC_DIR="$OUTDIR/$SAMPLE/fastqc"
SAMPLE_LOGDIR="$LOGDIR/$SAMPLE"
LOGFILE="$SAMPLE_LOGDIR/fastqc.log"

mkdir -p "$FASTQC_DIR"
mkdir -p "$SAMPLE_LOGDIR"

# ---------------------- Conda Init ---------------------------
source /home/kchewe/anaconda3/etc/profile.d/conda.sh
conda activate bissnp_benchmark

if ! command -v fastqc &> /dev/null; then
    echo "[ERROR] fastqc not found." | tee -a "$LOGFILE"
    exit 1
fi

# ---------------------- LOG HEADER ---------------------------
{
echo "------------------------------------------------------------"
echo " [FastQC] Sample:  $SAMPLE"
echo " Timestamp:         $(date)"
echo " R1:                $R1"
[[ -n "$R2" ]] && echo " R2:                $R2"
echo " Threads:           $THREADS"
echo " Output dir:        $FASTQC_DIR"
echo " Log file:          $LOGFILE"
echo "------------------------------------------------------------"
} | tee "$LOGFILE"

# ---------------------- Run FastQC ---------------------------

fastqc --threads "$THREADS" -o "$FASTQC_DIR" "$R1" 2>&1 | tee -a "$LOGFILE"

if [[ -n "$R2" ]]; then
    fastqc --threads "$THREADS" -o "$FASTQC_DIR" "$R2" 2>&1 | tee -a "$LOGFILE"
fi

{
echo "[FastQC] Completed for sample $SAMPLE"
echo " Outputs stored in: $FASTQC_DIR"
echo "------------------------------------------------------------"
} | tee -a "$LOGFILE"
