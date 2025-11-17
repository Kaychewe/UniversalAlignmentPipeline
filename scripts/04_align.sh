#!/usr/bin/env bash

# ------------------------------------------------------------
# 03_align.sh
# Perform alignment using bwa-mem2 or bismark (SE or PE)
#
# Usage:
#   bash 03_align.sh <OUTDIR> <LOGDIR> <SAMPLE> <ALIGNER> <THREADS> <REF> <REF_DIR> <R1> [R2]
#
# Output:
#   <OUTDIR>/<SAMPLE>/aligned/<SAMPLE>.unsorted.bam
# Log:
#   <LOGDIR>/<SAMPLE>/align.log
# ------------------------------------------------------------

# TODO: write function to looks through current Snakefile dir and tmp dir 

set -euo pipefail

# ---------------------- USAGE CHECK -------------------------
if [[ $# -lt 9 ]]; then
    echo "Usage: bash 03_align.sh <OUTDIR> <LOGDIR> <SAMPLE> <ALIGNER> <THREADS> <REF> <REF_DIR> <R1> [R2]"
    exit 1
fi

OUTDIR=$1
LOGDIR=$2
SAMPLE=$3
ALIGNER=$4
THREADS=$5
REF=$6
REF_DIR=$7
R1=$8
R2=${9:-}

ALIGN_DIR="$OUTDIR/$SAMPLE/aligned"
SAMPLE_LOGDIR="$LOGDIR/$SAMPLE"
LOGFILE="$SAMPLE_LOGDIR/align.log"

mkdir -p "$ALIGN_DIR"
mkdir -p "$SAMPLE_LOGDIR"

# ---------------------- Conda Init ---------------------------
source /home/kchewe/anaconda3/etc/profile.d/conda.sh
conda activate bissnp_benchmark

# ---------------------- Check Tools --------------------------
if ! command -v samtools &> /dev/null; then
    echo "[ERROR] samtools not found." | tee -a "$LOGFILE"
    exit 1
fi

if [[ "$ALIGNER" == "bwa-mem2" ]] && ! command -v bwa-mem2 &> /dev/null; then
    echo "[ERROR] bwa-mem2 not found." | tee -a "$LOGFILE"
    exit 1
fi

if [[ "$ALIGNER" == "bismark" ]] && ! command -v bismark &> /dev/null; then
    echo "[ERROR] bismark not found." | tee -a "$LOGFILE"
    exit 1
fi

# ---------------------- Output BAM ---------------------------
OUT_BAM="$ALIGN_DIR/${SAMPLE}.unsorted.bam"

# ---------------------- Log Header ---------------------------
{
echo "------------------------------------------------------------"
echo " [ALIGN] Starting alignment for sample: $SAMPLE"
echo " Timestamp: $(date)"
echo " Aligner: $ALIGNER"
echo " Threads: $THREADS"
echo " Reference FASTA: $REF"
echo " Reference DIR:   $REF_DIR"
echo " R1: $R1"
[[ -n "$R2" ]] && echo " R2: $R2"
echo " Output BAM: $OUT_BAM"
echo " Log: $LOGFILE"
echo "------------------------------------------------------------"
} | tee "$LOGFILE"

# ---------------------- Alignment Execution ------------------

if [[ "$ALIGNER" == "bwa-mem2" ]]; then
    # ---------------------------------------------------------
    # BWA-MEM2 alignment
    # ---------------------------------------------------------
    if [[ -n "$R2" ]]; then
        echo "[ALIGN] Running bwa-mem2 (paired-end)" | tee -a "$LOGFILE"

        bwa-mem2 mem -t "$THREADS" "$REF" "$R1" "$R2" \
            2> >(tee -a "$LOGFILE" >&2) \
            | samtools view -bS - > "$OUT_BAM"

    else
        echo "[ALIGN] Running bwa-mem2 (single-end)" | tee -a "$LOGFILE"

        bwa-mem2 mem -t "$THREADS" "$REF" "$R1" \
            2> >(tee -a "$LOGFILE" >&2) \
            | samtools view -bS - > "$OUT_BAM"
    fi

elif [[ "$ALIGNER" == "bismark" ]]; then
    # ---------------------------------------------------------
    # BISMARK WGBS alignment (with proper multithreading)
    # ---------------------------------------------------------
    echo "[ALIGN] Running Bismark with --parallel $THREADS" | tee -a "$LOGFILE"

    if [[ -n "$R2" ]]; then
        # -------- Paired-end --------
        bismark \
            --genome "$REF_DIR" \
            -1 "$R1" -2 "$R2" \
            -o "$ALIGN_DIR" \
            --bam \
            --parallel "$THREADS" \
            2>&1 | tee -a "$LOGFILE"

        mv "$ALIGN_DIR"/*_pe.bam "$OUT_BAM"

    else
        # -------- Single-end --------
        bismark \
            --genome "$REF_DIR" \
            "$R1" \
            -o "$ALIGN_DIR" \
            --bam \
            --parallel "$THREADS" \
            2>&1 | tee -a "$LOGFILE"

        mv "$ALIGN_DIR"/*_se.bam "$OUT_BAM"
    fi


else
    echo "[ERROR] Unknown aligner: $ALIGNER" | tee -a "$LOGFILE"
    exit 1
fi

# ---------------------- Finish -------------------------------
{
echo "[ALIGN] Completed alignment for $SAMPLE"
echo " Output BAM: $OUT_BAM"
echo "------------------------------------------------------------"
} | tee -a "$LOGFILE"
