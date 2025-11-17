#!/usr/bin/env bash

# ------------------------------------------------------------
# 02_fastp.sh
# Run fastp for filtering/trimming (SE or PE)
#
# Usage:
#   bash 02_fastp.sh <OUTDIR> <LOGDIR> <SAMPLE> <THREADS> <FASTQ_R1> [FASTQ_R2]
#
# Output:
#   <OUTDIR>/<SAMPLE>/trimmed/
# Logs:
#   <LOGDIR>/<SAMPLE>/fastp.log
# ------------------------------------------------------------

set -euo pipefail

# ---------------------- USAGE CHECK -------------------------
if [[ $# -lt 5 ]]; then
    echo "Usage: bash 02_fastp.sh <OUTDIR> <LOGDIR> <SAMPLE> <THREADS> <FASTQ_R1> [FASTQ_R2]"
    exit 1
fi

OUTDIR=$1       # root output directory
LOGDIR=$2       # root log directory
SAMPLE=$3       # sample name
THREADS=$4      # threads from config
R1=$5           # FASTQ R1
R2=${6:-}       # FASTQ R2 (optional)

TRIM_DIR="$OUTDIR/$SAMPLE/trimmed"
SAMPLE_LOGDIR="$LOGDIR/$SAMPLE"
LOGFILE="$SAMPLE_LOGDIR/fastp.log"

mkdir -p "$TRIM_DIR"
mkdir -p "$SAMPLE_LOGDIR"

# ---------------------- Conda Init ---------------------------
source /home/kchewe/anaconda3/etc/profile.d/conda.sh
conda activate bissnp_benchmark

if ! command -v fastp &> /dev/null; then
    echo "[ERROR] fastp not found." | tee -a "$LOGFILE"
    exit 1
fi

# ---------------------- Output Naming ------------------------
clean_name() {
    local fname
    fname=$(basename "$1")

    # Remove known FASTQ extensions
    fname="${fname%.fastq.gz}"
    fname="${fname%.fq.gz}"
    fname="${fname%.fastq}"
    fname="${fname%.fq}"

    # Remove optional subsampling tag
    fname="${fname%_sub}"

    # Add unified clean suffix
    echo "${fname}_clean.fastq.gz"
}


R1_OUT="$TRIM_DIR/$(clean_name "$R1")"
JSON="$TRIM_DIR/fastp.json"
HTML="$TRIM_DIR/fastp.html"

# ---------------------- LOG HEADER ---------------------------
{
echo "------------------------------------------------------------"
echo " [fastp] Start trimming for sample: $SAMPLE"
echo " Timestamp: $(date)"
echo " R1: $R1"
[[ -n "$R2" ]] && echo " R2: $R2"
echo " Threads: $THREADS"
echo " Output: $TRIM_DIR"
echo " Logs:   $LOGFILE"
echo "------------------------------------------------------------"
} | tee "$LOGFILE"

# ---------------------- Paired-end Mode -----------------------
if [[ -n "$R2" ]]; then
    R2_OUT="$TRIM_DIR/$(clean_name "$R2")"

    fastp \
        -i "$R1" -I "$R2" \
        -o "$R1_OUT" -O "$R2_OUT" \
        --thread "$THREADS" \
        --detect_adapter_for_pe \
        --html "$HTML" \
        --json "$JSON" 2>&1 | tee -a "$LOGFILE"

    {
    echo "[fastp] Completed paired-end trimming for $SAMPLE"
    echo "  Cleaned R1: $R1_OUT"
    echo "  Cleaned R2: $R2_OUT"
    echo "  Reports: $HTML, $JSON"
    } | tee -a "$LOGFILE"
    exit 0
fi

# ---------------------- Single-end Mode -----------------------
fastp \
    -i "$R1" \
    -o "$R1_OUT" \
    --thread "$THREADS" \
    --html "$HTML" \
    --json "$JSON" 2>&1 | tee -a "$LOGFILE"

{
echo "[fastp] Completed single-end trimming for $SAMPLE"
echo "  Cleaned R1: $R1_OUT"
echo "  Reports: $HTML, $JSON"
} | tee -a "$LOGFILE"
