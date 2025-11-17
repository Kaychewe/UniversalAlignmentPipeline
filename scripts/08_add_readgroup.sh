#!/usr/bin/env bash

# ------------------------------------------------------------
# 07_add_readgroup.sh
# Add or replace read groups in a BAM file using Picard
#
# Usage:
#   bash 06_add_readgroup.sh <OUTDIR> <LOGDIR> <SAMPLE> <THREADS> <INPUT_BAM>
#
# Output:
#   <OUTDIR>/<SAMPLE>/rgfixed/<SAMPLE>.rg.bam
# Log:
#   <LOGDIR>/<SAMPLE>/readgroup.log
# ------------------------------------------------------------

set -euo pipefail

# ---------------------- USAGE CHECK -------------------------
if [[ $# -lt 5 ]]; then
    echo "Usage: bash 06_add_readgroup.sh <OUTDIR> <LOGDIR> <SAMPLE> <THREADS> <INPUT_BAM>"
    exit 1
fi

OUTDIR=$1
LOGDIR=$2
SAMPLE=$3
THREADS=$4
IN_BAM=$5

RG_DIR="$OUTDIR/$SAMPLE/rgfixed"
SAMPLE_LOGDIR="$LOGDIR/$SAMPLE"
LOGFILE="$SAMPLE_LOGDIR/readgroup.log"

OUT_BAM="$RG_DIR/${SAMPLE}.rg.bam"

mkdir -p "$RG_DIR"
mkdir -p "$SAMPLE_LOGDIR"

# ---------------------- Conda Init ---------------------------
source /home/kchewe/anaconda3/etc/profile.d/conda.sh
conda activate bissnp_benchmark

# ---------------------- Validate Input -----------------------
if [[ ! -f "$IN_BAM" ]]; then
    echo "[ERROR] Input BAM does not exist: $IN_BAM" | tee "$LOGFILE"
    exit 1
fi

if ! command -v picard &> /dev/null; then
    echo "[ERROR] picard not found." | tee -a "$LOGFILE"
    exit 1
fi

# ---------------------- Check RG Presence ---------------------
HAS_RG=$(samtools view -H "$IN_BAM" | grep -c "^@RG" || true)

{
echo "------------------------------------------------------------"
echo " [READ GROUP] Checking read groups for sample: $SAMPLE"
echo " Timestamp: $(date)"
echo " Input BAM:  $IN_BAM"
echo " Output BAM: $OUT_BAM"
echo " Log:        $LOGFILE"
echo "------------------------------------------------------------"
} | tee "$LOGFILE"

if [[ "$HAS_RG" -gt 0 ]]; then
    echo "[READ GROUP] BAM already has read groups — copying without changes." | tee -a "$LOGFILE"
    cp "$IN_BAM" "$OUT_BAM"
    exit 0
fi

echo "[READ GROUP] No @RG found — adding read group header…" | tee -a "$LOGFILE"

# ---------------------- Construct Read Group Fields ----------
RGID="$SAMPLE"
RGSM="$SAMPLE"
RGLB="lib1"
RGPL="ILLUMINA"
RGPU="${SAMPLE}.unit1"

PICARD_CMD=(
    picard AddOrReplaceReadGroups
    I="$IN_BAM"
    O="$OUT_BAM"
    RGID="$RGID"
    RGLB="$RGLB"
    RGPL="$RGPL"
    RGSM="$RGSM"
    RGPU="$RGPU"
    VALIDATION_STRINGENCY=SILENT
)

# ---------------------- Show Full Command ---------------------
echo "[READ GROUP] Running Picard with command:" | tee -a "$LOGFILE"
echo "  ${PICARD_CMD[@]}" | tee -a "$LOGFILE"

# ---------------------- Execute Picard ------------------------
"${PICARD_CMD[@]}" \
    2>&1 | tee -a "$LOGFILE"

# ---------------------- Index BAM -----------------------------
echo "[LOG] Indexing BAM file: " | tee -a "$LOGFILE"
samtools index "$OUT_BAM" -@ "$THREADS"  | tee -a "$LOGFILE"

EXIT_CODE=${PIPESTATUS[0]}

if [[ $EXIT_CODE -ne 0 ]]; then
    echo "[ERROR] Picard failed with exit code $EXIT_CODE" | tee -a "$LOGFILE"
    exit 1
fi

echo "[READ GROUP] Successfully added read groups." | tee -a "$LOGFILE"
echo "------------------------------------------------------------" | tee -a "$LOGFILE"
