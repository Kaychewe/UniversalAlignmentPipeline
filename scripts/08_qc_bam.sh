#!/usr/bin/env bash

# ------------------------------------------------------------
# 07_qc_bam.sh
# Generate BAM QC metrics: stats, flagstat, idxstats, depth summary
#
# Usage:
#   bash 08_qc_bam.sh <OUTDIR> <LOGDIR> <SAMPLE> <THREADS> <INPUT_BAM>
#
# Outputs:
#   qc/<SAMPLE>.stats.txt
#   qc/<SAMPLE>.flagstat.txt
#   qc/<SAMPLE>.idxstats.txt
#   qc/<SAMPLE>.depth.txt            (optional: full depth)
#   qc/<SAMPLE>.depth_summary.txt    (average depth, % covered)
# ------------------------------------------------------------

set -euo pipefail

if [[ $# -lt 5 ]]; then
    echo "Usage: bash 08_qc_bam.sh <OUTDIR> <LOGDIR> <SAMPLE> <THREADS> <INPUT_BAM>"
    exit 1
fi

OUTDIR=$1
LOGDIR=$2
SAMPLE=$3
THREADS=$4
IN_BAM=$5

QC_DIR="$OUTDIR/$SAMPLE/qc"
SAMPLE_LOGDIR="$LOGDIR/$SAMPLE"
LOGFILE="$SAMPLE_LOGDIR/qc.log"

mkdir -p "$QC_DIR"
mkdir -p "$SAMPLE_LOGDIR"

STATS_OUT="$QC_DIR/${SAMPLE}.stats.txt"
FLAGSTAT_OUT="$QC_DIR/${SAMPLE}.flagstat.txt"
IDXSTATS_OUT="$QC_DIR/${SAMPLE}.idxstats.txt"
DEPTH_OUT="$QC_DIR/${SAMPLE}.depth.txt"       # optional, can disable
DEPTH_SUMMARY="$QC_DIR/${SAMPLE}.depth_summary.txt"

# ---------------------- Conda Init ---------------------------
source /home/kchewe/anaconda3/etc/profile.d/conda.sh
conda activate bissnp_benchmark

# ---------------------- Check samtools -----------------------
if ! command -v samtools &> /dev/null; then
    echo "[ERROR] samtools not found." | tee -a "$LOGFILE"
    exit 1
fi

# ---------------------- Log Header ---------------------------
{
echo "------------------------------------------------------------"
echo " [QC] BAM QC for sample: $SAMPLE"
echo " Timestamp: $(date)"
echo " Input BAM:  $IN_BAM"
echo " Threads:    $THREADS"
echo " Output dir: $QC_DIR"
echo " Log:        $LOGFILE"
echo "------------------------------------------------------------"
} | tee "$LOGFILE"

# ---------------------- 1. samtools stats ---------------------
echo "[QC] Running samtools stats..." | tee -a "$LOGFILE"
samtools stats "$IN_BAM" > "$STATS_OUT" 2>>"$LOGFILE"

# ---------------------- 2. samtools flagstat ------------------
echo "[QC] Running samtools flagstat..." | tee -a "$LOGFILE"
samtools flagstat "$IN_BAM" > "$FLAGSTAT_OUT" 2>>"$LOGFILE"

# ---------------------- 3. samtools idxstats ------------------
echo "[QC] Running samtools idxstats..." | tee -a "$LOGFILE"
samtools idxstats "$IN_BAM" > "$IDXSTATS_OUT" 2>>"$LOGFILE"

# ---------------------- 4. Coverage Summary -------------------
echo "[QC] Computing average depth and coverage…" | tee -a "$LOGFILE"

# Depth summary only (fast, no huge file)
samtools depth -@ "$THREADS" "$IN_BAM" \
    | tee "$QC_DIR/tmp.depth" \
    | awk '
        {sum += $3; if ($3 > 0) covered++} 
        END { 
            print "Average_Depth\t"sum/NR; 
            print "Covered_Positions\t"covered;
            print "Total_Positions\t"NR;
            print "Fraction_Covered\t"covered/NR 
        }
    ' > "$DEPTH_SUMMARY"

mv "$QC_DIR/tmp.depth" "$DEPTH_OUT"

# ---------------------- Finish -------------------------------
{
echo "[QC] Completed QC for $SAMPLE"
echo " Outputs:"
echo "   Stats:           $STATS_OUT"
echo "   Flagstat:        $FLAGSTAT_OUT"
echo "   Idxstats:        $IDXSTATS_OUT"
echo "   Depth:           $DEPTH_OUT"
echo "   Depth Summary:   $DEPTH_SUMMARY"
echo "------------------------------------------------------------"
} | tee -a "$LOGFILE"
