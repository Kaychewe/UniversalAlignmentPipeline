#!/usr/bin/env bash
# ------------------------------------------------------------
# Script: sort_bam.sh
# Purpose: Sort and index BAM safely for Snakemake pipelines.
# Usage:   ./sort_bam.sh <input.bam> <output.bam> <threads> <logfile>
# ------------------------------------------------------------

inbam=$1
outbam=$2
threads=${3:-8}
log=${4:-sort.log}

echo "[LOG] Sorting BAM..." | tee -a "$log"

if [ ! -s "$inbam" ]; then
    echo "[ERROR] Input BAM does not exist or is empty: $inbam" | tee -a "$log"
    exit 1
fi

# Sort BAM
samtools sort -@ "$threads" -o "$outbam" "$inbam" >>"$log" 2>&1
if [ $? -ne 0 ] || [ ! -s "$outbam" ]; then
    echo "[ERROR] Sorting failed or output BAM missing: $outbam" | tee -a "$log"
    exit 1
fi

# Index BAM
samtools index "$outbam" >>"$log" 2>&1
if [ $? -ne 0 ]; then
    echo "[WARN] BAM indexing encountered issues, continuing..." | tee -a "$log"
fi

# Create sentinel OK file
touch "${outbam%.bam}.ok"
echo "[LOG] Sorting completed successfully: $outbam" | tee -a "$log"
