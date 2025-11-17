#!/usr/bin/env bash
# ------------------------------------------------------------
# Script: build_reference_index.sh
# Purpose: Prepare reference genome indexes for BWA-MEM2,
#          samtools (.fai), Picard (.dict), and Bismark.
# Usage:
#   ./build_reference_index.sh <reference_fasta> <reference_dir> <threads> <log_file>
#
# Example:
#   ./build_reference_index.sh \
#       /home/kchewe/projects/02.BSSNP/genomes/hg38/hg38.fa \
#       /home/kchewe/projects/02.BSSNP/genomes/hg38 \
#       8 \
#       logs/hg38.index.log
# ------------------------------------------------------------


REF_FA=${1:? "Error: reference FASTA path required"}
REF_DIR=${2:? "Error: reference directory required"}
THREADS=${3:-4}
LOG=${4:-"reference_index.log"}

echo "[LOG] Checking reference index for ${REF_FA}" | tee -a "$LOG"
echo "[LOG] Using ${THREADS} threads" | tee -a "$LOG"

# ------------------------------------------------------------
# 1) samtools faidx
# ------------------------------------------------------------
if [ ! -f "${REF_FA}.fai" ]; then
    echo "[LOG] samtools faidx ${REF_FA}" | tee -a "$LOG"
    samtools faidx  "${REF_FA}" -@ $THREADS >>"$LOG" 2>&1
else
    echo "[LOG] FAI present." | tee -a "$LOG"
fi

# ------------------------------------------------------------
# 2) Picard CreateSequenceDictionary
# ------------------------------------------------------------
DICT="${REF_FA%.*}.dict"
if [ ! -f "$DICT" ]; then
    echo "[LOG] Creating sequence dictionary: $DICT" | tee -a "$LOG"
    if [ -n "${PICARD:-}" ]; then
        PIC="java -Xmx32g -jar $PICARD"
    else
        PIC="picard"
    fi
    $PIC CreateSequenceDictionary R="${REF_FA}" O="$DICT" >>"$LOG" 2>&1
else
    echo "[LOG] .dict present." | tee -a "$LOG"
fi

# ------------------------------------------------------------
# 3) BWA-MEM2 index
# ------------------------------------------------------------
if [ ! -e "${REF_FA}.0123" ] && [ ! -e "${REF_FA}.bwt.2bit.64" ]; then
    echo "[LOG] Building BWA-MEM2 index with ${THREADS} threads..." | tee -a "$LOG"
    bwa-mem2 index -t "${THREADS}" "${REF_FA}" >>"$LOG" 2>&1
else
    echo "[LOG] BWA-MEM2 index present." | tee -a "$LOG"
fi

# ------------------------------------------------------------
# 4) Bismark genome preparation
# ------------------------------------------------------------
if [ ! -d "${REF_DIR}/Bisulfite_Genome" ]; then
    echo "[LOG] Preparing Bismark genome (${THREADS} threads): ${REF_DIR}" | tee -a "$LOG"
    bismark_genome_preparation --parallel "${THREADS}" --verbose "${REF_DIR}" >>"$LOG" 2>&1
else
    echo "[LOG] Bismark genome present." | tee -a "$LOG"
fi

# ------------------------------------------------------------
# Done
# ------------------------------------------------------------
echo "[LOG] Reference indexing complete for ${REF_FA}" | tee -a "$LOG"
touch "${REF_DIR}/index.ok"
