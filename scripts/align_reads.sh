#!/usr/bin/env bash
# ------------------------------------------------------------
# Script: align_reads.sh
# Purpose: End-to-end alignment pipeline supporting both
#          single-end and paired-end reads using either
#          Bismark (WGBS) or BWA-MEM2 (WGS).
#
# Features:
#   • Auto-detects SE vs PE
#   • Runs FastQC before and after trimming (unless --skip-qc)
#   • Adapter & quality trimming via fastp (unless --skip-qc)
#   • Clean logging and error handling
#
# Usage:
#   ./align_reads.sh [--skip-qc] <aligner> <fastq_R1> <fastq_R2|none> \
#       <reference_fasta> <reference_dir> <threads> \
#       <output_bam> <log_file>
#
# Example:
#   ./align_reads.sh bwa-mem2 \
#       data/A549_R1.fastq.gz data/A549_R2.fastq.gz \
#       /home/kchewe/projects/02.BSSNP/genomes/hg38/hg38.fa \
#       /home/kchewe/projects/02.BSSNP/genomes/hg38 \
#       16 \
#       results/A549.unsorted.bam \
#       logs/A549.align.log
# ------------------------------------------------------------

set -euo pipefail # might need to be removed

# ------------------------------------------------------------
# Parse --skip-qc flag
# ------------------------------------------------------------
SKIP_QC=false
for arg in "$@"; do
    if [[ "$arg" == "--skip-qc" ]]; then
        SKIP_QC=true
        shift
        break
    fi
done

# ------------------------------------------------------------
# Positional arguments
# ------------------------------------------------------------
ALIGNER=${1:? "Error: aligner (bismark|bwa-mem2) required"}
FASTQ_R1=${2:? "Error: FASTQ R1 required"}
FASTQ_R2=${3:-"none"}
REF_FA=${4:? "Error: reference FASTA required"}
REF_DIR=${5:? "Error: reference directory required"}
THREADS=${6:-8}
OUT_BAM=${7:? "Error: output BAM required"}
LOG=${8:-"alignment.log"}

OUT_DIR=$(dirname "${OUT_BAM}")
QC_DIR="${OUT_DIR}/qc"
TRIM_DIR="${OUT_DIR}/trimmed"
mkdir -p "${OUT_DIR}" "${QC_DIR}" "${TRIM_DIR}"

echo "====================================================" | tee -a "$LOG"
echo "[LOG] Starting alignment using ${ALIGNER}" | tee -a "$LOG"
echo "[LOG] FASTQ R1: ${FASTQ_R1}" | tee -a "$LOG"
echo "[LOG] FASTQ R2: ${FASTQ_R2}" | tee -a "$LOG"
echo "[LOG] Reference: ${REF_FA}" | tee -a "$LOG"
echo "[LOG] Threads: ${THREADS}" | tee -a "$LOG"
echo "[LOG] Skip QC: ${SKIP_QC}" | tee -a "$LOG"
echo "====================================================" | tee -a "$LOG"


# ------------------------------------------------------------
# Skip rerun if BAM already exists
# ------------------------------------------------------------
if [ -s "${OUT_BAM}" ]; then
    echo "[LOG] Output BAM already exists (${OUT_BAM}). Skipping alignment." | tee -a "$LOG"
    exit 0
fi


# ------------------------------------------------------------
# Auto-detect paired vs single-end
# ------------------------------------------------------------
if [ "${FASTQ_R2}" = "none" ] || [ ! -s "${FASTQ_R2}" ]; then
    MODE="SE"
    echo "[LOG] Detected single-end mode." | tee -a "$LOG"
else
    MODE="PE"
    echo "[LOG] Detected paired-end mode." | tee -a "$LOG"
    
    # Extract base names
    BASENAME_R1=$(basename "${FASTQ_R1}" .fastq.gz)
    BASENAME_R2=$(basename "${FASTQ_R2}" .fastq.gz)
fi

# ------------------------------------------------------------
# Optional QC + trimming steps
    # TODOs
    # add check if trimmed fastq already exists 
    # check integrity 
    # re-run only if corrput or truncated 
# ------------------------------------------------------------
if [ "$SKIP_QC" = false ]; then

    echo "[LOG] Running QC and trimming (disabled by default)" | tee -a "$LOG"
    if [ "$MODE" = "PE" ]; then
        fastqc -t "${THREADS}" -o "${QC_DIR}" "${FASTQ_R1}" "${FASTQ_R2}" >>"$LOG" 2>&1
    else
        fastqc -t "${THREADS}" -o "${QC_DIR}" "${FASTQ_R1}" >>"$LOG" 2>&1
    fi

    echo "[LOG] Running fastp for adapter/quality trimming..." | tee -a "$LOG"
    if [ "$MODE" = "PE" ]; then


        # Paired-end trimming
        fastp -i "${FASTQ_R1}" -I "${FASTQ_R2}" \
            -o "${TRIM_DIR}/${BASENAME_R1}.trimmed.fq.gz" \
            -O "${TRIM_DIR}/${BASENAME_R2}.trimmed.fq.gz" \
            -w "${THREADS}" \
            -h "${QC_DIR}/fastp.html" -j "${QC_DIR}/fastp.json" >>"$LOG" 2>&1

        CLEAN_R1="${TRIM_DIR}/${BASENAME_R1}.trimmed.fq.gz"
        CLEAN_R2="${TRIM_DIR}/${BASENAME_R2}.trimmed.fq.gz"
    else
        BASENAME_R1=$(basename "${FASTQ_R1}" .fastq.gz)
        fastp -i "${FASTQ_R1}" \
            -o "${TRIM_DIR}/${BASENAME_R1}.trimmed.fq.gz" \
            -w "${THREADS}" \
            -h "${QC_DIR}/fastp.html" -j "${QC_DIR}/fastp.json" >>"$LOG" 2>&1

        CLEAN_R1="${TRIM_DIR}/${BASENAME_R1}.trimmed.fq.gz"
        CLEAN_R2="none"
    fi

    echo "[LOG] Running FastQC on cleaned reads..." | tee -a "$LOG"
    if [ "$MODE" = "PE" ]; then
        fastqc -t "${THREADS}" -o "${QC_DIR}" "${CLEAN_R1}" "${CLEAN_R2}" >>"$LOG" 2>&1
    else
        fastqc -t "${THREADS}" -o "${QC_DIR}" "${CLEAN_R1}" >>"$LOG" 2>&1
    fi
else
    echo "[LOG] --skip-qc flag detected: skipping FastQC + fastp" | tee -a "$LOG"
    CLEAN_R1="${FASTQ_R1}"
    CLEAN_R2="${FASTQ_R2}"
    # CLEAN_R1="${TRIM_DIR}/${BASENAME_R1}.trimmed.fq.gz"
    # CLEAN_R2="${TRIM_DIR}/${BASENAME_R2}.trimmed.fq.gz"
fi

# ------------------------------------------------------------
# Alignment step
# ------------------------------------------------------------
echo "[LOG] Beginning alignment (${ALIGNER}, ${MODE})..." | tee -a "$LOG"

if [ "${ALIGNER}" = "bismark" ]; then
    # Cap parallel threads to avoid excessive spawning
    if [ "$THREADS" -gt 30 ]; then
        PARALLEL=20
    else
        PARALLEL=$THREADS
    fi

    # Run Bismark alignment
    if [ "$MODE" = "PE" ]; then
        bismark --genome "${REF_DIR}" \
                -1 "${FASTQ_R1}" -2 "${FASTQ_R2}" \
                --output_dir "${OUT_DIR}" \
                --parallel "${PARALLEL}" >>"$LOG" 2>&1
        BAM_FILE=$(find "${OUT_DIR}" -maxdepth 1 -type f -name "*_pe.bam" | head -n 1)
    else
        bismark --genome "${REF_DIR}" \
                "${FASTQ_R1}" \
                --output_dir "${OUT_DIR}" \
                --parallel "${PARALLEL}" >>"$LOG" 2>&1
        BAM_FILE=$(find "${OUT_DIR}" -maxdepth 1 -type f -name "*_se.bam" | head -n 1)
    fi

    # Rename unpredictable BAM filename to Snakemake-expected output
    if [ -s "${BAM_FILE}" ]; then
        mv -v "${BAM_FILE}" "${OUT_BAM}" | tee -a "$LOG"
        echo "[LOG] Alignment complete (renamed to expected output): ${g}" | tee -a "$LOG"
    else
        echo "[ERROR] Bismark produced no BAM output! Check ${LOG}" | tee -a "$LOG"
        exit 1
    fi

elif [ "${ALIGNER}" = "bwa-mem2" ]; then
    if [ "$MODE" = "PE" ]; then
        bwa-mem2 mem -t "${THREADS}" "${REF_FA}" "${FASTQ_R1}" "${FASTQ_R2}" \
            | samtools view -b - > "${OUT_BAM}" 2>>"$LOG"
    else
        bwa-mem2 mem -t "${THREADS}" "${REF_FA}" "${FASTQ_R1}" \
            | samtools view -b - > "${OUT_BAM}" 2>>"$LOG"
    fi

    if [ -s "${OUT_BAM}" ]; then
        echo "[LOG] Alignment complete: ${OUT_BAM}" | tee -a "$LOG"
    else
        echo "[ERROR] Alignment failed: BAM is empty." | tee -a "$LOG"
        exit 1
    fi
else
    echo "[ERROR] Unknown aligner specified: ${ALIGNER}" | tee -a "$LOG"
    exit 1
fi

echo "[LOG] Done." | tee -a "$LOG"