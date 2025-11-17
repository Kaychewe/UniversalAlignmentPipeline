#!/usr/bin/env bash
set -euo pipefail
IN_BAM=${1:? "Input BAM required"}
OUT_BAM=${2:? "Output BAM required"}
SAMPLE=${3:? "Sample name required"}
LOG=${4:-"addRG.log"}

echo "[LOG] Adding read groups for ${SAMPLE}" | tee -a "$LOG"
if [ -n "${PICARD:-}" ]; then
    PIC="java -Xmx8g -jar $PICARD"
else
    PIC="picard"
fi

$PIC AddOrReplaceReadGroups \
    I="${IN_BAM}" \
    O="${OUT_BAM}" \
    ID="${SAMPLE}" \
    LB="lib_${SAMPLE}" \
    PL="ILLUMINA" \
    PU="unit1" \
    SM="${SAMPLE}" \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT >>"$LOG" 2>&1

touch "${OUT_BAM%.bam}.rg.ok"
