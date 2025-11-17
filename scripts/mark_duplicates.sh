#!/usr/bin/env bash
set -euo pipefail
IN_BAM=${1:? "Input BAM required"}
OUT_BAM=${2:? "Output BAM required"}
METRICS=${3:? "Metrics file required"}
LOG=${4:-"dedup.log"}

echo "[LOG] Marking duplicates" | tee -a "$LOG"
if [ -n "${PICARD:-}" ]; then
    PIC="java -Xmx12g -jar $PICARD"
else
    PIC="picard"
fi

$PIC MarkDuplicates \
    I="${IN_BAM}" \
    O="${OUT_BAM}" \
    M="${METRICS}" \
    REMOVE_DUPLICATES=false \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT >>"$LOG" 2>&1

touch "${OUT_BAM%.bam}.dedup.ok"
