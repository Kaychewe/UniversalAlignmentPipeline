#!/usr/bin/env bash
set -euo pipefail
IN_BAM=${1:? "Input BAM required"}
STATS=${2:? "Stats output required"}
DEPTH=${3:? "Depth output required"}
COV=${4:? "Coverage output required"}
LOG=${5:-"qc.log"}

echo "[LOG] Running samtools stats and depth" | tee -a "$LOG"
samtools stats "${IN_BAM}" > "${STATS}" 2>>"$LOG"
samtools depth -a "${IN_BAM}" > "${DEPTH}" 2>>"$LOG"

awk '
  BEGIN{cov1=0; cov5=0; cov10=0; total=0; sum=0}
  {d=$3; if(d>=1) cov1++; if(d>=5) cov5++; if(d>=10) cov10++; total++; sum+=d}
  END{
    if(total>0){
      printf("≥1x\t%.4f\n≥5x\t%.4f\n≥10x\t%.4f\nMeanDepth\t%.2f\n",
      cov1/total*100, cov5/total*100, cov10/total*100, sum/total)
    } else {
      print "NoDepthData\t0"
    }
  }' "${DEPTH}" > "${COV}"

touch "${IN_BAM%.bam}.qc.ok"
