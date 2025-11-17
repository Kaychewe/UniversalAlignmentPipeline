#!/usr/bin/env bash
set -euo pipefail

IN=$1
OUT=$2

source /home/kchewe/anaconda3/etc/profile.d/conda.sh
conda activate bissnp_benchmark

samtools flagstat "$IN" > "$OUT"
