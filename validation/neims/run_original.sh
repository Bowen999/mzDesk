#!/bin/bash
# 用原代码 + 原权重在 Docker 里预测 molecules.csv 和 extra_test_set_smiles.txt，结果写到 original_output.json（体积较大，不入库）。
# 用法：./run_original.sh [原仓库目录] [权重目录]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ORIG="${1:-${NEIMS_ORIG:-/tmp/deep-molecular-massspec}}"
WEIGHTS="${2:-$HOME/.mzdesk/models/neims/massspec_weights}"
N_FP_ONLY="${N_FP_ONLY:-20000}"

if [ ! -d "$ORIG" ]; then
  git clone https://github.com/brain-research/deep-molecular-massspec.git "$ORIG"
  git -C "$ORIG" checkout f82884a076a00fd45185dc303cd67efbc0c26bac     # 验证时用的提交
fi
docker info >/dev/null 2>&1 || { echo "Docker is not running"; exit 1; }
docker build --platform linux/amd64 -t neims-orig "$HERE/docker"

# 指纹对比用的 SMILES：原仓库训练集的前 N 行
head -n "$N_FP_ONLY" "$ORIG/training_splits/train_set_smiles.txt" > "$HERE/fp_smiles.txt"

docker run --rm --platform linux/amd64 \
  -v "$ORIG":/orig:ro -v "$WEIGHTS":/weights:ro -v "$HERE":/work \
  neims-orig python /work/original_predict.py /work/original_output.json /weights "$N_FP_ONLY"
rm -f "$HERE/fp_smiles.txt"
