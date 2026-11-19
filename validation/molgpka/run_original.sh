#!/bin/bash
# 用原代码 + 原权重（PyTorch + torch_geometric）预测 molecules.csv 和 extra_test_set_smiles.txt，结果写到 original_output.json（不入库）。
# 用一个临时 venv（--system-site-packages 借用已装的 torch / pandas / RDKit，只把 torch_geometric 装进 venv），不动基础环境。
# 用法：./run_original.sh [工作目录]       环境变量 PY 指定装有 torch、pandas、RDKit 的 Python（默认 python3）
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="${1:-${MOLGPKA_WORK:-/tmp/molgpka_orig}}"
PY="${PY:-python3}"
COMMIT=4dc835219baad18c2e705274ac54d49597703966          # 验证时用的提交（和 pka_setup.py 下载权重的提交一致）

mkdir -p "$WORK"
if [ ! -d "$WORK/repo" ]; then
  git clone https://github.com/Xundrug/MolGpKa.git "$WORK/repo"
  git -C "$WORK/repo" checkout "$COMMIT"
fi
if [ ! -x "$WORK/venv/bin/python" ]; then
  "$PY" -m venv --system-site-packages "$WORK/venv"
  "$WORK/venv/bin/python" -m pip install --quiet "torch_geometric==2.4.0"   # 2.4 仍带原代码用的 GlobalAttention 和 data.DataLoader
fi
MOLGPKA_COMMIT="$COMMIT" "$WORK/venv/bin/python" -I "$HERE/run_original.py" "$WORK/repo" "$HERE/original_output.json"
