# JESTR / MIST 冒烟测试报告（CASMI2022）

日期：2026-10-29 · 工作路径：`benchmark/`
目的：确认两个代谢组学模型能在本机（macOS 26 / Apple M4 / arm64，无 NVIDIA GPU）装得上、跑得通，
用 CASMI2022 的 10 张 MS2 谱走完"谱图 → 模型输出"全链路，并记录踩坑点。

---

## 1. 结论

| 模型 | 状态 | 输入 | 输出 | 耗时（10 谱，CPU） | 冒烟指标 |
|---|---|---|---|---|---|
| **JESTR** | ✅ 跑通 | 10 × MS2（50 候选/谱） | 候选排序（rank + cosine） | **0.7 s** | rank-1 **3/10**，rank≤6 **7/10** |
| **MIST** | ✅ 跑通（仅正离子） | 10 × MS2 `[M+H]+` | 4096 位指纹 + 256 维对比嵌入 | **0.19 s**（指纹） | 指纹 Tanimoto **0.242**；Top-1 **6/10** |

两个模型都是"代码 + 预训练权重"型仓库，**不含 CASMI2022 适配**：谱图抽取、候选集构造、`.ms` 转换都需要自己补。
本次补的胶水代码全部在 `benchmark/tools/` 与各自的 `smoke_test/` 下，未改动论文实验代码逻辑
（MIST 仅 1 个文件 4 行 import 保护，见 §5）。

关键结论：**链路通了，但这不是准确率评测**——JESTR 的 decoy 是随机 CASMI 结构，MIST 的候选池只有 10 个真值，
数字只说明"有信号、方向对"。

---

## 2. 冒烟数据（10 张 CASMI2022 谱）

CASMI2022 原始数据是 145 个 DDA `.mzML`，500 个 challenge 挂在 `metadata/casmi2022_challenges.csv`。
每个 challenge 对应 run 里的一张 MS2，需要按 (precursor m/z, RT) 反查扫描。

筛选条件：`[M+H]+`、单电荷、中性质量 150–500 Da、一文件一 challenge、分子式不重复 → 10 张：

| challenge | formula | precursor m/z | RT (min) | 匹配 RT | 峰数 |
|---|---|---|---|---|---|
| 184 | C10H16O2 | 169.1223 | 2.01 | 2.010 | 54 |
| 130 | C6H9N2O2P | 173.0474 | 0.47 | 0.493 | 38 |
| 26 | C12H17NO | 192.1383 | 3.35 | 3.347 | 25 |
| 369 | C11H16O3 | 197.1172 | 1.62 | 1.619 | 78 |
| 163 | C7H6BrNO2 | 215.9655 | 2.05 | 2.049 | 49 |
| 342 | C11H10O5 | 223.0601 | 0.99 | 1.037 | 58 |
| 310 | C12H19NO3 | 226.1438 | 2.64 | 2.666 | 70 |
| 454 | C11H14O5 | 227.0914 | 1.11 | 1.153 | 69 |
| 285 | C14H14N2O | 227.1179 | 4.31 | 4.313 | 25 |
| 219 | C12H20O4 | 229.1434 | 5.45 | 5.457 | 55 |

抽取脚本：`tools/extract_casmi_smoke.py`（pyteomics 逐文件扫描，20 mDa 前体窗 + RT 最近邻），
产出 `tools/smoke10pos/{smoke.mgf, labels.tsv, selected.csv}`。

另有一组 **5 正 + 5 负**（`tools/smoke10/`，`[M+H]+` / `[M-H]-`，质量 <500），用于验证 JESTR 的负离子通路；
MIST 因离子词表限制未跑负离子（§5）。

---

## 3. JESTR

- 仓库：`HassounLab/JESTR1`（MIT），权重随仓附带（NPLIB1，24 MB）。
- 环境：Python 3.12 + torch 2.3.0 + **dgl 2.2.0** + dgllife 0.3.2（`.venv`，922 MB）。
- 流程：MGF → `norm_mzi`（强度归一 999，m/z<1000）→ 1000 bin 分箱（`log10over3`）→ MolEnc/SpecEnc 取嵌入 → 余弦排序。
- 结果（10 谱 × 50 候选）：

| rank-1 | rank≤6 | 最差 | 总耗时 |
|---|---|---|---|
| 3/10 | 7/10 | 32/50（casmi130） | 0.7 s |

- **自检**：用仓库自带 `small_data.mgf` + notebook 候选表复现作者 demo，得 **rank 1/305**
  （true cos 0.956 vs 最优 decoy 0.878）→ 说明适配后模型行为与作者一致，上面的 30% 是数据/候选集问题，不是代码问题。
- 混合正负离子组（5+5）：rank-1 2/10，负离子可正常推理。

详见 [`JESTR/README_SMOKE.md`](JESTR/README_SMOKE.md)。

---

## 4. MIST

- 仓库：`samgoldman97/mist`，权重来自 [Zenodo 8316682](https://zenodo.org/record/8316682)（FP 58 MB + contrastive 11 MB）。
- 环境：Python 3.10（uv 装独立 CPython）+ torch 2.3.0 + **pytorch-lightning 1.8.6**（`.venv`，732 MB）。
- 流程：MGF → 每谱一个 SIRIUS `.ms` → `assign_subformulae.py`（碎片子式）→ `pred_fp.py`（指纹）→ `embed_contrast.py`（对比嵌入）。

| 指标 | 值 |
|---|---|
| 指纹输出 | 10 × 4096 |
| 与真值 Morgan4096 的 Tanimoto（pred>0.5） | 均值 **0.242**（中位 0.233） |
| Top-1 检索（10 真值里选最近） | **6/10** |
| 耗时 | 子式 ≈0.3 s + 指纹 **0.19 s** + 嵌入 0.18 s |

命中 026 / 163 / 310 / 342 / 130 / 369。详见 [`MIST/README_SMOKE.md`](MIST/README_SMOKE.md)。

---

## 5. 踩坑与解决

| # | 问题 | 原因 | 解决 |
|---|---|---|---|
| 1 | conda 建不了环境 | `/opt/anaconda3/envs` 与 `~/.conda/environments.txt` 不可写 | JESTR 用 venv；MIST 用 `uv venv --python 3.10` |
| 2 | DGL 导入失败 | 需写 `~/.dgl`；graphbolt 找不到匹配 torch 的 dylib | `DGLBACKEND=pytorch`；torch 固定 2.3.0（wheel 只到 2.3） |
| 3 | torchdata 导入失败 | `DILL_AVAILABLE` 在 torch≥2.3 被删 | venv 内 `sitecustomize.py` 补属性 |
| 4 | **JESTR 权重维度不匹配（78 vs 74）** | `params.yaml` 的 14 元素表比训练时多 4 个 | 覆盖为 10 元素 `H C O N P S Cl F Br I`，恢复 74 维 |
| 5 | PL 1.9 导入即崩 | 1.9 删了 `loggers.base.rank_zero_experiment` | 降到 **1.8.6** |
| 6 | `pkg_resources` 缺失 | setuptools≥70 移除了它 | 装 `setuptools<70` |
| 7 | MIST 导入要 ray | `utils/__init__` 无条件 `from .tune_utils import *` | 包 `try/except ImportError`（唯一上游改动，4 行） |
| 8 | MIST 找不到谱图 | 只认每谱一个的 `.ms` | `tools/mgf_to_ms.py` 转换 |
| 9 | 负离子 `KeyError: '[M-H]-'` | 预训练离子词表只有 7 种正离子 | 冒烟改用 10 张 `[M+H]+`；负离子需另训模型 |

---

## 6. 复现

```bash
cd benchmark

# 1) 抽 10 张 CASMI2022 谱（一次即可，产物已提交在 tools/smoke10pos/）
python3 tools/extract_casmi_smoke.py --casmi casmi2022 --out-dir tools/smoke10pos --n-pos 10 --n-neg 0

# 2) JESTR（约 1 s）
bash JESTR/smoke_test/run.sh

# 3) MIST（约 1 s）
bash MIST/smoke_test/run.sh
```

两个 `run.sh` 都已在 2026-10-29 04:16 端到端复跑通过。

---

## 7. 产物清单

| 路径 | 内容 |
|---|---|
| `tools/extract_casmi_smoke.py` | mzML → MGF/labels 抽取（按 precursor+RT 匹配扫描） |
| `tools/mgf_to_ms.py` | MGF → SIRIUS `.ms`（MIST 输入格式） |
| `tools/smoke10pos/` | 10 张 `[M+H]+`：`smoke.mgf`、`labels.tsv`、`selected.csv` |
| `tools/smoke10/` | 5 正 + 5 负（JESTR 混合集） |
| `JESTR/smoke_test/run.sh`、`run_jestr_smoke.py` | JESTR 批量化冒烟脚本 |
| `JESTR/smoke_test/out_smoke10pos/jestr_ranks.csv` | 逐谱 rank / cosine |
| `MIST/smoke_test/run.sh` | MIST 四步流水线 |
| `MIST/smoke_test/out/fp_preds_casmi2022pos10.p` | 指纹预测 |
| `MIST/smoke_test/out/contrastive_embeds.p` | 对比嵌入 |
| `MIST/smoke_test/out/retrieval_top1.csv` | Top-1 检索逐条结果 |

---

## 8. 下一步建议

1. **换成"正规"候选集再评一次**：JESTR 用同分子式 PubChem 候选（notebook 里的 `get_from_pubchem` 路径），
   MIST 用 `form_subsets.py` + PubChem 子集建 HDF5 做全库检索——这才能和论文数字对话。
2. **负离子**：MIST 需要带负离子词表的权重（或微调）；JESTR 可直接跑但权重按 `[M+H]+` 训练，负离子结果需谨慎。
3. **放大到全量 CASMI2022**：145 个 mzML / 500 challenge 抽全后跑 500 张谱，JESTR≈0.07 s/谱、MIST≈0.02 s/谱（CPU），
   主要成本在 mzML 解析与子式分配，可并行。
4. **接入 mzDesk benchmark**：把 `run_jestr_smoke.py` / MIST 三步封装成 `benchmark/models/{jestr,mist}.py`，
   与 `std/`（SIRIUS）走同一套 evaluation（Top-k、MRR、指纹 Tanimoto），形成统一对比表。
5. **环境固化**：把两套 `.venv` 的 `pip freeze` 与 `sitecustomize.py`、PL/元素表补丁写进 `environment/`，
   避免下次重装再踩一遍。
