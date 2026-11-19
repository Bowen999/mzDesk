# NEIMS（EI 质谱预测）验证：numpy 重写 vs 原代码

`mzdesk/ms2/neims_engine.py` 用 numpy + RDKit 重写了 NEIMS 的推理。这里把它和**原代码 + 原权重**（Docker 里的 TF 1.13）在同一批分子上逐项对照。

- 原代码：<https://github.com/brain-research/deep-molecular-massspec>，提交 `f82884a`（Apache-2.0，已归档）
- 原权重：官方 `massspec_weights.zip`（Wayback Machine 存档，sha256 `3c80ee1b…aa52`），检查点 `model.ckpt-100000`
- 论文：Wei et al., ACS Cent. Sci. 2019, 5, 700-708

## 复现

```bash
python -m mzdesk.ms2.neims_setup          # 下载权重并转成 ~/.mzdesk/models/neims/neims_weights.npz（需要一次 TensorFlow）
validation/neims/run_original.sh          # Docker（linux/amd64）里跑原代码，写 original_output.json（约 10 MB，不入库）
python validation/neims/compare.py        # 对比，写 comparison.csv、summary.json，并更新 tests/fixtures/neims_reference.json
```

原代码环境（`docker/Dockerfile`）：Python 3.6.15、TensorFlow 1.13.2、RDKit **2019.03.4**、numpy 1.16.6。宿主机：RDKit 2024.09.6、numpy 1.26。

## 原代码的做法（numpy 版照此实现）

| 项 | 内容 |
|---|---|
| 特征 | `GetHashedMorganFingerprint(mol, 2, nBits=4096)`：Morgan **计数**指纹，半径 2，4096 位；无其他特征 |
| 分子量 | `CalcExactMolWt`（精确质量）→ float32 → 四舍五入（偶数优先），记作 `model_mass`，用作反向预测的锚点 |
| 网络 | Dense 2000 + ReLU；7 个残差块 `h += Dense2000(ReLU(BN(Dense1000(ReLU(BN(h))))))`；BN + ReLU；三个线性头（1000 维）：正向、反向、门控 |
| 合并 | 正向谱只保留 m/z <= M+5；反向谱 `rev[j] = bwd[M-j+5]`（越界为 0）；`out = ReLU(sigmoid(gate)·fwd + (1-sigmoid(gate))·rev)` |
| 输出 | m/z 0-999 的 1000 个 bin；除以最大值后乘 999 |
| 写出 | 所有非零 bin 以 `'%d %d'` 写入 SDF：强度被**截断**为整数，没有阈值（强度 < 1 的 bin 会写成 0） |
| 超参数 | 来自 `spectra_predictor._DEFAULT_HPARAMS`；`.meta` 图确认：18 个 Dense、15 个 BN（eps 1e-3）、1 个 Sigmoid、2 个 Round、偏移常数 5 |

权重 96 个张量（42.2 M 参数，float32，npz 169 MB），优化器的 Adam 状态不转换。

## 引擎的取舍

- **显示阈值** `MIN_REL = 0.1 %`：原代码没有阈值；0.1 % 约等于原代码整数化后不为 0 的最小强度（1/999）。ReLU 输出本来就稀疏，这个阈值在 35 个分子里每个只去掉 0-4 个峰。
- **标称质量** `molecular_ion_mz`：同位素最丰同位素的整数质量之和（M+• 真正出现的位置）。模型内部用的 `model_mass`（精确质量四舍五入）在大的富氢分子上会多 1（例：三硬脂酸甘油酯 890 vs 891），多溴化物上会少 1（C12H3Br7O：716 vs 715）；测试集 11600 个分子里有 26 个不同，25 个是 C32 以上的大分子（多为长链脂类）。
- 盐 / 多组分：取最大组分并给警告（EI 测的是中性分子）；带净电荷：警告；元素不在训练集常见元素内（训练集 260306 个分子里出现不足 100 次）：警告；重原子 > 100（原代码训练集的上限，实际最大 91）：警告；标称质量 > 999：报错。

## 结果

分子：`molecules.csv` 的 35 个（芳香、脂肪、卤代、含 N / S / P / Si / Se、咖啡因、甾体、糖、M = 890 的三硬脂酸甘油酯等）+ `extra_test_set_smiles.txt`（原仓库测试集随机抽 300 个，seed 20241024）。谱的比较用 999 刻度（最大峰 = 999）。

| | 35 个指定分子 | 300 个测试集分子 |
|---|---|---|
| 指纹完全一致 | 35 / 35 | 299 / 300 |
| `model_mass` 一致 | 35 / 35 | 300 / 300 |
| 谱的最大绝对差（指纹一致的分子） | 0.0013 | 0.0016 |
| 原代码 SDF 峰文本逐字相同 | 35 / 35 | 297 / 300 |
| 显示峰表（>= 0.1 %，1 位小数）相同 | 35 / 35 | 297 / 300 |

- 指纹一致的分子，谱只差 float32 数值噪声（最大 0.0016 / 999）。其中 4 个分子的峰表不完全相同：2 个是整数截断刚好落在边界上（某个峰差 1），2 个是 1 位小数四舍五入落在边界上（差 0.1），没有任何峰只出现在一边。
- 另 1 个是指纹不同（见下），那个分子的谱完全不同（`test_set_038`）。
- 把原代码的指纹和质量直接喂给 numpy 前向，最大差仍是 0.0016，说明网络部分本身是等价的。

### RDKit 版本带来的差异

对原仓库训练集的前 20000 个分子，只比较指纹：容器（RDKit 2019.03.4）和宿主机（2024.09.6）有 **10 个**不同（0.05 %），全是大的稠合多环 / 杂原子桥联体系，原因是新版 RDKit 的芳香性判定不同（芳香原子数不同），另有一个含 Rh-I 的金属配合物。它们的 SMILES 在 `summary.json`。这些分子在新版 RDKit 下的预测会和原代码不同；常见分子不受影响。

### 咖啡因（只报告，没有调参）

预测：M+• 194 是基峰；其余较强的峰 m/z（相对强度 %）：67 (73.6)、193 (66.0)、82 (65.6)、109 (55.9)、68 (47.6)、195 (44.9)、56 (44.1)、81 (43.7)、83 (41.1)、42 (38.8)、137 (38.4)、55 (34.1)、165 (32.7)、180 (31.5)。这些数与 Docker 里的原代码一致。

实测咖啡因 EI 谱（凭记忆，未在本机核对）也是 194 为基峰，109、82、67、55 为主要碎片：这几个峰在预测里都出现了，且都 >= 30 %。但 193 和 195 在预测里很强，实测谱里它们是小峰，这是模型本身的特点（和原代码一致），不是移植误差。

### 速度

本机（arm64 Mac）：导入 0.25 s，权重载入 0.04 s（文件缓存热时），页面 runner 首次调用（含权重载入）约 0.35 s，之后每次预测 8-17 ms；进程峰值内存约 250 MB。

## 文件

| 文件 | 内容 |
|---|---|
| `molecules.csv`、`extra_test_set_smiles.txt` | 输入分子 |
| `docker/Dockerfile`、`run_original.sh`、`original_predict.py` | 在 Docker 里运行原代码 |
| `compare.py` | 对比脚本；还会写 `tests/fixtures/neims_reference.json`（5 个分子，来自原代码的指纹、质量和谱） |
| `comparison.csv`、`summary.json` | 对比结果（每个分子一行 / 汇总） |
