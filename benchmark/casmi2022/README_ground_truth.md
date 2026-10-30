# CASMI 2022：官方标准答案

`casmi2022_ground_truth.csv` 一道题一行，共 14 列 500 行：Priority / Bonus 各 250。它只包含官方最终答案表直接提供的字段，不包含分析参数、项目计算结果或处理后的谱图。

## 列

| 列 | 内容 |
|---|---|
| `challenge_id` | 官方题号，1–500 |
| `formula`、`adduct_official`、`inchikey`、`smiles_official`、`monoisotopic_mass` | 官方中性分子式、加合物、完整 InChIKey、SMILES 和单同位素质量 |
| `priority`、`file`、`rt_min_given`、`precursor_mz_given` | 官方题目类别、原始文件名、保留时间及给定前体 m/z |
| `collection`、`rack_id`、`well_id`、`mix_or_set` | 官方样本信息 |

这些列分别对应官方 `casmi2022_key_final.xlsx` 中的 14 个字段；整数型样本编号以不带 `.0` 的文本保存。`precursor_mz_given` 是官方给定值，不是项目从 MS1 谱中提取的前体值。

CSV 包含答案，**不可将整行或答案字段送入盲测模型**。它不再是完整 mzDesk 工作流输入；谱图和运行设置需要另外提供。可用 Python `csv.DictReader`、pandas 或 `build_ground_truth.load_records()` 读取。

结构评估按项目基准方案比较 InChIKey 前 14 位；官方 SMILES 中有 277 题缺少官方 InChIKey 的立体信息，不能把 SMILES 衍生的完整 InChIKey 当作完整立体真值。`adduct_official` 保留原文。

## 数据与处理

- 答案直接取自 `metadata/source/casmi2022_key_final.xlsx`，包括两份更正表对应的 4 题；测试逐题核对官方最终答案、更正表和 `metadata/casmi2022_challenges.csv`。
- 项目提取的谱图另存于 `sirius_run/ms/challenge_NNN.ms`，不在此 CSV 中。提取方法见 `sirius_run/scripts/01_extract_ms.py`；源 mzML 属于 MassIVE 数据集 MSV000089239。
- 如需旧的谱图及盲测表单打包逻辑，可使用 `build_ground_truth.build_extracted_records()`，但其输出包含项目处理结果和默认设置，并非纯官方字段；需要 500 个本地 `.ms` 文件。
- 500 题均有 MS2 和 MS1；498 题提取状态为 `OK`，第 18 题使用相邻 MS1 扫描，第 393 题的 RT 窗扩至 ±0.5 min，这两题状态为 `OK_WITH_NOTES`。第 251 题仅有一个 MS1 峰，不能作为同位素模式使用。

## 运行设置与评测注意事项

运行参数不是官方答案的一部分。正式评测应事先统一设定并记录；不同软件应使用同一提取谱和可比容差。

文件不会自动移除训练集或谱库中的同结构，也不提供统一的候选结构库；不能直接据此声称无数据泄漏或公平比较。正式评测应按 `admin/README.md` 先冻结候选库、谱库过滤和训练集分层等规则。

## 重建与检查

```bash
python3 benchmark/casmi2022/build_ground_truth.py
python3 -m pytest -q tests/test_casmi2022_ground_truth.py
```

重建 CSV 需要 `openpyxl` 和 `metadata/source/casmi2022_key_final.xlsx`，不需要本地 `.ms` 文件。
