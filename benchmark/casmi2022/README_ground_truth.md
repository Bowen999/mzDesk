# CASMI 2022：输入与标准答案

`casmi2022_ground_truth.jsonl` 是一份自包含的数据文件，每行一题，共 500 题：Priority / Bonus 各 250，正离子 304、负离子 196。它记录的是**官方标准答案**，不是 mzDesk 的预测结果。

## 每行内容

| 字段 | 内容 |
|---|---|
| `challenge_id` | 官方题号，1–500 |
| `inputs` | 可直接交给 mzDesk `esi_ms2_annotation` 工作流的完整表单；`precursor.adduct` 为 `Auto`，不泄露真值 |
| `ground_truth` | 官方中性分子式、加合物、SMILES、完整 InChIKey、InChIKey 前 14 位和单同位素质量；另附理论 m/z 与校验备注 |
| `provenance` | 官方原始前体 m/z 和 RT、样本信息、mzML 路径、提取状态、扫描 RT、提取谱文件的 SHA-256 |
| `ms2_by_collision_energy` | 合并前各碰撞能量的 MS2 峰表，`peaks` 每项为 `[m/z, intensity]` |

注意：文件同时包含输入和答案。盲测时**只传 `inputs` 对象**，不要把整条记录或 `ground_truth` 送入模型。例如：

```python
import json
from mzdesk.annotation.pipeline import parse_inputs

with open("benchmark/casmi2022/casmi2022_ground_truth.jsonl", encoding="utf-8") as f:
    record = json.loads(next(f))
query = parse_inputs(record["inputs"])  # 只验证表单，不启动分析
```

结构评估按项目基准方案比较 InChIKey 前 14 位；官方 SMILES 中有 277 题缺少官方 InChIKey 的立体信息，不能把 SMILES 衍生的完整 InChIKey 当作完整立体真值。`adduct_official` 保留原文；`adduct_mzdesk` 将 `[M+FA-H]-` 统一为 `[M+HCOO]-`。

## 数据与处理

- 答案取自 `metadata/source/casmi2022_key_final.xlsx`，包括两份更正表对应的 4 题；项目现有 `metadata/casmi2022_challenges.csv` 和 `sirius_run/results/truth.tsv` 与最终答案一致。生成时再次逐题核对这两个表及提取报告。
- 谱图取自项目已有的盲法提取结果 `sirius_run/ms/challenge_NNN.ms`，提取方法见 `sirius_run/scripts/01_extract_ms.py`：从原始 mzML 按题目 m/z 和 RT 查找扫描，保留 MS1 同位素峰及最多 3 个碰撞能量的 MS2。源 mzML 属于 MassIVE 数据集 MSV000089239，路径保留在记录中。
- mzDesk 一次接收一张 MS2，`inputs.ms2` 因此使用项目 `benchmark/tools/jm/prep.py` 的已有规则：按 m/z 排序，距簇中首峰不超过 10 ppm 的峰分为一簇，保留簇中最强峰的 m/z 和强度，不对强度求和或跨能量归一化。合并前各能量峰表仍完整保留。此谱与现有 JESTR 合并谱逐题核对一致；原始提取输出的 m/z 和强度分别保留 6 位和 2 位小数。
- `inputs.precursor.precursor_mz` 是提取谱的 MS1 单同位素实测值，保留 5 位小数；官方给定值另存于 `provenance.precursor_mz_given`，两者不可混淆。RT 仅用于溯源，不是当前工作流表单输入。
- 500 题均有 MS2 和 MS1；498 题提取状态为 `OK`，第 18 题使用相邻 MS1 扫描，第 393 题的 RT 窗扩至 ±0.5 min，这两题状态为 `OK_WITH_NOTES`。第 251 题仅有一个 MS1 峰，不能作为同位素模式使用。

## 运行设置与评测注意事项

`inputs` 使用当前页面的默认示例设置：MassBank 和 GNPS、10 ppm 前体容差、0.02 Da 碎片容差、组学类型未知、`halogen=false`、不额外比较数据库同分异构体、流动相添加剂和碱金属未知、pH 2.7。**这些是表单设置，不是 CASMI 实验条件或冻结的基准方案**。尤其 `halogen=false` 会限制分子式搜索空间；正式评测应事先统一设定并记录。不同软件应使用同一提取谱和可比容差。

文件不会自动移除训练集或谱库中的同结构，也不提供统一的候选结构库；不能直接据此声称无数据泄漏或公平比较。正式评测应按 `admin/README.md` 先冻结候选库、谱库过滤和训练集分层等规则。

## 重建与检查

```bash
python3 benchmark/casmi2022/build_ground_truth.py
python3 -m pytest -q tests/test_casmi2022_ground_truth.py
```

重建需要 `sirius_run/ms/` 中 500 个本地 `.ms` 文件；它们被 git 忽略。若缺失，需先取得 `metadata/casmi2022_challenges.csv` 所列原始 mzML，并用 `sirius_run/scripts/01_extract_ms.py` 重新提取（依赖 pyteomics、numpy 和 pandas）。JSONL 本身已包含分析所需峰表，使用时不需要这些本地文件。
