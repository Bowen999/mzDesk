# 公开 LC-MS/MS 研究

| 研究 | 样本 | 非靶向 / 靶向原始数据 | 处理后数据 | 实测 CCS |
|---|---|---|---|---|
| [Hep-G2 化合物处理（2025）](https://doi.org/10.1007/s11306-025-02309-0) | 276 份 Hep-G2 样本（23 个条件 × 4 个生物学重复 × 3 针）+ 37 个 QC | [MTBLS12332](https://www.ebi.ac.uk/metabolights/MTBLS12332)：4 个 assay 的 `.d.zip`（RP18 / HILIC × 正负离子），共 626 个文件、90.1 GB | 4 个 MAF（747 / 2251 / 3030 / 3880 行），**只有已注释化合物**；`reliability`、`database`、`search_engine` 三列为空 | 有：timsTOF 的 1/K0 在原始数据里，可换算 CCS；MAF 无 CCS 列 |

- 处理软件：MetaboScape 2021b（T-ReX 提取 + 谱库注释），注释要求 MS/MS 证据且 CCS 偏差 ≤ 5%；下游 R 4.4.2 + MetaboAnalyst 6.0。
- **未注释特征**：补充材料只有化合物浓度、LC 梯度和三种富集方法的覆盖 CSV，没有峰表；未注释峰只能从 raw 重建。
- 本地已下载：仅 `RP18_pos` 的 `B9_Comp`（Mevastatin）12 个文件，1.58 GB，位于 `/Volumes/G-TECH/MetabolomicsWorkbench/MTBLS12332/RP18_pos/`。
- B9 板对照为 `B9_DMSO` / `B9_PBS`，未下载。
- 已核对论文、补充材料、MAF 行数、原始文件清单与 zip 完整性；未逐个解析原始文件。
