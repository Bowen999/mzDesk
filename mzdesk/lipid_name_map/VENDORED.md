# lipid_name_map（拷进来的代码）

- 来源：`~/Desktop/lipid_plus_main/benchmark/name_map/lipid_name_map/`（LIPID+ 项目的名字规范化包）
- 版本：0.4.0（`__version__`）
- 复制日期：2026-10-31
- 作者：Bowen Yang（内部代码）
- 改动：**没有**。除了 `__pycache__` 和 `.pyc`（不拷）以及本说明文件之外，目录和源包逐字节相同，`data/` 里的 json / csv 一并拷入。
  - 包内没有绝对 `import lipid_name_map`（只在文档字符串和命令行说明里出现）；子模块用相对导入，`semantics/registry.py` 用 `__package__` 拼家族模块名；
    数据文件都按 `Path(__file__).parent / "data"` 找，所以放在 `mzdesk.lipid_name_map` 下照常工作。
  - 只用标准库。源包写的是 Python >= 3.12；mzDesk 声明 >= 3.10，包里没有 3.11 / 3.12 才有的语法（用 `ast.parse(feature_version=(3, 10))` 检查过），但只在 3.12 上跑过测试。已在 Python 3.10.19 上验证（审核员在 conda 环境 `fragment` 上跑了 1,500 个名字，结果与 3.12 一致）。
- 用法：`from mzdesk.lipid_name_map import normalize`。Lipid Short Name 页（`mzdesk/api/runners/lipid.py`）和示意结构（`mzdesk/lipid/`）用它。
- 打包：`pyproject.toml` 的 package-data 含 `lipid_name_map/data/**/*`，安装后数据文件在。
- 升级：重新拷源目录覆盖本目录，更新上面的版本和日期，跑 `tests/test_lipid_short_name.py`（版本检查和 golden 对照会提醒，新增的原因码 / flag 要补 `mzdesk/lipid/texts.py`）。
