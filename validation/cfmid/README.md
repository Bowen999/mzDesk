> **归档说明**：这是最初的 CFM-ID Docker 包装脚本（原路径 `~/Documents/software/CFM-ID`），仅作验证记录保留。
> 应用内的 ESI MS/MS 页面现在使用 `mzdesk/ms2/cfmid_engine.py`（同一镜像，同样先 MSRB、再回退 cfm-predict）。
> `test_cfmid.sh` 仍可在本目录直接运行（路径取自脚本所在目录），输出在 `test_output/`。

# CFM-ID Docker 预测脚本

## 前置条件
- Docker Desktop(脚本会自动启动并拉取镜像 `wishartlab/cfmid:msml4.4.10_msrb1.1.15`)

## 用法
1. 编辑 `input.csv`,四列(表头必填):`id, smiles, adduct_type, note`
   - `note` 可为空;含逗号的字段用双引号括起(如 `"PI(18:1/11,12-EpETE)"`)
2. 运行:

```bash
./test_cfmid.sh          # 默认:MSRB 规则碎片引擎(推荐,脂质专用)
./test_cfmid.sh msrb     # 同上
./test_cfmid.sh cfm      # 通用 ML 引擎 cfm-predict
```

- **MSRB**:覆盖约 20 类脂质及多酚等,支持 `[M+H]+`、`[M-H]-`、`[M+Na]+`、`[M+NH4]+`、`[M+Cl]-` 等加合物;不支持的化合物自动回退到 cfm-predict
- **cfm-predict**:仅支持 `[M+H]+` / `[M-H]-`,其余 adduct 报错

## 输出
- 全部结果合并为一个 MSP 文件:`test_output/output.msp`(每条记录含 `Adduct:` 与 `Note:` 标注,每个化合物 3 个碰撞能量)
- 中间文件在 `test_output/msrb/`(单化合物原始谱与运行日志)
- 退出码 0 = 全部通过;1 = 有失败
