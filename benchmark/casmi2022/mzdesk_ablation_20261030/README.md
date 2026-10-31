# CASMI 2022：mzDesk 全部 500 题——一次重型运行 + 离线重放（完整流程、消融、去除泄漏）

目的：ESI MS/MS Annotation 的第 1–5 步（分子式、谱库、MS-FINDER、FRIGID）只在 500 题上各跑一遍，原始输出存盘；
第 6–8 步（合并、重排、同分异构体、答案）由流程自己的 `run_query` 在存下的输出上重放，因此各种配置（完整流程、不搜库的无泄漏流程、各层消融）不需要再跑引擎，
只有第 7 步（对新的 top-1 再跑一次 MS-FINDER）需要真实运行，并按参数缓存。

## 文件

| 文件 | 作用 |
|---|---|
| `common.py` | 冻结代码根、盲测输入 `prepare()`、跨进程文件锁、并行/可续跑/日志的进程池 |
| `collect.py` | 重型一遍：500 题的第 1–5 步无条件运行（MS-FINDER 和 FRIGID 即使流程会跳过也运行；FRIGID 只对正离子加合物有效，负离子由引擎自己记为跳过），谱库候选不截断。不读答案 |
| `replay.py` | 10 个配置的重放，写 `replay/<配置>/challenge_NNN.json`；第 7 步缓存在 `step7_cache/`。不读答案 |
| `verify.py` | 抽查：重放 == 正常 `run_query`；消融里省略第 7 步不改排名；full 的第 7 步未失败；leakfree_full 无谱库输出 |
| `standalone.py` | 单软件默认设置（谱库、BUDDY、MS-FINDER、FRIGID 各自单独运行，不经 mzDesk 的加合物/分子式/合并步骤），写 `standalone/<工具>/challenge_NNN.json`；`--verify-frigid` 核对 FRIGID 复用 |
| `buddy_record.py` | 做出 BUDDY 缓存（`cNNN.pkl`：每题每个加合物的 BUDDY 输出），`standalone.py` 的 BUDDY 默认行与 FRIGID 默认分子式由它而来（缓存路径 `MZDESK_BUDDY_CACHE`） |
| `score.py` | 逐题复算排名（含对照工具，并与已有评分表核对）→ `scores.json`（机器可读的逐题结果） |
| `write_report.py` | `scores.json` → `benchmark/result/CASMI2022_mzdesk_full500.md`（只含结果表） |
| `leak_lists/` | 各工具"已见"结构名单（IK14），复现报告 (c) 所需；来源脚本在 `leak_lists/sources/` 与 `overlap_source.py` |

## 配置（replay.py）

- `full`：流程原样，闸门/跳过判断由 `run_query` 决定，含第 7 步。`leakfree_full`（mzDesk 无泄漏）：同上，但完全不搜谱库——第 3 步无输出、闸门不通过，直接从头注释；第 4–5 步按流程自己的判断（MS-FINDER；MS-FINDER 不够好且为正离子时才跑 FRIGID），再合并、重排、第 7 步、答案。不读真值。
  谱库层没有任何候选的题（存盘的谱库步骤无错误且 0 个候选）上，去掉谱库不改变任何输入，`leakfree_full` 直接取 `full` 的结果（`copied_from_full: true`）；已用真实重放核对：两种做法都跑过的此类题，排名、候选、答案与第 4–7 步完全一致。有谱库候选的题才重放。
- 消融（只看排名；第 7 步换成"跳过"的桩，各层一律并入，不论闸门/"够好"）：`all_layers`（参照）、`library_only`、`msfinder_only`、`frigid_only`（仅正离子）、
  `no_library`、`no_msfinder`、`no_frigid`。这些消融都是强制并入各层，与 `leakfree_full` 不同。
- `buddy_only`：分子式层，报告真值分子式在第 2 步分子式列表中的位置。
- 第 7 步只给出同分异构体判断，不改变任何候选的名次或级别，所以消融里省略它是干净的（`verify.py` 逐题核对）；代价：消融结果里的 `answer.type` / `isomer` 不可用。

## 命令与顺序

先冻结代码（不要用正在被其他会话修改的工作区）：

```bash
cd /path/to/mzDesk
rm -rf /tmp/mzdesk_code_v0.2.40 && mkdir -p /tmp/mzdesk_code_v0.2.40
git archive <v0.2.40 提交> | tar -x -C /tmp/mzdesk_code_v0.2.40
export MZDESK_CODE_ROOT=/tmp/mzdesk_code_v0.2.40
D=benchmark/casmi2022/mzdesk_ablation_20261030
```

```bash
python3 $D/collect.py --workers 3 --heavy-slots 2          # 1. 重型一遍，可中断后重启续跑
python3 $D/replay.py  --workers 3 --heavy-slots 2          # 2. 10 个配置的重放（collect 全部完成后）
python3 $D/verify.py  --ids 9,16,116,227,370               # 3. 抽查（每题另跑一次正常流程，数分钟）
python3 $D/score.py                                        # 4. -> scores.json（已有则加 --replace）
python3 $D/write_report.py                                 # 5. -> benchmark/result/CASMI2022_mzdesk_full500.md（已有则加 --replace）
```

单软件默认设置（第 4 步之前；需要 BUDDY 缓存，`MZDESK_FRIGID_DEVICE=cpu`）：

```bash
export MZDESK_BUDDY_CACHE=/path/to/ppm/cache                # buddy_record.py 的输出（cNNN.pkl）
export MZDESK_FRIGID_DEVICE=cpu
export MZDESK_MSFINDER_DIR=/path/to/intact/MSFINDER_v3.73   # 本次用的是 ~/tools/msfinder 的副本：原安装里 MsfinderFormulaDB-VS13_bin.efd 被截成 0 字节，MS-FINDER 无法运行；副本里换回发行包原件
python3 $D/standalone.py --verify-frigid 16,10,15          # a. 重跑 FRIGID，与 collect 存盘的输出逐项比较（全部一致才复用存盘输出）
python3 $D/standalone.py --workers 3 --heavy-slots 2       # b. 全部：谱库 500 题；BUDDY / MS-FINDER / FRIGID 只在真实加合物 = 默认加合物的 285 题（FRIGID 其中 177 题正离子）
python3 $D/standalone.py --tools library,buddy             #    只跑不需要重型引擎的部分（秒级）
```

- 默认加合物 [M+H]+ / [M-H]-；真实加合物是别的 215 题不运行，由 `score.py` 直接判错（留在分母中）。真值在 `standalone.py` 里只用来选这 285 题，`--all-adducts` 则全部运行。
- 谱库默认：前体 10 ppm、加权谱熵 ≥ 0.70、匹配峰 ≥ min(4, 峰数)、按 IK14 去重、按相似度排序，无闸门检查，不截断。
- MS-FINDER 默认：`msfinder_step.run_msfinder_layer` 不给分子式列表（MS-FINDER 自己找分子式），其余同 mzDesk 的 MS-FINDER 层。FRIGID 默认：分子式 = BUDDY 在 [M+H]+ 下 10 ppm 内的第 1 名（缓存），单个目标；存盘中同分子式 + 加合物的 cpu 运行被复用（其余目标分子式不同，合并没有丢弃它的候选）。
- 同样需要同一个 `MZDESK_CODE_ROOT`；日志 `standalone.log`；重新启动会续跑，出错的题（`*.error.json`）重跑。

- `collect.py` 与 `replay.py` 必须用同一个 `MZDESK_CODE_ROOT`：存盘的是 mzdesk 的对象，`replay.py` 比对代码指纹，不同则拒绝（`--allow-code-mismatch` 可强行）。
  不设 `MZDESK_CODE_ROOT` 则两者都拒绝启动，永远不会 import 工作区里的 `mzdesk`。
- 并行、续跑、锁与日志见 `common.py`：每个工作进程用基于 `flock` 的文件锁替换进程内的重型锁（同时最多 `--heavy-slots` 个 MS-FINDER/FRIGID），
  `max_tasks_per_child=None`（Python 3.12 的 `ProcessPoolExecutor` 在更换工作进程时会卡住），每题超时（collect 3600 秒，replay 每配置 1800 秒）记为错误，
  工作进程被杀则重启进程池、同一题两次被杀记为错误。日志 `collect.log` / `replay.log`：每题一行，出错另有 `ERROR` 行，结束写 `WALL_SECONDS` 与 `ALL DONE n_ok n_err`；各进程输出在 `worker_logs/`。
- 续跑：`collect.py` 跳过已有 `store/challenge_NNN.pkl` + `.meta.json` 的题（出错的 `.error.json` 会重跑）；`replay.py` 逐配置跳过已有 `replay/<配置>/challenge_NNN.json` 的题。
- 子集：`--ids 9,116,227`，`--out DIR`（输出到别处），`replay.py --configs full,leakfree_full`，`replay.py --store DIR`（collect 的输出目录）。`score.py --partial --dir DIR` 和 `write_report.py --dir DIR --target FILE` 可检查这类输出。
- 每题记录：代码指纹（`mzdesk/**/*.py` 的哈希）、每步耗时、FRIGID 实际设备（mps / cpu，每次调用）、各层状态（`store/challenge_NNN.meta.json`）。
- `collect.py` 的 `manifest.json` 记录设置与每题输入校验值，并拒绝在代码指纹变化后续跑（换冻结代码请换 `--out`）。

## 输出

```
store/          collect：challenge_NNN.pkl（重放所需的全部原始输出）、.meta.json、.error.json
replay/<配置>/  challenge_NNN.json（紧凑结果：排序候选、各步状态、答案、第 7 步缓存命中数）
step7_cache/    第 7 步 MS-FINDER 结果（按参数）
scores.json     逐题、逐配置、逐工具的名次，子集定义，运行时间，核对问题
```

## 去除泄漏的定义（leak_lists/）

`score.py` 由 `leak_lists/overlap.json` 与 `lib_ids.json` 生成每个工具的"已见"集合：
CSI:FingerID = SIRIUS 同极性指纹预测器训练结构（`sirius_pos` / `sirius_neg`）；FRIGID = `msg_train` ∪ `msg4m`；MIST = `mist_real_train_plus_aug`；
JESTR = `jestr_nplib1_train`；谱库 = `lib_any`（本地 MassBank/GNPS 中有该结构的任何一条谱，同 `q1.json` 的 `any`）。
FRIGID 的扩散语言模型预训练语料与 MSNovelist 的训练结构无法排除。所有工具都未见的题共 202 题（`meta.subsets.unseen_common`）。
