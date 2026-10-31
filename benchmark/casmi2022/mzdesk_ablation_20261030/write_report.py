#!/usr/bin/env python3
"""Make the readable report (results only) from scores.json -> benchmark/result/CASMI2022_mzdesk_full500.md.

  python3 benchmark/casmi2022/mzdesk_ablation_20261030/write_report.py [--dir DIR] [--target FILE] [--replace]
"""
from __future__ import annotations

import argparse
import json
import statistics
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TARGET = ROOT / "benchmark/result/CASMI2022_mzdesk_full500.md"
COMPARATORS = ("CSI:FingerID", "MSNovelist", "JESTR", "MIST-contrastive")
POSITIVE_ONLY = ("frigid_only", "default_frigid")
LABELS = {
    "full": "mzDesk 完整流程 (full)", "leakfree_full": "mzDesk 无泄漏（不搜库，直接从头注释）(leakfree_full)",
    "all_layers": "三层合并 (all_layers)", "library_only": "只用谱库 (library_only)", "msfinder_only": "只用 MS-FINDER (msfinder_only)",
    "frigid_only": "只用 FRIGID (frigid_only，仅正离子)", "no_library": "去掉谱库 (no_library)", "no_msfinder": "去掉 MS-FINDER (no_msfinder)",
    "no_frigid": "去掉 FRIGID (no_frigid)", "buddy_only": "只用 BUDDY 分子式 (buddy_only)",
}
ABLATIONS = ("all_layers", "library_only", "no_msfinder", "no_frigid", "msfinder_only", "frigid_only", "no_library")
LF = "mzDesk 无泄漏（不搜库）"


class Data:
    def __init__(self, scores):
        self.meta, self.cases = scores["meta"], scores["cases"]
        self.by_id = {c["id"]: c for c in self.cases}
        self.subsets = {k: set(v) for k, v in self.meta["subsets"].items()}

    def rank(self, case, key):
        """The rank dict of one case for key = comparator name or configuration name (None when not applicable)."""
        if key in COMPARATORS:
            return case["comparators"][key]
        r = case["configs"].get(key)
        if r is not None and key in POSITIVE_ONLY and case["polarity"] != "positive":
            return None
        return r

    def rows(self, key, subset="all", polarity=None):
        out = []
        for case in self.cases:
            if case["id"] not in self.subsets[subset]:
                continue
            if polarity and case["polarity"] != polarity:
                continue
            r = self.rank(case, key)
            if r is not None:
                out.append(r)
        return out


def stats(rows):
    n = len(rows)
    top = lambda k: sum(r["worst"] is not None and r["worst"] <= k for r in rows)
    return {"n": n, "wrong": sum(bool(r.get("wrong_adduct")) for r in rows), "no_output": sum(not r["n"] and not r.get("wrong_adduct") for r in rows), "present": sum(r["worst"] is not None for r in rows),
            "t1": top(1), "t3": top(3), "t10": top(10)}


def cell(k, n):
    return f"{k}/{n} ({100 * k / n:.1f}%)" if n else "-"


def line(label, rows, extra=()):
    s = stats(rows)
    return f"| {label} | {s['n']} | {s['no_output']} | {s['present']} | {cell(s['t1'], s['n'])} | {cell(s['t3'], s['n'])} | {cell(s['t10'], s['n'])} |" + "".join(f" {x} |" for x in extra)


HEAD = ["| 方法 | 题数 | 无输出 | 真值在列表中 | Top-1 | Top-3 | Top-10 |", "|---|---:|---:|---:|---:|---:|---:|"]


def hms(seconds):
    if seconds is None:
        return "无记录"
    return f"{int(seconds // 3600)} 小时 {int(seconds % 3600 // 60)} 分" if seconds >= 3600 else f"{seconds / 60:.1f} 分钟"


def dist(values):
    values = [v for v in values if v is not None]
    return f"平均 {statistics.mean(values):.1f}，中位数 {statistics.median(values):.1f}，最长 {max(values):.1f}" if values else "无记录"


# ------------------------------------------------------------------ sections
def main_tables(d):
    text = ["## (a) 主结果：mzDesk 完整流程与对照工具", "",
            "mzDesk (full) 为完整流程（盲测，库内含同结构谱时未去除）；\"无泄漏\"一行完全不搜谱库（第 3 步无输出、闸门不通过），直接从头注释："
            "第 4–5 步按流程自己的判断（MS-FINDER；MS-FINDER 不够好且为正离子时才跑 FRIGID），再合并、重排、同分异构体检查，不读真值。", ""]
    for title, pol in (("全部", None), ("正离子", "positive"), ("负离子", "negative")):
        n = len([c for c in d.cases if pol is None or c["polarity"] == pol])
        text += [f"### {title}（{n} 题）", ""] + HEAD
        for label, key in (("mzDesk (full)", "full"), (LF, "leakfree_full")) + tuple((c, c) for c in COMPARATORS):
            rows = d.rows(key, "all", pol)
            if rows:
                text.append(line(label, rows))
        text.append("")
    return text


def ablation_tables(d):
    text = ["## (b) 消融（仅看排名，不含第 7 步）", "",
            "每个配置只保留所列层的候选，各层结果一律并入（不论谱库闸门或 MS-FINDER \"够好\"的判定），再走流程自己的合并与重排。"
            "前两行（含第 7 步）是流程的实际运行，不是强制并入。", "",
            "| 配置 | 题数 | 无输出 | 真值在列表中 | Top-1 | Top-3 | Top-10 | 正离子 Top-1 | 负离子 Top-1 |", "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]

    def abl(label, key):
        rows = d.rows(key)
        s = stats(rows)
        pos, neg = stats(d.rows(key, "all", "positive")), stats(d.rows(key, "all", "negative"))
        return (f"| {label} | {s['n']} | {s['no_output']} | {s['present']} | {cell(s['t1'], s['n'])} | {cell(s['t3'], s['n'])} | {cell(s['t10'], s['n'])} | "
                f"{cell(pos['t1'], pos['n'])} | {cell(neg['t1'], neg['n']) if neg['n'] else '不适用'} |")
    text.append(abl(LABELS["full"] + "，含第 7 步", "full"))
    text.append(abl(LABELS["leakfree_full"] + "，含第 7 步", "leakfree_full"))
    for name in ABLATIONS:
        text.append(abl(LABELS[name], name))
    text += ["", "### 分子式层（BUDDY）", "",
             "第 2 步的分子式列表（各加合物按 BUDDY 排名交错），判定：真值的中性分子式在列表中的位置（相同分子式只算一次）。不涉及结构。", "",
             "| 子集 | 题数 | 无输出 | 真值分子式在列表中 | Top-1 | Top-3 | Top-10 |", "|---|---:|---:|---:|---:|---:|---:|"]
    for title, pol in (("全部", None), ("正离子", "positive"), ("负离子", "negative")):
        text.append(line(title, d.rows("buddy_only", "all", pol)))
    text.append("")
    return text


def leakfree_tables(d):
    sub = d.subsets
    common = "unseen_common"
    pos = [c for c in d.cases if c["id"] in sub[common] and c["polarity"] == "positive"]
    neg = [c for c in d.cases if c["id"] in sub[common] and c["polarity"] == "negative"]
    seen = d.meta["seen"]
    text = ["## (c) 泄漏控制：所有被比较工具都没见过真值结构的题", "",
            f"共 {len(sub[common])} 题（正离子 {len(pos)}、负离子 {len(neg)}）：真值 IK14 不在下列任何一项中——SIRIUS 指纹预测器训练集（同极性）、"
            f"FRIGID 的 MSG 训练集 ∪ MSG-4M、MIST 的真实+增强训练集、JESTR 的 NPLIB1 训练集、本地 MassBank/GNPS 中任何一条该结构的谱。"
            "名单与来源见 `leak_lists/`。MSNovelist 与 FRIGID 的扩散语言模型预训练语料无法排除，其结果在此子集上仅供参考（标 †）。", "",
            "各工具已见题数（真值结构出现在其训练/参考集中）：" + "；".join(f"{k} {len(v)}" for k, v in seen.items()) + "。", ""]
    for title, pol in (("全部", None), ("正离子", "positive"), ("负离子", "negative")):
        n = len([c for c in d.cases if c["id"] in sub[common] and (pol is None or c["polarity"] == pol)])
        text += [f"### {title}（{n} 题）", ""] + HEAD
        for label, key in (("mzDesk (full)", "full"), (LF, "leakfree_full"), ("CSI:FingerID", "CSI:FingerID"), ("MSNovelist †", "MSNovelist"),
                           ("JESTR", "JESTR"), ("MIST-contrastive", "MIST-contrastive")):
            rows = d.rows(key, common, pol)
            if rows:
                text.append(line(label, rows))
        text.append("")
    text += ["### 同一子集上的 mzDesk 各配置（排名，不含第 7 步；该子集无谱库同结构命中）", ""] + HEAD
    for name in ("all_layers", "library_only", "msfinder_only", "frigid_only", "no_library", "no_msfinder", "no_frigid"):
        rows = d.rows(name, common)
        if rows:
            text.append(line(LABELS[name], rows))
    text += ["", "### 每个工具在自己的\"未见\"子集上", "",
             "各工具只排除自己见过的题；子集大小不同，数值不能横向严格比较。\"全部 500 题\"列为同一工具在全部题上的 Top-1，作参照。MIST 只含正离子题。", "",
             "| 工具 | 未见子集 | 子集题数（评分） | 无输出 | Top-1 | Top-3 | Top-10 | 全部 500 题 Top-1 |", "|---|---|---:|---:|---:|---:|---:|---:|"]
    own = (("mzDesk (full)", "full", "unseen_mzDesk", "不在谱库、不在 FRIGID 训练集；MS-FINDER 无可排除的训练集；FRIGID DLM 无法排除"),
           (LF, "leakfree_full", "unseen_FRIGID", "不搜库，只需不在 FRIGID 训练集；MS-FINDER 无可排除的训练集；FRIGID DLM 无法排除"),
           ("CSI:FingerID", "CSI:FingerID", "unseen_CSI:FingerID", "不在 SIRIUS 同极性训练集"),
           ("MSNovelist", "MSNovelist", "all", "无法排除，用全部 500 题"),
           ("JESTR", "JESTR", "unseen_JESTR", "不在 NPLIB1 训练集"),
           ("MIST-contrastive", "MIST-contrastive", "unseen_MIST", "不在 MIST 训练+增强"),
           ("FRIGID 层 (frigid_only)", "frigid_only", "unseen_FRIGID", "不在 MSG 训练 ∪ MSG-4M；DLM 预训练无法排除"))
    for label, key, subset, why in own:
        rows = d.rows(key, subset)
        s = stats(rows)
        all_rows = stats(d.rows(key, "all"))
        text.append(f"| {label} | {subset.replace('unseen_', '未见 ')} ({why}) | {s['n']} | {s['no_output']} | {cell(s['t1'], s['n'])} | {cell(s['t3'], s['n'])} | "
                    f"{cell(s['t10'], s['n'])} | {cell(all_rows['t1'], all_rows['n'])} |")
    text.append("")
    return text


def default_section(d):
    """(d) single tools with default settings next to the matching mzDesk-input ablation and mzDesk full."""
    sa = d.meta["standalone"]
    n_wrong = len(d.subsets["default_adduct_wrong"])
    n285 = len(d.subsets["default_adduct"])
    wrong_pos = sum(d.by_id[i]["polarity"] == "positive" for i in d.subsets["default_adduct_wrong"])
    pos_total = sum(c["polarity"] == "positive" for c in d.cases)
    text = ["## (d) 单软件默认设置：不经过 mzDesk 的加合物、分子式与合并步骤", "",
            "每个工具单独、用默认设置运行，与\"mzDesk 输入\"的对应消融（同一工具，但加合物与分子式来自 mzDesk 的第 1–2 步，再经 mzDesk 的合并与重排）以及 mzDesk 完整流程并列，"
            "用来衡量 mzDesk 的总贡献（加合物 + 分子式 + 合并）。",
            "",
            f"- 默认加合物：正离子 [M+H]+，负离子 [M-H]-。真实加合物是别的（[M+NH4]+、[M+Na]+、[M+FA-H]-、[M+Cl]-）的题**直接判错**，不运行、仍在分母中："
            f"{n_wrong} 题（正离子 {wrong_pos}、负离子 {n_wrong - wrong_pos}）；只有其余 {n285} 题（正离子 {n285 - sum(1 for i in d.subsets['default_adduct'] if d.by_id[i]['polarity'] == 'negative')}、"
            f"负离子 {sum(1 for i in d.subsets['default_adduct'] if d.by_id[i]['polarity'] == 'negative')}）真正运行。真值只用来选题和判分，不作为任何工具的输入。",
            "- 谱库：不用加合物，500 题全部运行；前体 10 ppm，MassBank + GNPS，加权谱熵 ≥ 0.70 且匹配峰 ≥ 4，按 IK14 去重，按相似度排序；没有 mzDesk 的闸门检查（加合物/分子式/同位素），不截断。",
            "- BUDDY：不重新运行，取缓存的 BUDDY 输出中默认加合物、10 ppm 以内的分子式，按 BUDDY 自己的顺序（不按质量误差重排），判定真值中性分子式的位置。"
            "MS-FINDER：默认加合物，分子式由 MS-FINDER 自己生成（mzDesk MS-FINDER 层的\"无分子式列表\"模式，10 ppm / 0.01 Da、同一结构库、TreeDepth 与超时重试、top_n 10 都与 mzDesk 的 MS-FINDER 层相同；卤素按 MS-FINDER 默认关闭，MS1 同位素峰交给 MS-FINDER）。",
            f"- FRIGID：只有正离子 [M+H]+（分母为全部正离子题 {pos_total} 题，与 frigid_only 相同）；分子式取 BUDDY 在 [M+H]+ 下 10 ppm 以内的第 1 名，设置同 mzDesk 的 FRIGID 层（top_n 5，种子 42，cpu）。"
            "单软件默认行的排序是各工具自己的分数（谱库相似度、MS-FINDER 总分、FRIGID 指纹相似度），不经 mzDesk 的合并与重排。",
            ""]
    COLS = ["| 方法 | 题数 | 无输出 | 加合物错误（判错） | 真值在列表中 | Top-1 | Top-3 | Top-10 | 正离子 Top-1 | 负离子 Top-1 |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]

    def row(label, key, subset="all", positive_only=False):
        rows = d.rows(key, subset, "positive" if positive_only else None)
        if not rows:
            return None
        s = stats(rows)
        pos, neg = stats(d.rows(key, subset, "positive")), stats(d.rows(key, subset, "negative") if not positive_only else [])
        return (f"| {label} | {s['n']} | {s['no_output']} | {s['wrong']} | {s['present']} | {cell(s['t1'], s['n'])} | {cell(s['t3'], s['n'])} | {cell(s['t10'], s['n'])} | "
                f"{cell(pos['t1'], pos['n'])} | {cell(neg['t1'], neg['n']) if neg['n'] else '不适用'} |")

    def table(subset):
        spec = [("谱库，**默认**", "default_library", False)]
        spec += [("谱库，mzDesk 输入 (library_only)", "library_only", False)]
        spec += [("MS-FINDER，**默认**", "default_msfinder", False), ("MS-FINDER，mzDesk 输入 (msfinder_only)", "msfinder_only", False),
                 ("FRIGID（仅正离子），**默认**", "default_frigid", True), ("FRIGID（仅正离子），mzDesk 输入 (frigid_only)", "frigid_only", True),
                 ("BUDDY 分子式，**默认**", "default_buddy", False), ("BUDDY 分子式，mzDesk 输入 (buddy_only)", "buddy_only", False),
                 ("mzDesk 完整流程 (full)", "full", False), ("mzDesk 完整流程 (full)，仅正离子", "full", True)]
        spec += [("mzDesk 无泄漏（不搜库，直接从头注释）(leakfree_full)", "leakfree_full", False)]
        out = list(COLS)
        for label, key, pos_only in spec:
            line_ = row(label, key, subset, pos_only)
            if line_:
                out.append(line_)
        return out

    text += ["### 全部 500 题（FRIGID 行的分母是 %d 题正离子）" % pos_total, "",
             "BUDDY 行判定的是中性分子式的位置（不涉及结构），\"mzDesk 完整流程\"没有对应的分子式行。", ""] + table("all") + [""]
    sub = d.subsets["unseen_common"]
    text += [f"### {len(sub)} 题无泄漏子集（同 (c)，该子集无谱库同结构命中）", ""] + table("unseen_common") + [""]

    text += [f"### 只看默认加合物就是真实加合物的 {n285} 题（[M+H]+ / [M-H]-）", "",
             "这 285 题上默认加合物等于真实加合物，没有\"加合物错误\"的判错，对比的是分子式来源与合并步骤的差别；同一题集上给出对应消融与完整流程。FRIGID 行的分母是其中的正离子题。", ""]
    h = ["| 方法 | 题数 | 无输出 | 真值在列表中 | Top-1 | Top-3 | Top-10 |", "|---|---:|---:|---:|---:|---:|---:|"]
    spec = [("谱库，默认", "default_library", False), ("谱库，mzDesk 输入 (library_only)", "library_only", False),
            ("MS-FINDER，**默认**", "default_msfinder", False), ("MS-FINDER，mzDesk 输入 (msfinder_only)", "msfinder_only", False),
            ("FRIGID（仅正离子），**默认**", "default_frigid", True), ("FRIGID（仅正离子），mzDesk 输入 (frigid_only)", "frigid_only", True),
            ("BUDDY 分子式，**默认**", "default_buddy", False), ("BUDDY 分子式，mzDesk 输入 (buddy_only)", "buddy_only", False),
            ("mzDesk 完整流程 (full)", "full", False), ("mzDesk 完整流程 (full)，仅正离子", "full", True)]
    text += h
    for label, key, pos_only in spec:
        rows = d.rows(key, "default_adduct", "positive" if pos_only else None)
        if rows:
            text.append(line(label, rows))
    text.append("")

    # how the outputs were obtained
    ver = sa.get("frigid_verify") or []
    src = sa.get("frigid_sources") or {}
    text += ["### 运行与核对", "",
             f"- 运行：谱库 500 题；BUDDY（缓存）、MS-FINDER 与 FRIGID 只在默认加合物等于真实加合物的 {n285} 题（FRIGID 其中的正离子题）上；墙钟 {hms(sa.get('wall_seconds'))}（各次启动之和，含 FRIGID 复用核对的重跑和冒烟测试）。",
             "- FRIGID 输出来源：" + "；".join(f"{ {'run': '本次实际运行', 'none': 'BUDDY 无分子式，未运行'}.get(k, '沿用 collect 存盘（同分子式、加合物、设置）')} {v} 题" for k, v in sorted(src.items())) + "。"
             + (f" 复用核对：重跑 {len(ver)} 题，与 collect 存盘输出完全一致 {sum(bool(v.get('identical')) for v in ver)} 题。" if ver else ""),
             "- BUDDY 缓存核对：" + "；".join(f"{k} {v}" for k, v in sorted(sa.get("buddy_cache_checks", {}).items())) + "（与 collect 存盘的第 2 步分子式列表逐项比较）。", ""]
    return text


def runtime_section(d):
    cases = d.cases
    col = lambda f: [f(c) for c in cases]
    text = ["## (e) 运行时间", "",
            f"- 重型一遍（collect.py，第 1–5 步无条件运行）：总墙钟 {hms(d.meta['collect_wall_seconds'])}（3 个并行进程，同时最多 2 个 MS-FINDER/FRIGID 任务）。",
            f"- 每题（秒）：分子式+加合物 {dist(col(lambda c: c['collect']['seconds'].get('adduct+formula')))}",
            f"  - 谱库 {dist(col(lambda c: c['collect']['seconds'].get('library')))}",
            f"  - MS-FINDER {dist(col(lambda c: c['collect']['seconds'].get('msfinder')))}",
            f"  - FRIGID（仅正离子有实际运行） {dist(col(lambda c: c['collect']['seconds'].get('frigid') if c['polarity'] == 'positive' else None))}",
            f"  - 合计（含等待重型引擎） {dist(col(lambda c: c['collect']['seconds'].get('total')))}"]
    wf = [c["workflow_seconds"]["total"] for c in cases if c.get("workflow_seconds")]
    if wf:
        text.append(f"- 按 full 的实际跳过判断折算，用户单题运行流程的计算时间（秒；各步用 collect 时的单次计时，第 7 步为引擎时间，不含等待）：{dist(wf)}。")
        skip4 = sum(c["config_info"]["full"]["step_status"]["4"] == "skipped" for c in cases if "full" in c["config_info"])
        skip5 = sum(c["config_info"]["full"]["step_status"]["5"] == "skipped" for c in cases if "full" in c["config_info"])
        text.append(f"- full 中第 4 步（MS-FINDER）被跳过（谱库闸门通过）：{skip4} 题；第 5 步（FRIGID）状态为跳过：{skip5} 题（含闸门通过、MS-FINDER 够好、以及负离子引擎不支持）。")
    wl = [c["workflow_seconds_leakfree"]["total"] for c in cases if c.get("workflow_seconds_leakfree")]
    if wl:
        skip5l = sum(c["config_info"]["leakfree_full"]["step_status"]["5"] == "skipped" for c in cases if "leakfree_full" in c["config_info"])
        text.append(f"- leakfree_full（不搜库，第 4 步一定运行）同样折算：{dist(wl)}；第 5 步（FRIGID）状态为跳过 {skip5l} 题（MS-FINDER 够好或负离子）。")
    dev = Counter("+".join(c["collect"]["frigid_devices"]) for c in cases if c["polarity"] == "positive")
    text.append("- FRIGID 实际设备（正离子题，由当时可用内存决定）：" + "，".join(f"{k} {v}" for k, v in sorted(dev.items())) + "。")
    text.append(f"- 重放（replay.py）：总墙钟 {hms(d.meta['replay_wall_seconds'])}（各次启动之和，含第 7 步的 MS-FINDER，按参数缓存，full 与 leakfree_full 共用）。")
    per = []
    for name in d.meta["configs"]:
        sec = [c["config_info"][name]["replay_seconds"] for c in cases if name in c["config_info"] and "replay_seconds" in c["config_info"][name]]
        if sec:
            per.append(f"{name} {statistics.mean(sec):.2f}")
    text.append("  - 每题平均秒数（full、leakfree_full 含第 7 步的 MS-FINDER，leakfree_full 中取 full 结果的题记 0；其余为合并+重排）：" + "；".join(per))
    runs = [c["config_info"][n]["step7_cache"] for c in cases for n in ("full", "leakfree_full") if n in c["config_info"] and c["config_info"][n].get("step7_cache")]
    text.append(f"- 第 7 步（对 top-1 及同分子式候选再跑一次 MS-FINDER）：full 与 leakfree_full 的最终结果中真实运行 {sum(r['misses'] for r in runs)} 次，共 "
                f"{sum(r['miss_seconds'] for r in runs) / 60:.1f} 分钟，缓存命中 {sum(r['hits'] for r in runs)} 次（两个配置顶位及同分子式候选相同时共用）。")
    for n in ("full", "leakfree_full"):
        st = Counter(c["config_info"][n]["step_status"]["7"] for c in cases if n in c["config_info"] and "step_status" in c["config_info"][n])
        text.append(f"  - {n} 第 7 步状态：" + "，".join(f"{k} {v}" for k, v in sorted(st.items())) + "。")
    text.append("")
    return text


def settings_section(d):
    m = d.meta

    def step7_failed(name):
        return sum(c["config_info"].get(name, {}).get("step_status", {}).get("7") == "failed" for c in d.cases)
    problems = m["problems"]
    ids_with_problems = "、".join(str(i) for i in sorted({int(x.split()[0]) for x in problems})) or "无"
    lf = [c["config_info"]["leakfree_full"] for c in d.cases if "leakfree_full" in c["config_info"]]
    copied = sum(bool(x.get("copied_from_full")) for x in lf)
    replayed, proved = len(lf) - copied, m.get("leakfree_equivalence", {}).get("compared", 0)
    errors = [c["id"] for c in d.cases for k, r in c["configs"].items() if r.get("error")]
    layer_fail = Counter(k for c in d.cases for k, v in c["collect"]["layers"].items() if v["status"] == "failed")
    return ["## (f) 设置与限制", "",
            "- mzDesk 设置：实测前体及 MS1/MS2，加合物自动，10 ppm 前体 / 0.01 Da 碎片，允许卤素，MassBank + GNPS 全库，"
            "Metabolomics，甲酸，碱金属未知，pH 2.7，不查数据库同分异构体；碰撞能谱按 10 ppm 合并。BUDDY 内部使用 10 ppm MS1 / 20 ppm MS2。",
            f"- 两遍结构：第 1–5 步在冻结代码（`{m.get('code_root')}`，指纹 {', '.join(m['code_fingerprints'])}）上无条件运行一次；第 6–8 步由流程自己的 `run_query` "
            "在存下的输出上重放（引擎入口返回缓存），所以闸门/跳过判断、合并、重排、答案逻辑与流程完全相同。`verify.py` 核对过重放与正常运行一致（见 README）。",
            "- 消融不含第 7 步：同分异构体检查不改变任何排名或级别（`verify.py` 逐题核对），消融配置里它被替换为\"跳过\"，答案类型不可用。",
            "- mzDesk 无泄漏 (leakfree_full)：不搜谱库，第 3 步无输出、闸门不通过，直接用第 4–5 步的存盘输出（第 4 步 MS-FINDER 一定运行；第 5 步 FRIGID 仅在 MS-FINDER 不够好且为正离子时运行），"
            "再合并、重排、第 7 步、答案，全部由流程自己的 `run_query` 决定，不读真值。它的第 4–5 步输出与 full 用的是同一份存盘输出，因此与 full 的差别只来自谱库层。",
            f"- leakfree_full：存盘谱库步骤无错误且 0 个候选的 {copied} 题直接取 full 的结果（此时去掉谱库不改变任何输入；两种做法都真实重放过的 {proved} 题逐项一致），其余 {replayed} 题真实重放。",
            "- 谱库命中在流程里被截到 20 个化合物，重放同样截断到 20。",
            "- 第 7 步与 FRIGID 结果依赖设备与引擎，个别题可能不可复现；FRIGID 固定种子 42，但设备（GPU/CPU）不同会导致分数微小差异。",
            "- 泄漏控制只能排除已知训练/参考集；MSNovelist 与 FRIGID 扩散语言模型的预训练语料无法排除；MS-FINDER 无可查的训练集。\"未见\"子集依据 IK14，不考虑立体化学。",
            "- 判对：IK14 一致，不判断立体化学；平局用最差名次；无输出的题在分母中；MIST 只支持正离子，其分母只含正离子题。",
            "- CSI:FingerID 和 MSNovelist 依赖在线服务，属场景外参考；各工具的候选库、谱图处理和算法不统一，不能作严格公平比较。仅用于探索。",
            f"- 运行情况：collect 阶段各层失败数 {dict(layer_fail) or '无'}；配置出错 {len(errors)} 次；核对不一致 {len(problems)} 项（见 scores.json 的 meta.problems；涉及题号 {ids_with_problems}，原因是引擎给出的某个 SMILES 无法被 RDKit 解析，这些题改用显示名次）；"
            f"第 7 步失败 full {step7_failed('full')} 题、leakfree_full {step7_failed('leakfree_full')} 题。",
            "- 单软件默认设置（(d)）：MS-FINDER 使用 `~/tools/msfinder` 的一份副本运行——该安装里的 `MsfinderFormulaDB-VS13_bin.efd` 在 2026-10-30 17:52 变成 0 字节，"
            "此后任何 MS-FINDER 运行都会\"finished without writing a result file\"；副本中该文件换回随发行包的原件（与 `.backup` 相同），其余原样。BUDDY 默认行来自缓存的 BUDDY 输出，不重新运行。",
            "- 复现命令与顺序见 [README](../casmi2022/mzdesk_ablation_20261030/README.md)；[逐题机器可读结果](../casmi2022/mzdesk_ablation_20261030/scores.json)。", ""]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=str(HERE), help="folder holding scores.json (default: this folder)")
    ap.add_argument("--target", default=str(TARGET), help="report file (smoke tests must not use the default)")
    ap.add_argument("--replace", action="store_true", help="replace the previous report")
    args = ap.parse_args()
    d = Data(json.loads((Path(args.dir) / "scores.json").read_text()))
    if not d.meta["partial"]:
        assert len(d.cases) == 500
    n_pos = sum(c["polarity"] == "positive" for c in d.cases)
    text = ["# CASMI 2022：mzDesk 全部 500 题——完整流程、消融与泄漏控制", "",
            f"共 {len(d.cases)} 题，其中正离子 {n_pos} 题、负离子 {len(d.cases) - n_pos} 题。仅用于探索，不足以判断工具优劣。", ""]
    text += main_tables(d) + ablation_tables(d) + leakfree_tables(d) + default_section(d) + runtime_section(d) + settings_section(d)
    target = Path(args.target)
    if target.exists() and not args.replace:
        raise FileExistsError(target)
    target.write_text("\n".join(text), encoding="utf-8")
    print(target)


if __name__ == "__main__":
    main()
