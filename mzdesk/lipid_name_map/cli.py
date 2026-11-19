"""命令行入口（W3）：``python3 -m lipid_name_map``。

对外接口
--------
``main(argv=None) -> int``       解析参数、处理、返回退出码（``__main__`` 用 ``raise SystemExit(main())`` 调用）
``build_parser() -> ArgumentParser``
``CSV_COLUMNS``                  ``--out *.csv`` 的列：``input status reason level class_full L1 L2 L3 LM_species LM_molecular flags``
``to_row(result) -> dict``       ``normalize`` 返回的字典 -> CSV 行（``flags`` 用 ``;`` 连接）

用法（在项目根目录，或包已安装）::

    python3 -m lipid_name_map "PC 16:0/18:1"                 # 单个名字：竖排显示
    python3 -m lipid_name_map "PC 16:0/18:1" --json          # 单个名字：完整 JSON（多个名字或批量时每行一个 JSON）
    python3 -m lipid_name_map "PC 18:1" --source pnnl        # 来源提示（目前只有 pnnl 影响解析）
    python3 -m lipid_name_map --file in.txt --out out.csv    # 每行一个名字 -> CSV（扩展名 .csv/.tsv/.jsonl 决定格式）
    cat in.txt | python3 -m lipid_name_map --out out.jsonl   # 无名字、无 --file 时读 stdin（``--file -`` 同义）

输入文件按 UTF-8 读（可带 BOM；坏字节替换成 U+FFFD，由 ``normalize`` 拒绝）。**输入的每一行对应输出的一行**，
空行输出 ``rejected/empty``，不会被跳过，方便和原文件逐行对齐。

输出格式（``--format``，缺省按 ``--out`` 的扩展名；没有 ``--out`` 时：命令行只给一个名字 = ``block``，其余（多个名字、``--file``、stdin）= ``csv``）：
``block`` 竖排可读 · ``csv`` · ``tsv`` · ``jsonl``（每行完整 ``normalize`` 字典，含 ``meta``） · ``json``（单个名字缩进 JSON）。
``--json`` 是 ``--format json``（多个名字时退化成 ``jsonl``）的简写。

退出码：``0`` 成功（名字作为参数给出时，全部 ``ok``；批量/stdin/``--file`` 时只要处理完成就是 0，拒绝只是数据）；
``1`` 名字作为参数给出时，至少一个不是 ``ok``；``2`` 用法错误或读写文件失败。
批量模式结束时向 stderr 打印一行汇总（``-q`` 关闭）。
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections import Counter
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import IO, Any

from . import __version__
from .config import DEFAULT, SOURCES, Config

__all__ = ["main", "build_parser", "CSV_COLUMNS", "to_row", "read_names"]

CSV_COLUMNS: tuple[str, ...] = ("input", "status", "reason", "level", "class_full", "L1", "L2", "L3", "LM_species",
                                "LM_molecular", "flags")
_FORMATS = ("block", "csv", "tsv", "jsonl", "json")
_EXT_FORMAT = {".csv": "csv", ".tsv": "tsv", ".jsonl": "jsonl", ".ndjson": "jsonl"}


# ----------------------------------------------------------------------------------------------------------------
# 输入 / 输出
# ----------------------------------------------------------------------------------------------------------------
def _decode(raw: bytes) -> str:
    return raw.decode("utf-8-sig", errors="replace")


def _split_lines(text: str) -> list[str]:
    """按行切（``\\r\\n`` 与 ``\\n``）；末尾的换行不产生多余的空行。"""
    lines = text.replace("\r\n", "\n").split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    return lines


def read_names(path: str | None, stdin: IO[Any] | None = None) -> list[str]:
    """读名字列表：``path`` 为 ``None`` 或 ``"-"`` 时读 stdin；否则读文件。每行一个名字，行数不变。"""
    if path is None or path == "-":
        stream = stdin if stdin is not None else sys.stdin
        buf = getattr(stream, "buffer", None)
        text = _decode(buf.read()) if buf is not None else str(stream.read())
    else:
        text = _decode(Path(path).read_bytes())
    return _split_lines(text)


def to_row(result: dict) -> dict[str, Any]:
    """``normalize`` 的返回字典 -> CSV 行（``CSV_COLUMNS``）；``flags`` 用 ``;`` 连接。"""
    return {"input": result["input"] if isinstance(result["input"], str) else repr(result["input"]),
            "status": result["status"], "reason": result["reason"], "level": result["level"],
            "class_full": result["class_full"], "L1": result["L1"], "L2": result["L2"], "L3": result["L3"],
            "LM_species": result["LM_species"], "LM_molecular": result["LM_molecular"],
            "flags": ";".join(result["flags"])}


def _format_block(r: dict) -> str:
    """单个名字的竖排显示。"""
    lines = [f"input   {r['input']}"]
    if r["status"] == "ok":
        lines.append(f"status  ok (level {r['level']})")
        lines.append(f"class   {r['class_full']}")
        lines.append(f"L1      {r['L1']}")
        lines.append(f"L2      {r['L2'] or '-'}")
        lines.append(f"L3      {r['L3'] or '-'}")
        lines.append(f"LM      {r['LM_species'] or '-'}")
        lines.append(f"LM_mol  {r['LM_molecular'] or '-'}")
    else:
        lines.append(f"status  {r['status']} ({r['reason']})")
        detail = r["meta"].get("detail")
        if detail:
            lines.append(f"detail  {detail}")
        hint = r["meta"].get("class_hint")
        if hint:
            lines.append(f"head    {hint}")
    if r["flags"]:
        lines.append(f"flags   {', '.join(r['flags'])}")
    return "\n".join(lines)


def _write_results(results: Sequence[dict], fmt: str, out: IO[str]) -> None:
    if fmt == "block":
        out.write("\n\n".join(_format_block(r) for r in results) + "\n")
    elif fmt in ("csv", "tsv"):
        w = csv.DictWriter(out, fieldnames=CSV_COLUMNS, delimiter="," if fmt == "csv" else "\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        for r in results:
            w.writerow(to_row(r))
    elif fmt == "json" and len(results) == 1:
        out.write(json.dumps(results[0], ensure_ascii=False, indent=2) + "\n")
    else:                                                    # jsonl（以及多个名字的 json）
        for r in results:
            out.write(json.dumps(r, ensure_ascii=False) + "\n")


def _summary(results: Iterable[dict], seconds: float) -> str:
    res = list(results)
    c = Counter(r["status"] for r in res)
    rate = len(res) / seconds if seconds > 0 else float("inf")
    return (f"{len(res)} names: ok {c['ok']}, rejected {c['rejected']}, unsupported {c['unsupported']} "
            f"({seconds:.2f}s, {rate:,.0f} names/s)")


# ----------------------------------------------------------------------------------------------------------------
# 参数
# ----------------------------------------------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="python3 -m lipid_name_map",
        description="脂质 shorthand 名 -> 标准化 L1 / L2 / L3（benchmark 去 sn 写法）。输入支持到哪一级就输出到哪一级，其余留空。",
        epilog="退出码：0 成功；1 作为参数给出的名字里有非 ok；2 用法或读写错误。详见 README.md。")
    ap.add_argument("names", nargs="*", metavar="NAME", help="要标准化的名字（带空格时加引号）；不给则读 --file / stdin")
    ap.add_argument("--source", choices=sorted(SOURCES), default=None,
                    help="来源提示；目前只有 pnnl 影响解析（多链类只写一个 C:DB 时当溶血），其余只记入 meta.source_used")
    ap.add_argument("--file", "-f", metavar="PATH", default=None, help="每行一个名字的文本文件（UTF-8）；'-' = stdin")
    ap.add_argument("--out", "-o", metavar="PATH", default=None,
                    help="结果文件（.csv/.tsv/.jsonl 按扩展名定格式；'-' = stdout）；缺省打印到 stdout")
    ap.add_argument("--format", choices=_FORMATS, default=None, help="覆盖输出格式（默认见模块说明）")
    ap.add_argument("--json", action="store_true", help="输出完整 JSON（单个名字缩进，多个名字每行一个）")
    ap.add_argument("--allow-unverified-alias", action="store_true", help="接受 unverified 的类别名（默认拒绝）")
    ap.add_argument("--emit-mecpr", action="store_true", help="名字显式写了 Me/cPr 支链时，在 meta.L3_with_MeCpr 另给带支链的 L3")
    ap.add_argument("--no-oxylipin", action="store_true", help="不展开介质缩写 HETE/HEPE/HDoHE/HODE/EET/EpETE/EpDPE（按俗名拒绝）")
    ap.add_argument("--no-base-oxygen", action="store_true",
                    help="鞘脂名字没给氧数（SM 30:1、SPB 18:1）时不按 d 型补 ;O2，拒绝 missing_oxygen_count（C26 原设计）")
    ap.add_argument("--pnnl-ether-heuristic", action="store_true",
                    help="PNNL 的 PE-P/PC-P 等 +1 落点用 benchmark 的“放饱和链”启发式（不推荐；需配 --source pnnl）")
    ap.add_argument("-q", "--quiet", action="store_true", help="批量模式不向 stderr 打印汇总")
    ap.add_argument("--version", action="version", version=f"lipid_name_map {__version__}")
    return ap


def _config(a: argparse.Namespace) -> Config:
    kw: dict[str, Any] = {}
    if a.allow_unverified_alias:
        kw["allow_unverified_alias"] = True
    if a.emit_mecpr:
        kw["emit_mecpr_field"] = True
    if a.no_oxylipin:
        kw["allow_oxylipin_abbrev"] = False
    if a.no_base_oxygen:
        kw["assume_base_oxygen"] = False
    if a.pnnl_ether_heuristic:
        kw["pnnl_ether_heuristic"] = True
    return Config(**kw) if kw else DEFAULT


def _pick_format(a: argparse.Namespace) -> str:
    if a.format:
        return a.format
    if a.json:
        return "json"
    if a.out and a.out != "-":
        return _EXT_FORMAT.get(Path(a.out).suffix.lower(), "csv")
    return "block" if len(a.names) == 1 else "csv"


# ----------------------------------------------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------------------------------------------
def main(argv: Sequence[str] | None = None) -> int:
    ap = build_parser()
    a = ap.parse_args(argv)
    if a.names and a.file is not None:
        ap.error("give names as arguments or --file, not both")
    if a.json and a.format and a.format not in ("json", "jsonl"):
        ap.error("--json conflicts with --format " + a.format)

    batch = not a.names
    if batch and a.file is None and sys.stdin.isatty():
        ap.print_usage(sys.stderr)
        print("error: no names given (pass NAME..., --file PATH, or pipe names on stdin)", file=sys.stderr)
        return 2

    from . import normalize_many                                # 惰性 import：--help / --version 不加载整套数据
    try:
        names = list(a.names) if a.names else read_names(a.file)
    except OSError as e:
        print(f"error: cannot read input: {e}", file=sys.stderr)
        return 2

    cfg = _config(a)
    t0 = time.perf_counter()
    results = normalize_many(names, source=a.source, config=cfg)
    seconds = time.perf_counter() - t0

    fmt = _pick_format(a)
    try:
        if a.out and a.out != "-":
            with open(a.out, "w", encoding="utf-8", newline="") as fh:
                _write_results(results, fmt, fh)
        else:
            _write_results(results, fmt, sys.stdout)
            sys.stdout.flush()
    except OSError as e:
        print(f"error: cannot write output: {e}", file=sys.stderr)
        return 2

    if batch and not a.quiet:
        print(_summary(results, seconds), file=sys.stderr)
    if not batch and any(r["status"] != "ok" for r in results):
        return 1
    return 0
