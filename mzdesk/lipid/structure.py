"""脂质名字 -> 示意结构（Lipid Short Name 页）。

输入是 mzdesk.lipid_name_map.normalize() 的返回字典（或直接给名字）。名字写到哪一级，结构就按这一级拼：写不到的部分按默认规则猜，
猜出来的原子打上标记，画图时高亮（draw.py）。结构只是"一个与名字相符的示意"，不是实测结构；分子式与名字一致。

SMILES 里的原子序号（atom map）是位标记，读回来后清零：
  FLAG_DB    1  默认 C=C 位置（名字没写位置）
  FLAG_FG    2  默认官能团位置（OH / oxo / OOH / COOH / 环氧）
  FLAG_SPLIT 4  按默认规则拆出来的链（名字只给了总量或没有顺序的链）

支持的类：甘油磷脂（PC PE PS PG PI PA，含 lyso / ether）、DG TG MG（含 ether）、CL、Cer SM SPB Cer-Hex Cer-Hex2、
FA FOH CAR FA-etn FA-amide、胆固醇型固醇酯（ST-oacyl，CE / SE）、FAHFA。其余类没有示意结构（NoStructureForClass）。

出图前有一道分子式核对（illustrate）：画出来的 C、H、N、O、P、S 个数必须等于"类常数 + LIPID+ 解析出的总量 (C, DB, O)"，对不上就不出图（Unsupported），
这样链太短丢了 C=C、未分型的氧没画上或多画了之类的拼装错误不会变成一张和名字不符的图。
"""
import math
import re
from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from mzdesk.lipid_name_map import normalize

FLAG_DB, FLAG_FG, FLAG_SPLIT = 1, 2, 4
MAX_SUM_C, MAX_SUM_DB = 120, 30            # 拆总量的动态规划随 C、DB 变大而变慢；超过这个范围不画


class Unsupported(Exception):
    """这个名字画不出示意结构（链太短、记号不支持……）。"""


class NoStructureForClass(Unsupported):
    """这个类没有示意结构。"""


@dataclass
class Illustration:
    mol: Chem.Mol
    smiles: str
    formula: str
    db_atoms: list = field(default_factory=list)       # 默认 C=C 位置的原子
    fg_atoms: list = field(default_factory=list)       # 默认官能团位置的原子
    split_atoms: list = field(default_factory=list)    # 按默认规则拆出来的链的原子
    assumptions: list = field(default_factory=list)    # 简短英文，给页面的 Assumptions 用
    level: int = 0
    cls: str = ""


# ------------------------------------------------------------------ 默认位置
def default_db_positions(C, DB, kind):
    """一条 C 个碳、DB 个 C=C 的链的默认双键位置 [(位置, 'Z'|'E'), ...]。kind: acyl / alkyl / vinyl / base（鞘氨醇碱基）。"""
    if DB == 0:
        return []
    if kind == "vinyl":
        if DB == 1:
            return [(1, "Z")]
        rest = default_db_positions(C, DB - 1, "alkyl")
        rest = [(max(p, 4 + 3 * i), z) for i, (p, z) in enumerate(rest)]
        return [(1, "Z")] + rest
    if kind == "base":
        return [(4, "E")] + [(8 + 3 * k, "Z") for k in range(DB - 1)][:DB]
    lo = 2
    if DB == 1:
        p = 9 if C <= 18 else C - 9
        return [(max(lo, min(p, C - 2)), "Z")]
    anchor = 6 if DB in (2, 4) else 3
    pos = [C - anchor - 3 * k for k in range(DB)][::-1]
    if pos[0] < 3:
        for step in (3, 2):
            pos = [3 + step * k for k in range(DB)]
            if pos[-1] <= C - 2:
                break
        else:
            raise Unsupported("chain too short for %d C=C" % DB)
    return [(p, "Z") for p in pos]


# ------------------------------------------------------------------ 一条链的 SMILES 片段
def build_fragment(kind, C, dbs, subs, split):
    """从连接碳写到甲基端的片段。kind: acyl（C1 是 C(=O)）、alkyl（C1 接 O）、vinyl（C1=C2 接 O）、base（鞘氨醇，从 C4 起）。
    dbs {位置: Z/E}；subs {位置: OH / oxo / OOH / Ep / COOH / ('OX', 酰基片段)}；split：整条链标成"默认拆分"。"""
    fsplit = FLAG_SPLIT if split else 0
    start = {"acyl": 2, "alkyl": 1, "vinyl": 1, "base": 4}[kind]
    if any(p < start or p + 1 > C for p in dbs):
        raise Unsupported("C=C does not fit in a chain of %d carbons" % C)
    ep = {p for p, s in subs.items() if s == "Ep"}
    ep2 = {p + 1 for p in ep}
    # bonds[j] = 原子 j 和 j+1 之间的键符号；立体标记只在双键两端都有外侧邻居时写
    bonds = {p: "=" for p in dbs}
    for p, z in sorted(dbs.items()):
        has_before = p - 1 >= start - 1
        has_after = p + 2 <= C
        if kind == "acyl" and p == 1:
            continue
        if has_before and has_after and (p - 1) not in bonds and (p + 1) not in bonds:
            bonds[p - 1] = "/"
            bonds[p + 1] = "\\" if z == "Z" else "/"
    # 醚键旁的乙烯基双键（P- 链的 C1=C2）是名字规定的，不算猜的
    guessed = {p for p in dbs if not (kind == "vinyl" and p == 1)}
    out = []
    if kind == "acyl":
        out.append("[C:%d](=O)" % fsplit if fsplit else "C(=O)")
    for i in range(start, C + 1):
        s = subs.get(i)
        prev_bo = 2 if bonds.get(i - 1) == "=" else 1
        next_bo = 0 if i == C else (2 if bonds.get(i) == "=" else 1)
        flag = fsplit
        if i in guessed or (i - 1) in guessed:
            flag |= FLAG_DB
        branch, sub_bo = "", 0
        phyto_oh = s == "OH" and kind == "base" and i == 4     # 植物鞘氨醇的 4-OH 是碱基本身规定的，不是猜的
        if phyto_oh:
            branch, sub_bo = "(O)", 1
        elif s == "OH":
            branch, sub_bo, flag = "([OH:2])", 1, flag | FLAG_FG
        elif s == "oxo":
            branch, sub_bo, flag = "(=[O:2])", 2, flag | FLAG_FG
        elif s == "OOH":
            branch, sub_bo, flag = "(O[OH:2])", 1, flag | FLAG_FG
        elif isinstance(s, tuple) and s[0] == "OX":          # 酯化的酰基（FAHFA）
            branch, sub_bo, flag = "(O" + s[1] + ")", 1, flag | FLAG_FG
        elif s == "COOH":
            out.append(bonds.get(i - 1, "") + "[C:%d](=O)[OH:%d]" % (flag | FLAG_FG, flag | FLAG_FG))
            continue
        if i in ep:                                          # 环氧：两个相邻碳 + 一个 O
            flag |= FLAG_FG
            out.append(bonds.get(i - 1, "") + "[CH:%d]1" % flag)
            out.append("[O:2]")
            continue
        if i in ep2:
            out.append("[CH:%d]1" % (flag | FLAG_FG))
            continue
        h = 4 - prev_bo - next_bo - sub_bo
        if h < 0:
            raise Unsupported("valence at C%d" % i)
        hs = "" if h == 0 else ("H" if h == 1 else "H%d" % h)
        atom = ("[C%s:%d]" % (hs, flag)) if flag else "C"
        if phyto_oh:                                         # 天然构型 4R
            atom = "[C@H:%d]" % flag if flag else "[C@H]"
        out.append(bonds.get(i - 1, "") + atom + branch)
    return "".join(out)


def chain_fragment(kind, ch, split=False, extra_OH=0, oh_pos=None, fahfa_acyl=None):
    """一条链 ch = {C, DB, toks} -> SMILES 片段。"""
    C, DB = ch["C"], ch["DB"]
    toks = dict(ch.get("toks") or {})
    for k in toks:
        if k not in ("OH", "oxo", "OOH", "Ep", "COOH", "O"):
            raise Unsupported("chain modifier " + k)
    n_oh = toks.get("OH", 0) + toks.get("O", 0) + extra_OH
    n_oxo, n_ooh, n_ep, n_cooh = toks.get("oxo", 0), toks.get("OOH", 0), toks.get("Ep", 0), toks.get("COOH", 0)
    if n_cooh > 1 or n_ep > 1:
        raise Unsupported("more than one COOH / epoxide")
    # 数据集约定：标签里的 DB = C=C + oxo 个数 + 环氧个数 + COOH 个数
    cc = DB - n_ep - n_oxo - n_cooh
    if cc < 0:
        raise Unsupported("C=C count below the number of oxo / epoxide / COOH")
    dbp = default_db_positions(C, cc, kind)
    if kind == "base" and ch.get("t"):
        dbp = [(p, "Z") for p in (8, 11, 14, 17, 20)[:cc]]
    dbs = dict(dbp)
    taken = set(dbs) | {p + 1 for p in dbs}
    subs = {}
    first = 2 if kind != "base" else 4

    def free(p0):
        p = p0
        while p <= C - 1 and (p in taken or p in subs):
            p += 1
        return p

    if oh_pos:
        for p in oh_pos:
            if n_oh <= 0:
                break
            if p <= C - 1 and p not in taken and p not in subs:
                subs[p] = "OH"
                n_oh -= 1
    for group, n in (("OH", n_oh), ("oxo", n_oxo), ("OOH", n_ooh)):
        for _ in range(n):
            p = free(first)
            if p > C - 1:
                raise Unsupported("no room for a group")
            subs[p] = group
    if fahfa_acyl is not None:
        p = free(9 if C > 10 else first)
        if p > C - 1:
            raise Unsupported("no room for the ester")
        subs[p] = ("OX", fahfa_acyl)
    if n_cooh:
        subs[C] = "COOH"
    if n_ep:
        p = C - 3
        while p > 2 and (p in taken or p + 1 in taken or p in subs or p + 1 in subs):
            p -= 1
        if p <= 2:
            raise Unsupported("no room for the epoxide")
        subs[p] = "Ep"
    if kind == "acyl" and C < 2:
        raise Unsupported("acyl chain shorter than C2")
    return build_fragment(kind, C, dbs, subs, split)


# ------------------------------------------------------------------ 只知道总量时的默认拆分
# 常见链的权重（越大越典型）；表外的链：偶数碳优于奇数碳，DB 越多越罕见
COMMON = {(16, 0): 10, (18, 0): 9, (18, 1): 10, (16, 1): 6, (18, 2): 8, (20, 4): 8, (22, 6): 7, (20, 3): 5, (22, 5): 5,
          (14, 0): 5, (18, 3): 5, (20, 0): 4, (22, 0): 4, (24, 0): 4, (24, 1): 4, (22, 4): 4, (20, 5): 5, (20, 2): 3,
          (22, 1): 3, (20, 1): 3, (12, 0): 3, (15, 0): 3, (17, 0): 3, (17, 1): 2, (22, 2): 2, (14, 1): 2, (10, 0): 2, (8, 0): 2}


def chain_w(c, d):
    if (c, d) in COMMON:
        return math.log(COMMON[(c, d)])
    if d > max(0, (c - 3) // 3 + (1 if c > 6 else 0)):
        return -50.0
    w = 1.0 if c % 2 == 0 else 0.3
    return math.log(w) - 0.5 * d


def split_sum(tot_c, tot_db, n):
    """总量 (C, DB) 拆成 n 条链，取最典型的一组（对常见链表做动态规划）；返回按 (DB, C) 升序的 [(C, DB), ...]。"""
    if tot_c > MAX_SUM_C or tot_db > MAX_SUM_DB:
        raise Unsupported("sum composition too large to split")
    cand = [(c, d) for c in range(2, 45) for d in range(0, 13)]
    best = {(0, 0): (0.0, [])}
    for _ in range(n):
        nxt = {}
        for (c0, d0), (score, chains) in best.items():
            for c, d in cand:
                c1, d1 = c0 + c, d0 + d
                if c1 > tot_c or d1 > tot_db:
                    continue
                s = score + chain_w(c, d)
                if (c1, d1) not in nxt or s > nxt[(c1, d1)][0]:
                    nxt[(c1, d1)] = (s, chains + [(c, d)])
        best = nxt
    if (tot_c, tot_db) not in best:
        raise Unsupported("no default split")
    return sorted(best[(tot_c, tot_db)][1], key=lambda x: (x[1], x[0]))      # DB 少的、短的在前


def chains_of(r):
    return [dict(C=c["C"], DB=c["DB"], toks=dict(c["toks"]), ether=c.get("ether", "")) for c in r["meta"].get("chains") or []]


# ------------------------------------------------------------------ 头基
PL_HEAD = {
    "PC": "P(=O)([O-])OCC[N+](C)(C)C", "PE": "P(=O)(O)OCCN", "PS": "P(=O)(O)OC[C@H](N)C(=O)O",
    "PG": "P(=O)(O)OCC(O)CO", "PI": "P(=O)(O)OC1C(O)C(O)C(O)C(O)C1O", "PA": "P(=O)(O)O",
}
GLC = "O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
HEX2 = "OC1OC(CO)C(OC2OC(CO)C(O)C(O)C2O)C(O)C1O"
GLYCEROL_C2 = "[C@@H]"          # 天然 sn-甘油：磷脂的甘油 C2 是 (R)


def _total_o(chains):
    """链上已写明的氧个数（OH oxo Ep O 各 1，OOH COOH 各 2）。"""
    return sum(c["toks"].get("OH", 0) + c["toks"].get("oxo", 0) + c["toks"].get("Ep", 0) + 2 * c["toks"].get("OOH", 0)
               + 2 * c["toks"].get("COOH", 0) + c["toks"].get("O", 0) for c in chains)


def _distinct(chains):
    """不同的链有几种（链全相同时，sn 位置换不换都是同一个分子）。"""
    return len({(c["C"], c["DB"], tuple(sorted(c["toks"].items())), c["ether"]) for c in chains})


# ------------------------------------------------------------------ sn 顺序
_CHAIN_TOKEN = re.compile(r"(\d+):(\d+)")


def name_sn_order(raw):
    """名字里用 / 给出的链顺序 [(C, DB), ...]（0:0 占位也算）；没有用 / 连接的链（或不足两条）时为 None。"""
    text = str(raw).split("|")[0]
    toks = list(_CHAIN_TOKEN.finditer(text))
    if len(toks) < 2 or any("/" not in text[a.end():b.start()] for a, b in zip(toks, toks[1:])):
        return None
    return [(int(t.group(1)), int(t.group(2))) for t in toks]


def _sn_note(raw, drawn, assumed):
    """sn 位置的说明：名字用 / 给了顺序 -> 画的顺序和名字一致就不说，不一致说"名字里的顺序没用上"；没给顺序 -> assumed 为真时说"sn 位置是假设的"。
    drawn：画出来的 sn-1、sn-2、…… 位置上的 (C, DB)，空位 (0, 0)。"""
    given = name_sn_order(raw) if raw else None
    if given is None:
        return "sn positions assumed" if assumed else None
    drawn = list(drawn) + [(0, 0)] * (len(given) - len(drawn))
    return None if given == drawn[:len(given)] else "sn order in the name is not used"


# ------------------------------------------------------------------ 各类的拼装；每个函数返回 (SMILES, 假设列表, 是否用了默认拆分)
def _glycerolipid(cf, base, mods, chs, sm, assump, linkage="", raw=None):
    tot_c, tot_db, tot_o = sm.get("C"), sm.get("DB"), sm.get("O", 0)
    ether, lyso = "ether" in mods, "lyso" in mods
    if mods - {"ether", "lyso"}:
        raise NoStructureForClass("class " + cf)
    n = {"MG": 1, "DG": 2, "TG": 3}.get(base, 2)
    if base in PL_HEAD and lyso:
        n = 1
    plas = linkage == "P"                                    # 缩醛（P-）：醚链是 O-CH=CH-R，C1=C2 是名字规定的，标签里的 DB 含它
    split = False
    if chs is None:                                          # L1：只有总量，拆链
        split = True
        if ether:
            best = None
            ether_c = [(tot_c, tot_db)] if n == 1 else [(16, 0), (18, 0), (18, 1), (16, 1), (20, 0), (14, 0), (18, 2)]
            for ec in ether_c:
                if plas and n != 1:
                    ec = (ec[0], ec[1] + 1)
                if plas and ec[1] < 1:
                    continue
                if n == 1:
                    rest = []
                else:
                    try:
                        rest = split_sum(tot_c - ec[0], tot_db - ec[1], n - 1)
                    except Unsupported:
                        continue
                score = chain_w(ec[0], ec[1] - (1 if plas else 0)) + sum(chain_w(*x) for x in rest)
                if best is None or score > best[0]:
                    best = (score, ec, rest)
            if best is None:
                raise Unsupported("no default split for the ether chain")
            chs = [dict(C=best[1][0], DB=best[1][1], toks={}, ether="P" if plas else "O")] + [dict(C=c, DB=d, toks={}, ether="") for c, d in best[2]]
            if not plas:
                assump.append("Ether chain drawn as O-alkyl")
        else:
            chs = [dict(C=c, DB=d, toks={}, ether="") for c, d in split_sum(tot_c, tot_db, n)]
        if tot_o:                                            # 多出来的氧：放在最后一条（最不饱和的）酰基链上当 OH
            acyl_idx = [i for i, c in enumerate(chs) if not c["ether"]] or [0]
            chs[acyl_idx[-1]]["toks"] = {"OH": tot_o}
            assump.append("Untyped O drawn as OH")
    else:
        if len(chs) != n:
            raise Unsupported("chain count %d, class needs %d" % (len(chs), n))
        if ether and not any(c["ether"] for c in chs):       # L2 的醚类：不知道哪条是醚链，取最饱和、最短的（P- 的取带 C=C 里最饱和的）
            pool = [k for k in range(len(chs)) if chs[k]["DB"] >= 1] if plas else list(range(len(chs)))
            if not pool:
                raise Unsupported("no chain can carry the vinyl ether")
            i = min(pool, key=lambda k: (chs[k]["DB"], chs[k]["C"]))
            chs[i]["ether"] = "P" if plas else "O"
            split = True
            assump.append("Ether chain assumed to be the most saturated one" + ("" if plas else ", drawn as O-alkyl"))
        if (not ether) and any(c["ether"] for c in chs):
            raise Unsupported("ether chain in a non-ether class")
        if tot_o and _total_o(chs) < tot_o:                  # 链上没写明的氧 -> 酰基链上的 OH
            acyl_idx = [i for i, c in enumerate(chs) if not c["ether"]] or [0]
            last = chs[acyl_idx[-1]]
            last["toks"] = dict(last["toks"])
            last["toks"]["OH"] = last["toks"].get("OH", 0) + tot_o - _total_o(chs)
            assump.append("Untyped O drawn as OH")

    def frag(c):
        if c["ether"] == "P":
            return chain_fragment("vinyl", c, split)
        if c["ether"] == "O":
            return chain_fragment("alkyl", c, split)
        return chain_fragment("acyl", c, split)

    chs = sorted(chs, key=lambda c: (0 if c["ether"] else 1, c["DB"], c["C"]))        # 醚链在 sn-1；其余 DB 少的在前
    fr = [frag(c) for c in chs]
    head = "O" + PL_HEAD[base] if base in PL_HEAD else None
    c2 = GLYCEROL_C2
    if n == 1 and base in PL_HEAD:                           # 溶血磷脂：sn-1 链，sn-2 OH
        smi = "O(%s)C%s(O)C%s" % (fr[0], c2, head)
    elif base == "MG":
        smi = "O(%s)C%s(O)CO" % (fr[0], c2)
    elif base in PL_HEAD:
        smi = "O(%s)C%s(O%s)C%s" % (fr[0], c2, fr[1], head)
    elif base == "DG":
        smi = "O(%s)C%s(O%s)CO" % (fr[0], c2, fr[1])
    else:
        smi = "O(%s)C%s(O%s)CO%s" % (fr[0], c2, fr[1], fr[2])
    if base == "PI":
        assump.append("Inositol drawn without stereochemistry")
    elif base == "PG":
        assump.append("Head-group glycerol drawn without stereochemistry")
    acyl = [c for c in chs if not c["ether"]]                # 醚链（O- / P-）按约定在 sn-1，位置不是猜的
    assumed = (len(acyl) >= 2 and _distinct(acyl) > 1) or (n == 1 and not ether)
    drawn = [(c["C"], c["DB"] - (1 if c["ether"] == "P" else 0)) for c in chs]        # 名字里 P- 链的 DB 不含乙烯基 C=C
    drawn += [(0, 0)] * ({"MG": 3, "DG": 3, "TG": 3}.get(base, 2) - len(drawn))
    note = _sn_note(raw, drawn, assumed)
    if note:
        assump.append(note)
    return smi, assump, split


def _cardiolipin(chs, sm, assump, raw=None):
    split = chs is None
    if split:
        chs = [dict(C=c, DB=d, toks={}, ether="") for c, d in split_sum(sm.get("C"), sm.get("DB"), 4)]
    elif len(chs) != 4:
        raise Unsupported("CL needs 4 chains")
    tot_o = sm.get("O", 0)
    if tot_o > _total_o(chs):                                # 未分型的氧 -> 最后一条链上的 OH
        chs[-1]["toks"] = dict(chs[-1]["toks"])
        chs[-1]["toks"]["OH"] = chs[-1]["toks"].get("OH", 0) + tot_o - _total_o(chs)
        assump.append("Untyped O drawn as OH")
    fr = [chain_fragment("acyl", c, split) for c in chs]
    # 两个外侧甘油都是 sn-1,2-二酰基-sn-甘油-3-磷酸（都是 R；第二个从磷酸端写，标记要翻过来）；中间的甘油不标构型
    smi = "O(%s)C[C@@H](O%s)COP(=O)(O)OCC(O)COP(=O)(O)OC[C@H](O%s)CO%s" % tuple(fr)
    note = _sn_note(raw, [(c["C"], c["DB"]) for c in chs], _distinct(chs) > 1)
    if note:
        assump.append(note)
    return smi, assump, split


def _sphingolipid(cf, chs, sm, assump):
    tot_c, tot_db, tot_o = sm.get("C"), sm.get("DB"), sm.get("O", 0)
    has_acyl = cf != "SPB"
    split = False
    if chs is None:
        split = True
        if tot_o < 2:
            raise Unsupported("sphingolipid with fewer than 2 O")
        if has_acyl:
            bd = 1 if tot_db >= 1 else 0
            ac, ad = tot_c - 18, tot_db - bd
            if ac < 2:
                raise Unsupported("acyl chain too short")
            chs = [dict(C=18, DB=bd, toks={"OH": 2}, ether="B"), dict(C=ac, DB=ad, toks={"OH": tot_o - 2} if tot_o > 2 else {}, ether="")]
            assump.append("Sphingoid base assumed d18:%d, the rest as N-acyl chain" % bd)
        else:
            chs = [dict(C=tot_c, DB=tot_db, toks={"OH": tot_o}, ether="B")]
    else:
        if len(chs) != (2 if has_acyl else 1):
            raise Unsupported("sphingolipid chain count")
        if has_acyl:
            def base_score(c):
                return (c["toks"].get("OH", 0) >= 2) * 10 + (14 <= c["C"] <= 20) * 3 + (c["DB"] <= 2) * 3 + (c["DB"] >= 1) * 1 - c["C"] * 0.01
            if chs[0]["toks"].get("OH", 0) != chs[1]["toks"].get("OH", 0):
                ib = max(range(2), key=lambda k: chs[k]["toks"].get("OH", 0))          # 羟基多的那条是碱基
            else:
                ib = max(range(2), key=lambda k: base_score(chs[k]))
                split = True
                assump.append("Sphingoid base chain guessed")
            chs = [chs[ib], chs[1 - ib]]
            if chs[0]["toks"].get("OH", 0) < 2:              # L2：氧没写成 OH，碱基按 d 型（2 个 OH），多出来的氧放在酰基链上
                chs[0]["toks"] = dict(chs[0]["toks"])
                chs[0]["toks"]["OH"] = 2
                if tot_o - 2 > chs[1]["toks"].get("OH", 0):
                    chs[1]["toks"] = dict(chs[1]["toks"])
                    chs[1]["toks"]["OH"] = tot_o - 2
                assump.append("Sphingoid base drawn as d type (2 OH)")
        elif chs[0]["toks"].get("OH", 0) < 1:
            chs[0]["toks"] = {"OH": max(tot_o, 2)}
    b = chs[0]
    n_oh = b["toks"].get("OH", 0)
    if n_oh < 2 or n_oh > 3:
        raise Unsupported("sphingoid base with %d OH" % n_oh)
    bt = dict(b["toks"])
    bt.pop("OH", None)
    if n_oh == 3:                                            # 植物鞘氨醇：1,3-二醇在模板里，这里只加 4-OH
        bt["OH"] = 1
    fb = chain_fragment("base", dict(C=b["C"], DB=b["DB"], toks=bt, t=(n_oh == 3)), split)
    head = {"Cer": "O", "SM": "OP(=O)([O-])OCC[N+](C)(C)C", "SPB": "O", "Cer-Hex": GLC, "Cer-Hex2": HEX2}[cf]
    assump.append("Natural configuration at the sphingoid base")
    if cf == "Cer-Hex":
        assump.append("Hexose drawn as β-glucose")
    elif cf == "Cer-Hex2":
        assump.append("Sugars drawn without stereochemistry")
    if not has_acyl:
        return "N[C@@H](C%s)[C@H](O)%s" % (head, fb), assump, split
    fa = chain_fragment("acyl", chs[1], split, oh_pos=[2])
    return "N(%s)[C@@H](C%s)[C@H](O)%s" % (fa, head, fb), assump, split


def _single_chain(cf, chs, sm, assump):
    tot_c, tot_db, tot_o = sm.get("C"), sm.get("DB"), sm.get("O", 0)
    if chs is None:
        toks = {}
        if tot_o:
            toks = {"OH": tot_o}
            assump.append("Untyped O drawn as OH")
        chs, split = [dict(C=tot_c, DB=tot_db, toks=toks, ether="")], True
    else:
        split = False
        if len(chs) != 1:
            raise Unsupported("chain count")
        n_tok = _total_o(chs)
        if tot_o and n_tok < tot_o:
            chs[0]["toks"] = dict(chs[0]["toks"])
            chs[0]["toks"]["OH"] = chs[0]["toks"].get("OH", 0) + tot_o - n_tok
            assump.append("Untyped O drawn as OH")
    c = chs[0]
    if cf == "FOH":
        return "O" + chain_fragment("alkyl", c, split), assump, split
    f = chain_fragment("acyl", c, split)
    if cf == "FA":
        return "O" + f, assump, split
    if cf == "CAR":
        return "O(%s)C(CC(=O)[O-])C[N+](C)(C)C" % f, assump, split
    if cf == "FA-etn":
        return "N(%s)CCO" % f, assump, split
    return "N" + f, assump, split                            # FA-amide


# 侧链（按固醇的碳数、双键数）：27:1 胆固醇、27:2 去氢胆固醇、28:1 菜油甾醇、28:2（22E，C24 构型不定）、29:1 谷甾醇、29:2 豆甾醇；母核是胆甾-5-烯-3β-醇
STEROL_SIDE = {(27, 1): "[C@H](C)CCCC(C)C", (27, 2): "[C@H](C)CCC=C(C)C", (28, 1): "[C@H](C)CC[C@@H](C)C(C)C",
               (28, 2): "[C@H](C)/C=C/C(C)C(C)C", (29, 1): "[C@H](C)CC[C@@H](CC)C(C)C", (29, 2): "[C@H](C)/C=C/[C@@H](CC)C(C)C"}
STEROL_CORE = "O(%s)[C@H]1CC[C@@]2(C)C(=CC[C@H]3[C@@H]4CC[C@H](%s)[C@@]4(C)CC[C@@H]32)C1"        # 与 PubChem 的胆固醇逐个手性中心一致（tests 核对）


def _sterol_ester(chs, sm, assump):
    split = False
    if chs is None:
        chs = [dict(C=27, DB=1, toks={"OH": 1}, ether=""), dict(C=sm.get("C") - 27, DB=sm.get("DB") - 1, toks={}, ether="")]
        split = True
        assump.append("Sterol assumed cholesterol, the rest as acyl chain")
    if len(chs) != 2:
        raise Unsupported("sterol ester needs 2 chains")
    cands = [c for c in chs if 27 <= c["C"] <= 29 and 1 <= c["DB"] <= 2 and c["toks"].get("OH", 1) == 1 and not set(c["toks"]) - {"OH"}]
    if not cands:
        raise Unsupported("sterol part is not a C27-C29 sterol with one OH")
    # 两条链都像固醇（CE 28:2 = 28:2 酰基 + 胆固醇）：带 OH 的那条是固醇，再看 27:1（胆固醇），都不是取第一条
    with_oh = [c for c in cands if c["toks"].get("OH")]
    pool = with_oh if len(with_oh) == 1 else cands
    st = next((c for c in pool if (c["C"], c["DB"]) == (27, 1)), pool[0])
    acyl = [c for c in chs if c is not st][0]
    tot_o = sm.get("O", 0)
    if tot_o > 1 + _total_o([acyl]):                         # 固醇的 OH 是类隐含的 1 个氧；再多的未分型氧 -> 酰基链上的 OH
        acyl["toks"] = dict(acyl["toks"])
        acyl["toks"]["OH"] = acyl["toks"].get("OH", 0) + tot_o - 1 - _total_o([acyl])
        assump.append("Untyped O drawn as OH")
    f = chain_fragment("acyl", acyl, split)
    side = STEROL_SIDE[(st["C"], st["DB"])]
    smi = STEROL_CORE % (f, side)
    assump.append("Sterol skeleton from the carbon count (cholest-5-en-3-ol type)")
    return smi, assump, split


def _fahfa(chs, sm, assump):
    split = False
    if chs is None:
        chs = [dict(C=c, DB=d, toks={}, ether="") for c, d in split_sum(sm.get("C"), sm.get("DB"), 2)]
        split = True
    if len(chs) != 2:
        raise Unsupported("FAHFA needs 2 chains")
    hydroxy = [i for i, c in enumerate(chs) if c["toks"].get("OH", 0) >= 1]
    if not hydroxy:
        hydroxy = [max(range(2), key=lambda k: chs[k]["C"])]
        chs[hydroxy[0]]["toks"] = {"OH": 1}
        assump.append("Hydroxy fatty acid assumed to be the longer chain")
    if len(hydroxy) != 1:
        raise Unsupported("FAHFA hydroxy chain unclear")
    h, a = chs[hydroxy[0]], chs[1 - hydroxy[0]]
    tot_o = sm.get("O", 0)
    if tot_o > _total_o(chs):                                # 未分型的氧 -> 酰基链上的 OH（羟基脂肪酸链上的 OH 已经是酯的那个氧）
        a["toks"] = dict(a["toks"])
        a["toks"]["OH"] = a["toks"].get("OH", 0) + tot_o - _total_o(chs)
        assump.append("Untyped O drawn as OH")
    acyl = chain_fragment("acyl", a, split)
    ht = dict(h["toks"])
    ht["OH"] = ht.get("OH", 0) - 1
    if ht["OH"] <= 0:
        ht.pop("OH")
    hf = chain_fragment("acyl", dict(C=h["C"], DB=h["DB"], toks=ht), split, fahfa_acyl=acyl)
    return "O" + hf, assump + ["Ester position assumed at C9"], split


def build(r, raw=None):
    """normalize 结果 -> (带标记的 SMILES, 假设列表, 是否用了默认拆分)。raw 是用户输入的原文（只用来看名字有没有用 / 给 sn 顺序）。"""
    if r["status"] != "ok":
        raise Unsupported("name not parsed")
    cf = r["class_full"]
    base = cf.split("-")[0]
    mods = set(cf.split("-")[1:])
    sm = r["meta"].get("sum") or {}
    chs = chains_of(r) if r["level"] >= 2 else None
    chs = chs or None
    assump = []
    if base in PL_HEAD or base in ("DG", "TG", "MG"):
        return _glycerolipid(cf, base, mods, chs, sm, assump, r["meta"].get("ether_linkage", ""), raw)
    if cf == "CL":
        return _cardiolipin(chs, sm, assump, raw)
    if base in ("Cer", "SM", "SPB"):
        if cf not in ("Cer", "SM", "SPB", "Cer-Hex", "Cer-Hex2"):
            raise NoStructureForClass("class " + cf)
        return _sphingolipid(cf, chs, sm, assump)
    if cf in ("FA", "FOH", "CAR", "FA-etn", "FA-amide"):
        return _single_chain(cf, chs, sm, assump)
    if cf == "ST-oacyl":
        return _sterol_ester(chs, sm, assump)
    if cf == "FAHFA":
        return _fahfa(chs, sm, assump)
    raise NoStructureForClass("class " + cf)


# ------------------------------------------------------------------ 分子式核对
# 类常数：C=0、DB=0、O=0 的"空"类的元素个数 (C, H, N, O, P, S)。标签 (C, DB, O) 加上去：C += C，H += 2C - 2DB，O += O（每个 ;O 恰好加一个氧）。
# 数值由 LIPID+ 参考表按类拟合得到（每类 100 多到 600 个已知结构，拟合残差为 0），不是从本模块的画法推出来的；醚键和溶血各让 H +2、O -1。
CLASS_CONST = {
    "PC": (8, 16, 1, 8, 1, 0), "PE": (5, 10, 1, 8, 1, 0), "PS": (6, 10, 1, 10, 1, 0), "PG": (6, 11, 0, 10, 1, 0), "PI": (9, 15, 0, 13, 1, 0), "PA": (3, 5, 0, 8, 1, 0),
    "MG": (3, 6, 0, 4, 0, 0), "DG": (3, 4, 0, 5, 0, 0), "TG": (3, 2, 0, 6, 0, 0), "CL": (9, 14, 0, 17, 2, 0),
    "Cer": (0, 1, 1, 1, 0, 0), "SM": (5, 13, 2, 4, 1, 0), "SPB": (0, 3, 1, 0, 0, 0), "Cer-Hex": (6, 11, 1, 6, 0, 0), "Cer-Hex2": (12, 21, 1, 11, 0, 0),
    "FA": (0, 0, 0, 2, 0, 0), "FOH": (0, 2, 0, 1, 0, 0), "CAR": (7, 13, 1, 4, 0, 0), "FA-etn": (2, 5, 1, 2, 0, 0), "FA-amide": (0, 1, 1, 1, 0, 0),
    "ST-oacyl": (0, -8, 0, 1, 0, 0), "FAHFA": (0, -2, 0, 3, 0, 0),
}
_ELEMENTS = "CHNOPS"


def expected_elements(r):
    """normalize 结果 -> 按 LIPID+ 总量算出的元素个数 {C, H, N, O, P, S}；类没有常数时为 None。"""
    cf = r["class_full"]
    base, mods = cf.split("-")[0], cf.split("-")[1:]
    sm = r["meta"].get("sum") or {}
    if not {"C", "DB"} <= set(sm):
        return None
    if cf in CLASS_CONST:
        const, rest = CLASS_CONST[cf], []
    elif base in PL_HEAD or base in ("MG", "DG", "TG"):
        const, rest = CLASS_CONST[base], mods
    else:
        return None
    if set(rest) - {"ether", "lyso"} or ("lyso" in rest and base not in PL_HEAD):
        return None
    c0, h0, n0, o0, p0, s0 = const
    for _ in rest:                                           # ether、lyso 各 H +2、O -1
        h0, o0 = h0 + 2, o0 - 1
    C, DB, O = sm["C"], sm["DB"], sm.get("O", 0)
    return dict(zip(_ELEMENTS, (c0 + C, h0 + 2 * C - 2 * DB, n0, o0 + O, p0, s0)))


def drawn_elements(mol):
    """分子里的元素个数（含氢；非 CHNOPS 的元素也会出现，核对时就会对不上）。"""
    out = dict.fromkeys(_ELEMENTS, 0)
    for a in mol.GetAtoms():
        out[a.GetSymbol()] = out.get(a.GetSymbol(), 0) + 1
        out["H"] += a.GetTotalNumHs()
    return out


def illustrate(name_or_result, raw=None):
    """名字（或 normalize 的结果）-> Illustration；画不了抛 Unsupported / NoStructureForClass。raw：用户输入的原文（给 normalize 结果时才需要）。"""
    if isinstance(name_or_result, dict):
        r = name_or_result
    else:
        r, raw = normalize(name_or_result), raw or name_or_result
    smiles, assump, _ = build(r, raw)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise Unsupported("RDKit rejected the generated SMILES")
    flags = {FLAG_DB: [], FLAG_FG: [], FLAG_SPLIT: []}
    for a in mol.GetAtoms():
        m = a.GetAtomMapNum()
        if m:
            for bit in flags:
                if m & bit:
                    flags[bit].append(a.GetIdx())
            a.SetAtomMapNum(0)
    Chem.SanitizeMol(mol)
    want = expected_elements(r)                              # 出图前的核对：画出来的元素个数必须等于名字（LIPID+ 总量）说的
    if want is None or drawn_elements(mol) != want:
        raise Unsupported("the drawn formula does not match the name")
    assumptions = list(assump)                               # C=C、拆分的链在图下有色标，不在这里重复
    if flags[FLAG_FG]:
        assumptions.insert(0, "Group positions assumed")
    return Illustration(mol=mol, smiles=Chem.MolToSmiles(mol), formula=CalcMolFormula(mol), db_atoms=flags[FLAG_DB], fg_atoms=flags[FLAG_FG],
                        split_atoms=flags[FLAG_SPLIT], assumptions=_ordered(assumptions), level=r["level"], cls=r["class_full"])


def _ordered(items):
    """去重，保持顺序。"""
    return list(dict.fromkeys(items))


# 名字里写了位置时（normalize 会丢掉它们），图里画的是默认位置：页面要说一声
_POSITION_RE = re.compile(r"\d:\d+\(\d+[EZ]?(?:,\d+[EZ]?)*\)|;\d+(?:OH|OOH|Ep|oxo|COOH)\b")


def positions_in_name(raw):
    """名字里有没有写 C=C 位置（18:1(9Z)）或官能团位置（;12OH）。"""
    return bool(_POSITION_RE.search(str(raw)))
