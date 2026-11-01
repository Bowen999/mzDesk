"""M6 语义核心：``interpret(lx, cr, cfg, source=None, hi=None) -> Interp``（SYNTHESIS §2.2 / §3.4 / §4，AUDIT_A 的 D13-D20）。

流水线（每一步的规则都是**表驱动的通用规则**，家族例外走 ``Handler`` 钩子，见 ``families/__init__.py``）::

    select_hint(D13) -> Ox/O2x 前缀(D20) -> Ctx 构造 -> Handler.build -> 默认装配 -> Interp

Ctx 构造（``_build_ctx``）
    1. 醚种类合并（S2，``conventions.ether_kind``）：链标记 ``O- P- e`` / 原生 ``;P`` 与类标记 ``Plasmanyl- Plasmenyl- -O -P``；
       ``dO-`` / ``p`` -> ``unsupported_syntax``；不一致 -> ``class_chain_conflict``；类不允许醚 -> ``ether_not_allowed``；
       ``;P`` 落在 DB 为 0 的链上 -> ``malformed``（所有家族）。类是醚类时，``;P`` 贴在没有 ``O-`` 的链上（LIPID+ 解码输出 ``PE-ether 22:6_20:1;P``）
       当作 ``O-20:1;P``，flag ``ether_prefix_inferred_from_P``（``conventions.chain_mark(ether_class=…)``）。
       入口另有两道通用闸：嵌套 ``(FA c:d)`` / 括号 ``(O-c:d)`` 前缀只许 ``*-oacyl`` 类（否则 ``unsupported_syntax``）；相邻括号组
       ``(a)(b)`` 只许家族钩子 ``allows_adjacent_groups``（CL 对记法，D24）。
    2. 丢 ``0:0`` 占位（``placeholder_dropped``）；介质缩写链按 ``data/oxylipin.csv`` 展开（``oxylipin_expanded``；
       ``allow_oxylipin_abbrev=False`` -> ``oxylipin_disabled`` 后备码）；支链 ``iso/15me`` 只留 ``Me`` 零价 token（``branch_dropped``）。
    2b. LipidMatch 链后括号 ``(Ke) (CHO) (COOH)``（词法给 ``Atom.via``）：C=O 计入标签口径的 DB（``_work_chain`` 里 ``DB + lipidmatch_ox_db_extra``，
       ``18:1(Ke)`` -> ``18:2;oxo``）；``(CHO)`` 收成 ``oxo``，flag ``aldehyde_as_oxo``。
    3. 氧分配（§4.3，按 ``ClassRes.oxy_policy``）：
       ``typed_only``  已分型 token 进 L3；纯 ``;On`` -> ``tail_o``（L3 空）；链上已有分型 token 时 ``;O`` 记 ``O`` 型。
       ``sphingo_OH``  碱基前缀 ``d/t/m`` -> ``OH`` 2/3/1（S1）；其余链上的 ``;O/;nO`` -> ``OH``（打 ``oxygen_assumed_OH``，
                       ``Config.assume_sphingolipid_OH=False`` 时退回 ``typed_only``）。
       ``class_OH``    同 ``typed_only``；类隐含氧在装配时落成角色链（甾醇链 / 羟基脂肪酸链）上的一个 ``OH``
                       （``class_implied_OH``，``Config.class_implied_OH_in_L3=False`` 时进 ``tail_o``）。
    4. 类隐含氧（S8）：``implied_mode == "add"`` 加 ``implied_O``；``"max"`` 补到 ``implied_O``（名字氧已够则不加）；打 ``implied_O``。
    5. PNNL 类后缀 ``-d/-t/-m``（S4）：碱基羟基 2/3/1 个，碱基链位置未知（``pending`` 角色 ``base``），``ClassRes.pnnl_O`` 超出 2 的部分
       进 ``tail_o``（打 ``pnnl_suffix_O``）。
    6. 类带入链（``ClassRes.inject``，如 CE 的 ``27:1``）进 ``Ctx.injected``（``implied_sterol_chain``）。
    7. 尾部记号（``;(On)``、外层括号后的氧）是整分子氧 -> ``tail_o``，不挂链。
    8. MS-DIAL 规范化名（S16，``conventions.is_harmonized_name``）：``sphingo_OH`` 策略的类（含家族是 ``fa`` 的 ``SL``）打
       ``harmonized_name_oxygen_placement``——按氧策略判，不按家族判（原来只有 ``sl`` 钩子打，``SL`` 漏了）。
    9. 鞘脂缺氧数（C26，2026-10-20 改判）：``sphingo_OH`` 策略的类，第 1-5 步之后总氧仍为 0（名字一个氧都没给：无 ``;On``、无 ``d/t/m``、
       无 ``-d/-t/-m`` 后缀，类也不带隐含氧；LipidMatch ``SM(30:1)``、LipiDetective ``SPB 18:0``）-> 按 d 型碱基补 2 个羟基，做法同 PNNL
       后缀 ``-d``（``pending`` 角色 ``base``，碱基链按 ``role_candidates`` 定位），打 ``base_oxygen_assumed``。``Config.assume_base_oxygen=False``
       时不补，默认装配照旧拒绝 ``missing_oxygen_count``。

默认装配（``default_structs``，对应 SYNTHESIS §4.1-§4.2）
    * 链数核对（D16）：``n_groups + n_inject`` 必须 ∈ ``ClassRes.nchain``（空集 = 豁免），否则 ``chain_count``（不悄悄降级）。
      例外：名字只有**一个链组**而类需要更多链 -> 总量（``sum_mode``，``chains=None``，flag ``sum_only``）。
    * 单个带 ``d/t/m`` 的记号且 C<26 的总量 -> ``single_token_base``（C33）；**多链鞘脂类**（``sphingo_OH``）只给一个链 token、
      带显式氧数（``Cer 18:1;O2``、``SM 16:0;2O``、``Hex3Cer 25:4;2O``）且 C<26：按和式收（level 1，L1 照字面，flags ``sum_only`` +
      ``convention_unverified``；D27 撤销 D23——数据集原生 L1 串必须能回灌，约 180 个参考表和式名也靠它）；没写氧数的
      （``SM 30:1``、``Cer 15:1``）由 Ctx 构造第 9 步按 d 型补了氧，同样处理（C<26 的单 token 也打 ``convention_unverified``）；
      ``Config.assume_base_oxygen=False`` 时鞘脂（``sphingo_OH``）总氧为 0 -> ``missing_oxygen_count``（C26 原设计）；
      总氧 12 或 >13 -> ``o_not_expressible``（S15，对钩子返回的结构同样检查）。
    * 假设枚举（≤6 个，重复的合并）：醚链未知（类标记 / ``Ether*`` 无链标记）-> 每条非 ``N-`` 链各当一次醚链；
      ``P`` 种类同时是 ``+1`` 的落点（``p_landing_ambiguous``）；``pending`` 角色（碱基/羟基链/甾醇链）各枚举候选
      （``base_chain_unknown`` ``hfa_chain_unknown`` ``sterol_chain_unknown``）；超过 6 个 -> 退回总量。
      ``Config.pnnl_ether_heuristic`` 为真时醚/P 落点取 DB 最少（再碳最少）的链（不推荐）。
    * ``Interp.cls`` 是**最终类**（hint 后的类；``alts[*].class_full`` 与之一致，D13）。
"""
from __future__ import annotations

import itertools
import re
from typing import Final

from ..classes import class_info
from ..config import Config
from ..types import TOKEN_VALENCE, Chain, ClassRes, HeadInfo, Interp, Lexed, RawChain, Reject, Struct
from . import conventions as cv
from . import registry
from .families import BaseHandler
from .model import Ctx, WChain

__all__ = ["interpret", "default_structs", "totals", "MAX_ALTS", "SINGLE_TOKEN_BASE_MIN_C"]

#: ``Interp.alts`` 的上限
MAX_ALTS: Final = 6
#: 单个 ``d/t/m`` 记号的总量 C 低于它时是"只有碱基"，拒绝 ``single_token_base``（C33）
SINGLE_TOKEN_BASE_MIN_C: Final = 26

_OX_HEAD: Final = re.compile(r"O2?x[A-Z]")
_ROLE_FLAG: Final = {"base": "base_chain_unknown", "hfa": "hfa_chain_unknown", "sterol": "sterol_chain_unknown"}
_ROLE_OF_SLASH: Final = {"sterol_first": "sterol", "hfa_second": "hfa"}
_TAIL_TYPED: Final = frozenset({"OH", "oxo", "Ep", "COOH", "OOH", "OO"})


# ==========================================================================
# 入口
# ==========================================================================
def interpret(lx: Lexed, cr: ClassRes, cfg: Config, source: str | None = None,
              hi: HeadInfo | None = None) -> Interp:
    """语义核心，见模块文档。``hi`` 给了才会问家族 ``select_hint``（D13）。"""
    h = registry.handler(cr.family)
    flags0: set[str] = set()
    if hi is not None:
        new = h.select_hint(hi, lx, cr, cfg, source)
        if new is not None and new != cr:
            if new.lyso and not cr.lyso and lx.is_sum and cfg.single_token_is_lyso(source):
                flags0.add("source_hint_lyso")
            cr = new
            h = registry.handler(cr.family)
    if cv.adjacent_groups(lx) and not h.allows_adjacent_groups(cr):   # D24：``(a)(b)`` 相邻括号组只有 CL 对记法有定义
        raise Reject("malformed", f"adjacent parenthesised chain groups (a)(b) are not defined for {cr.class_full}")
    if (lx.nested_fa or lx.pre_acyl is not None) and not h.allows_nested_acyl(cr):
        # 嵌套 ``(FA c:d)`` / 括号 ``(O-c:d)`` 酰基前缀只有 *-oacyl 类有定义（``Handler.allows_nested_acyl``；``gl`` 另许 ``TG-oacyl2``）；其余类里它们会被当成
        # "又一条链"（``TG 16:0_18:1(FA 16:0)`` 数出 3 条链就被当普通 TG），必须拒绝（B2 请求；原来只有 GL 钩子查）。不用 class_not_in_vocab：头部本身是词表类
        raise Reject("unsupported_syntax", f"{cr.class_full}: an (FA c:d) / (O-c:d) acyl group is only defined for the *-oacyl classes")
    if _OX_HEAD.match(lx.head):                                       # D20：Ox/O2x 前缀已被 M2 剥掉的才会有氧标记
        raise Reject("unsupported_syntax", f"oxidised-lipid prefix without an oxygen mark: {lx.head!r}")
    native = lx.head == cr.class_full or (
        hi is not None and any(r.class_full == cr.class_full and "native" in r.dialects.split("|") for r in hi.rows))
    ctx = _build_ctx(lx, cr, cfg, source, h, native)
    ctx.flags |= flags0
    if ctx.cr.oxy_policy == "sphingo_OH" and cv.is_harmonized_name(lx, ctx.works, ctx.separators):
        ctx.flags.add("harmonized_name_oxygen_placement")             # S16：氧策略 sphingo_OH 的类都算（家族是 fa 的 SL 也一样）
    out = h.build(lx, cr, cfg, ctx)
    if out is None:
        structs = default_structs(ctx)
    elif isinstance(out, Struct):
        structs = [out]
    else:
        structs = list(out)
    return _finalize(ctx, structs)


def _finalize(ctx: Ctx, structs: list[Struct]) -> Interp:
    if not structs:
        raise ValueError("Handler.build returned no hypothesis")
    uniq: list[Struct] = []
    seen: set[tuple] = set()
    for st in structs:                                                 # 链序不同但多重集相同的假设合并（标签相同）
        key = _struct_key(st)
        if key not in seen:
            seen.add(key)
            uniq.append(st)
    if len(uniq) > MAX_ALTS:
        raise ValueError(f"{len(uniq)} hypotheses (max {MAX_ALTS})")
    cfs = {st.class_full for st in uniq}
    if len(cfs) != 1:
        raise ValueError(f"hypotheses disagree on the class: {sorted(cfs)}")
    cf = next(iter(cfs))
    cr = ctx.cr
    if cf != cr.class_full:                                            # 钩子改了类：Interp.cls 必须是最终类
        cr = class_info(cf).to_class_res(cr.alias_status, cr.inject)
    o = _struct_o(uniq[0])
    if o == 12 or o > 13:
        raise Reject("o_not_expressible", f"total oxygen {o} cannot be written in L1")
    flags = ctx.flags
    if ctx.cr.alias_status != "verified":
        flags.add("alias_inferred")
        if ctx.cr.alias_status == "unverified":
            flags.add("convention_unverified")
    return Interp(cls=cr, alts=tuple(uniq), flags=tuple(sorted(flags)), ether_linkage=ctx.ether_linkage)


def _struct_key(st: Struct) -> tuple:
    chains = None if st.chains is None else tuple(sorted(
        (c.C, c.DB, tuple(sorted(c.toks)), c.ether) for c in st.chains))
    return st.class_full, chains, st.tail_o, st.total


def _struct_o(st: Struct) -> int:
    if st.chains is None:
        return st.total[2] if st.total is not None else 0
    return cv.valence([t for c in st.chains for t in c.toks]) + st.tail_o


# ==========================================================================
# Ctx 构造
# ==========================================================================
def _work_chain(rc: RawChain, ctx: Ctx) -> WChain:
    cr, cfg, flags = ctx.cr, ctx.cfg, ctx.flags
    if rc.oxylipin:
        if not cfg.allow_oxylipin_abbrev:
            raise Reject("oxylipin_disabled", f"oxylipin abbreviation {rc.oxylipin} is disabled")
        c, db, toks = cv.expand_oxylipin(rc.oxylipin)
        flags.add("oxylipin_expanded")
        return WChain(c, db, toks, kind="main", raw=rc)
    if rc.branch == "methyl":                                         # (methyl) / methyl- ：甲基个数/位置不详，C 也不计入（C31）
        raise Reject("unsupported_syntax", f"unspecified methyl branch in {rc.raw!r}")
    if rc.base and cr.oxy_policy != "sphingo_OH":
        raise Reject("class_chain_conflict", f"sphingoid prefix {rc.base!r} in non-sphingolipid class {cr.class_full}")
    for a in rc.atoms:
        if a.kind.startswith("conj_") and not ctx.handler.allows_conjugate(cr):
            raise Reject("unsupported_chain_modifier", f"conjugate ;{a.kind[5:]} on {cr.class_full}")
    toks, tail = cv.atoms_to_tokens(rc.atoms, cr.oxy_policy, assume_oh=cfg.assume_sphingolipid_OH, base=rc.base)
    if cr.oxy_policy == "sphingo_OH" and cfg.assume_sphingolipid_OH:
        for a in rc.atoms:
            if a.kind == "O_count":
                flags.add("oxygen_assumed_OH")
                break
    mark, extra = cv.chain_mark(rc, ether_class=cr.ether)
    if rc.branch:
        flags.add("branch_dropped")
        toks.setdefault("Me", 1)
    db = rc.DB + cv.lipidmatch_ox_db_extra(rc.atoms)                  # LipidMatch (Ke)/(CHO)/(COOH)：C=O 计入标签口径的 DB（18:1(Ke) -> 18:2;oxo）
    if any(a.via == "CHO" for a in rc.atoms):
        flags.add("aldehyde_as_oxo")                                   # 醛基 CHO 收成 oxo（C=O，+O −2H）；醛 / 酮在标签里不分
    ctx.tail_o += tail
    if toks and cv.demote_overlimit(toks, db + extra):
        moved = cv.valence({k: v for k, v in toks.items() if k not in ("Me", "cPr")})
        toks = {k: v for k, v in toks.items() if k in ("Me", "cPr")}
        ctx.tail_o += moved
        flags.add("oxygen_untyped")
    role = "base" if rc.base else "acyl" if rc.n_acyl else ""
    # 输入显式写出的羟基（``;OH`` ``;OH2`` ``;(2OH)`` + ``d/t/m`` 换算）；``;O`` / ``;On`` 补成的 OH 不算（``Chain.oh_read``，LM 名的鞘脂碱基推断用）
    oh_read = min(sum(a.n for a in rc.atoms if a.kind in ("OH", "OH_pos")) + cv.base_oxygen(rc.base), toks.get("OH", 0))
    return WChain(rc.C, db, toks, mark=mark, db_extra=extra, role=role, kind="main", n_acyl=rc.n_acyl,
                  base=rc.base, raw=rc, oh_read=oh_read)


def _build_ctx(lx: Lexed, cr: ClassRes, cfg: Config, source: str | None, h: BaseHandler, native: bool = False) -> Ctx:
    ctx = Ctx(lx=lx, cr=cr, cfg=cfg, source=source, handler=h, native=native)
    flags = ctx.flags
    # 1 醚种类
    kind, kflags = cv.ether_kind(lx, cr)
    ctx.ether_kind = kind
    # linkage 状态（label_rules Q4）：P = 烯基（P- / -P / Plasmenyl- / ;P）；O = 类头明说烷基（Plasmanyl-）；
    # unspecified = 只写了 O- / e / -O / Ether*（不当成"已确认烷基"）
    ctx.ether_linkage = ("P" if kind == "P" else
                         "O" if kind and lx.head.lower().startswith("plasmanyl-") else "unspecified" if kind else "")
    flags.update(kflags)
    # 2 主链 / O-酰基链 / 类带入链
    for rc in lx.chains:
        if rc.placeholder:
            flags.add("placeholder_dropped")
            continue
        ctx.works.append(_work_chain(rc, ctx))
    extra_raw = ([lx.pre_acyl] if lx.pre_acyl is not None else []) + list(lx.nested_fa)
    for rc in extra_raw:
        ctx.extras.append(WChain(rc.C, rc.DB, {}, kind="oacyl", role="oacyl", raw=rc))
    for c, db in cv.parse_inject(cr.inject):
        role = "sterol" if cr.slash_role == "sterol_first" else "inject"
        ctx.injected.append(WChain(c, db, {}, kind="inject", role=role))
        flags.add("implied_sterol_chain")
    ctx.n_groups = len(ctx.works) + len(ctx.extras)
    ctx.n_inject = len(ctx.injected)
    ctx.separators = tuple(lx.seps)
    # 3 尾部记号：整分子氧
    for a in lx.tail_atoms:
        k = a.kind
        if k == "O_count":
            ctx.tail_o += a.n
        elif k.startswith("conj_"):
            if not h.allows_conjugate(cr):
                raise Reject("unsupported_chain_modifier", f"conjugate ;{k[5:]} on {cr.class_full}")
        elif k == "P":
            raise Reject("unsupported_syntax", "';P' marker outside a chain")
        elif k in _TAIL_TYPED:
            ctx.tail_o += TOKEN_VALENCE[k] * a.n
            flags.add("oxygen_untyped")
        elif k == "OH_pos" or k == "O_other":
            ctx.tail_o += a.n
            flags.add("oxygen_untyped")
    ctx.name_o = sum(cv.valence(w.toks) for w in ctx.all_chains()) + ctx.tail_o        # 名字自己写出的氧（D27 的 flag、CL 对记守卫用）
    # 4 PNNL 类后缀 -d/-t/-m
    suf = lx.feats.class_suffix
    if suf in ("-d", "-t", "-m"):
        if any(w.base for w in ctx.works):
            raise Reject("class_chain_conflict", f"class suffix {suf} together with an explicit sphingoid prefix")
        base_n = cv.base_oxygen(suf[1])
        rem = max(0, (cr.pnnl_O if cr.pnnl_O is not None else 2) - 2)
        flags.add("pnnl_suffix_O")
        if cr.oxy_policy == "sphingo_OH":
            ctx.pending.append(("base", base_n))
        else:
            ctx.tail_o += base_n
        ctx.tail_o += rem
    # 5 类隐含氧
    mode = cr.implied_mode
    if mode != "none":
        name_o = (sum(cv.valence(w.toks) for w in ctx.all_chains()) + ctx.tail_o + sum(n for _, n in ctx.pending))
        if mode == "add":
            extra = 0 if native else cr.implied_O                         # 原生写法的标签空间里已含类隐含氧（幂等）
        else:
            extra = max(0, cr.implied_O - name_o)
        if extra:
            flags.add("implied_O")
            slot = _ROLE_OF_SLASH.get(cr.slash_role)
            if cr.oxy_policy == "class_OH" and cfg.class_implied_OH_in_L3 and slot:
                ctx.pending.append((slot, extra))
                flags.add("class_implied_OH")
            else:
                ctx.tail_o += extra
    # 6 鞘脂缺氧数（C26 改判，模块文档 Ctx 构造第 9 步）：一个氧都没有 -> 按 d 型碱基补 2 个羟基，同 PNNL 后缀 -d
    if (cr.oxy_policy == "sphingo_OH" and cfg.assume_base_oxygen
            and sum(cv.valence(w.toks) for w in ctx.all_chains()) + ctx.tail_o + sum(n for _, n in ctx.pending) == 0):
        ctx.pending.append(("base", cv.base_oxygen("d")))
        flags.add("base_oxygen_assumed")
    # 7 总量还是拆链
    nominal = cr.nchain
    ctx.sum_mode = bool(nominal) and lx.is_sum and (ctx.n_groups + ctx.n_inject) < min(nominal)
    return ctx


# ==========================================================================
# 默认装配
# ==========================================================================
def totals(ctx: Ctx) -> tuple[int, int, int]:
    """``(C, DB, O)``：全部链（含 O-酰基链与类带入链）+ 待放置羟基 + ``tail_o``；DB 含 P 的 +1（只加一次）。"""
    allw = ctx.all_chains()
    c = sum(w.C for w in allw)
    db = sum(cv.label_db(w.DB, w.db_extra) for w in allw)
    if ctx.ether_kind == "P" and not any(w.mark for w in ctx.works):
        db += cv.P_DB_BONUS
    o = sum(cv.valence(w.toks) for w in allw) + ctx.tail_o + sum(n for _, n in ctx.pending)
    return c, db, o


def default_structs(ctx: Ctx) -> list[Struct]:
    """默认装配，见模块文档。"""
    cr, cfg = ctx.cr, ctx.cfg
    cf = cr.class_full
    nominal = cr.nchain
    works = ctx.works
    n_total = ctx.n_groups + ctx.n_inject
    if not ctx.sum_mode and nominal and n_total not in nominal:
        raise Reject("chain_count", f"{cf} takes {sorted(nominal)} chain(s), the name gives {n_total}")
    if ctx.sum_mode and len(works) == 1 and works[0].C < SINGLE_TOKEN_BASE_MIN_C:
        w0 = works[0]
        if w0.base:                                                        # C33：Cer d18:1
            raise Reject("single_token_base", f"a single {w0.base}-prefixed token {w0.C}:{w0.DB} "
                                              f"(C < {SINGLE_TOKEN_BASE_MIN_C}) is a base chain, not a total")
        if cr.oxy_policy == "sphingo_OH" and (ctx.name_o > 0 or "base_oxygen_assumed" in ctx.flags):
            ctx.flags.add("convention_unverified")                         # D27：多链鞘脂类单个 token + 显式氧数（Cer 18:1;O2）或补的氧（Cer 15:1）
                                                                           #     照字面按和式收（下面 sum_only）；碱基 / 总量二义性记在 flag 上
    c, db, o = totals(ctx)
    if cr.oxy_policy == "sphingo_OH" and o == 0:
        raise Reject("missing_oxygen_count", f"{cf} needs its oxygen count (d/t/m prefix or ;On)")
    if o == 12 or o > 13:
        raise Reject("o_not_expressible", f"total oxygen {o} cannot be written in L1")
    if ctx.sum_mode:
        ctx.flags.add("sum_only")
        return [Struct(cf, None, total=(c, db, o))]
    tail = ctx.tail_o

    axes: list[tuple[str, list[int]]] = []
    kind = ctx.ether_kind
    if kind and not any(w.mark for w in works):
        cands = [i for i, w in enumerate(works) if not w.n_acyl] or list(range(len(works)))
        if len(cands) > 1:
            ctx.flags.add("p_landing_ambiguous" if kind == "P" else "ether_chain_unknown")
            if cfg.pnnl_ether_heuristic:
                cands = [min(cands, key=lambda i: (works[i].DB, works[i].C, -i))]
        axes.append(("ether", cands))
    roles: dict[str, int] = {}
    for role, n in ctx.pending:
        roles[role] = roles.get(role, 0) + n
    placed: dict[str, int] = {}
    for role, n in roles.items():
        cands = ctx.role_candidates(role)
        if not cands:
            tail += n                                                      # 无处可放：退回未分型
            continue
        placed[role] = n
        if len(cands) > 1 and role in _ROLE_FLAG:
            ctx.flags.add(_ROLE_FLAG[role])
        axes.append((role, cands))
    size = 1
    for _, cands in axes:
        size *= len(cands)
    if size > MAX_ALTS:
        ctx.flags.add("sum_only")
        return [Struct(cf, None, total=(c, db, o))]
    if tail:
        ctx.flags.add("oxygen_untyped")
    out: list[Struct] = []
    for combo in itertools.product(*(cands for _, cands in axes)):
        assign = {name: idx for (name, _), idx in zip(axes, combo)}
        st = _build_struct(ctx, assign, placed, tail)
        if all(_struct_key(st) != _struct_key(prev) for prev in out):
            out.append(st)
    return out


def _build_struct(ctx: Ctx, assign: dict[str, int], placed: dict[str, int], tail: int) -> Struct:
    allw = ctx.all_chains()
    eth = assign.get("ether")
    add_oh: dict[int, int] = {}
    role_at: dict[int, str] = {}
    for role, n in placed.items():
        idx = assign[role]
        add_oh[idx] = add_oh.get(idx, 0) + n
        role_at[idx] = role
    chains: list[Chain] = []
    for i, w in enumerate(allw):
        toks = w.toks
        if i in add_oh:
            toks = dict(toks)
            toks["OH"] = toks.get("OH", 0) + add_oh[i]
        mark, extra = w.mark, w.db_extra
        if i == eth:
            mark, extra = ctx.ether_kind, (cv.P_DB_BONUS if ctx.ether_kind == "P" else 0)
        chains.append(Chain(w.C, cv.label_db(w.DB, extra), tuple(toks.items()), mark, w.role or role_at.get(i, ""), w.oh_read))
    return Struct(ctx.cr.class_full, tuple(chains), tail)
