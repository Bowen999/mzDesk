"""Separation 页的规则：由 logD、净电荷、TPSA、MW、立体中心数给出每种分离模式的推荐，以及对流动相 pH 的建议。

只有纯规则，不依赖模型权重。所有阈值都在 Thresholds 里（一处）。规则很粗，是起点不是结论：

  RP            logD >= 1 good；-1 <= logD < 1 possible；logD < -1 poor
  HILIC         logD < 0 good，或 TPSA >= 100 Å² 且 logD < 1 也算 good；logD > 2 poor；其余 possible
  ion-pair RP   分子在这个 pH 带电（完全不带电的比例 < 50 %）且 logD < 0 时 possible，否则 poor；离子对试剂不适合 MS
  mixed-mode    有永久电荷，或这个 pH 下最可能的净电荷 |q| >= 2 时 possible，否则 poor
  chiral        有 >= 1 个立体中心（RDKit 计数，未指定构型的也算；内消旋除外）：needed only for enantiomer separation，否则 not needed
  SEC           MW > 5000 good（蛋白、聚合物），否则 poor
  永久电荷      没有中性形式，logD 没有定义：RP 一律只给 possible（疏水部分还是可能保留，例如脂质），HILIC 当作 good，
                离子对当作带电且 logD < 0。
  pH 建议       每个位点在 pKa ± 2 之外 >99 % 是同一种形式（pKa ± 2 处少数形式 0.99 %）；所有位点都满足的 pH 范围就是稳健的流动相 pH。
"""
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Thresholds:
    rp_good: float = 1.0                  # RP: logD >= this -> good
    rp_poor: float = -1.0                 # RP: logD < this -> poor
    hilic_good: float = 0.0               # HILIC: logD < this -> good
    hilic_good_polar_logd: float = 1.0    # HILIC: ... or TPSA >= tpsa_high and logD < this -> good
    tpsa_high: float = 100.0              # Å²
    hilic_poor: float = 2.0               # HILIC: logD > this -> poor
    ion_pair_logd: float = 0.0            # ion-pair RP: charged and logD < this -> possible
    charged_below: float = 0.5            # "charged at this pH" = fraction of fully uncharged molecules below this
    mixed_mode_charge: int = 2            # mixed-mode / IEX: |dominant net charge| >= this -> possible
    sec_mw: float = 5000.0                # SEC: MW > this -> good (Da)
    window: float = 2.0                   # site is >99 % one form outside pKa ± window
    neutral_fraction: float = 0.9         # "neutral form dominates" for the pH hint
    max_site_notes: int = 10              # per-site pH lines shown (the sites closest to the chosen pH); the ranges use all sites


THRESHOLDS = Thresholds()
PH_MIN, PH_MAX = 0.0, 14.0
MINUS = "−"
ZWITTERION_NOTE = "zwitterion: logD likely strongly underestimated (zwitterions partition more than this model assumes)"


def charge_text(q):
    return "0" if q == 0 else (f"+{q}" if q > 0 else f"{MINUS}{-q}")


def logd_text(x):
    return "n/a" if x is None else (f"{x:.1f}".replace("-", MINUS))


def _num(x):
    """A threshold for display: -1 -> '−1', 0.5 -> '0.5'."""
    return f"{x:g}".replace("-", MINUS)


def mode_rows(*, ph, logd, tpsa, mw, charge, f_neutral, permanent, fixed_charge, n_stereo, n_unassigned=0, meso=False, th=THRESHOLDS):
    """-> [[mode, recommendation, reason], ...] in the order RP, HILIC, ion-pair RP, mixed-mode / ion exchange, chiral, SEC.

    logd is None when it is not defined (permanent charge). charge is the dominant net charge at ``ph``.
    """
    at = f"at pH {ph:g}"
    rows = []

    # ---- reversed phase
    if logd is None:
        rows.append(["RP", "possible", "permanently charged: no neutral form, logD is not defined; retention rests on the lipophilic part (ion-pairing or a low pH may help)"])
    elif logd >= th.rp_good:
        rows.append(["RP", "good", f"logD {logd_text(logd)} {at} is at least {_num(th.rp_good)}: well retained"])
    elif logd >= th.rp_poor:
        rows.append(["RP", "possible", f"logD {logd_text(logd)} {at} is between {_num(th.rp_poor)} and {_num(th.rp_good)}: weak retention"])
    else:
        rows.append(["RP", "poor", f"logD {logd_text(logd)} {at} is below {_num(th.rp_poor)}: elutes near the void volume"])

    # ---- HILIC
    if logd is None:
        rows.append(["HILIC", "good", "permanently charged: no neutral form, treated as polar"])
    elif logd < th.hilic_good:
        rows.append(["HILIC", "good", f"logD {logd_text(logd)} {at} is below {th.hilic_good:g}: polar"])
    elif tpsa >= th.tpsa_high and logd < th.hilic_good_polar_logd:
        rows.append(["HILIC", "good", f"TPSA {tpsa:.0f} Å² is at least {th.tpsa_high:g} and logD {logd_text(logd)} {at} is below {th.hilic_good_polar_logd:g}: polar"])
    elif logd > th.hilic_poor:
        rows.append(["HILIC", "poor", f"logD {logd_text(logd)} {at} is above {th.hilic_poor:g}: too lipophilic for HILIC retention"])
    else:
        rows.append(["HILIC", "possible", f"logD {logd_text(logd)} {at} is between {th.hilic_good:g} and {th.hilic_poor:g}; TPSA {tpsa:.0f} Å²"])

    # ---- ion-pair RP
    charged = permanent or f_neutral < th.charged_below
    if charged and (logd is None or logd < th.ion_pair_logd):
        rows.append(["Ion-pair RP", "possible", "charged " + at + (" and logD is not defined" if logd is None else f" with logD {logd_text(logd)} below {th.ion_pair_logd:g}")
                     + ": an ion-pair reagent adds retention, but ion-pair reagents are not MS-friendly"])
    elif not charged:
        rows.append(["Ion-pair RP", "poor", f"{100 * f_neutral:.0f} % of the molecules are uncharged {at}: nothing to pair"])
    else:
        rows.append(["Ion-pair RP", "poor", f"charged {at} but logD {logd_text(logd)} is at least {th.ion_pair_logd:g}: plain RP already retains it"])

    # ---- mixed-mode / ion exchange
    kind = "cation" if charge > 0 or (permanent and fixed_charge > 0) else "anion"
    if permanent:
        rows.append(["Mixed-mode / ion exchange", "possible", f"permanent charge ({charge_text(fixed_charge)}): {kind} exchange or mixed-mode"])
    elif abs(charge) >= th.mixed_mode_charge:
        rows.append(["Mixed-mode / ion exchange", "possible", f"net charge {charge_text(charge)} {at} (at least {th.mixed_mode_charge} in size): {kind} exchange or mixed-mode"])
    else:
        rows.append(["Mixed-mode / ion exchange", "poor", f"net charge {charge_text(charge)} {at}, not permanently charged: RP or HILIC is simpler"])

    # ---- chiral
    if n_stereo and not meso:
        extra = f" ({n_unassigned} not specified in the input)" if n_unassigned else ""
        rows.append(["Chiral", "needed only for enantiomer separation", f"{n_stereo} stereocentre{'s' if n_stereo > 1 else ''}{extra}"])
    elif n_stereo:
        rows.append(["Chiral", "not needed", f"{n_stereo} stereocentre{'s' if n_stereo > 1 else ''} but the molecule is meso (achiral)"])
    else:
        rows.append(["Chiral", "not needed", "no stereocentre"])

    # ---- SEC
    if mw > th.sec_mw:
        rows.append(["SEC", "good", f"MW {mw:.0f} Da is above {th.sec_mw:g}: proteins / polymers"])
    else:
        rows.append(["SEC", "poor", f"MW {mw:.0f} Da is not above {th.sec_mw:g}: SEC is for proteins and polymers"])
    return rows


def _intervals(mask, ph):
    """Consecutive True runs of ``mask`` on the pH grid -> [(lo, hi), ...]."""
    out, start = [], None
    for i, m in enumerate(mask):
        if m and start is None:
            start = i
        if start is not None and (not m or i == len(mask) - 1):
            end = i if m else i - 1
            out.append((float(ph[start]), float(ph[end])))
            start = None
    return out


def _range_text(ranges):
    return " or ".join(f"{lo:g}–{hi:g}" if hi > lo else f"{lo:g}" for lo, hi in ranges)


def site_window_text(site, th=THRESHOLDS):
    lo, hi = site.pka - th.window, site.pka + th.window
    below, above = ("neutral", "anion") if site.kind == "acid" else ("cation", "neutral")
    parts = []
    if lo > PH_MIN:
        parts.append(f"{below} at pH ≤ {lo:.1f}")
    if hi < PH_MAX:
        parts.append(f"{above} at pH ≥ {hi:.1f}")
    return f"{site.label} ({site.kind}, pKa {site.pka:.1f}): >99 % " + (" / ".join(parts) if parts else "one form over the whole range")


def robust_ph_ranges(sites, ph_grid, th=THRESHOLDS):
    """pH ranges in which every site is >99 % one form (outside pKa ± window for all sites)."""
    ok = np.ones(len(ph_grid), dtype=bool)
    for s in sites:
        ok &= (ph_grid <= s.pka - th.window + 1e-9) | (ph_grid >= s.pka + th.window - 1e-9)
    return _intervals(ok, ph_grid)


def neutral_ph_ranges(f_neutral, ph_grid, th=THRESHOLDS):
    return _intervals(f_neutral >= th.neutral_fraction, ph_grid)


def ph_notes(sites, ph, ph_grid, f_neutral, logp, permanent, th=THRESHOLDS):
    """Notes about the mobile-phase pH: the window of every site, the robust ranges, and whether ``ph`` is on a slope."""
    if not sites:
        return ["No ionizable site found: the molecule is uncharged at every pH, logD = logP, retention does not depend on pH."]
    shown = sorted(sorted(sites, key=lambda s: (abs(s.pka - ph), s.atom))[:th.max_site_notes], key=lambda s: (s.pka, s.atom))
    items = [site_window_text(s, th) for s in shown]
    if len(sites) > len(shown):
        items.append(f"{len(sites) - len(shown)} more sites not listed here (the ranges below use all {len(sites)}); see the table.")
    ranges = robust_ph_ranges(sites, ph_grid, th)
    if ranges:
        items.append(f"Robust mobile-phase pH (every site >99 % one form, pKa ± {th.window:g}): pH {_range_text(ranges)}.")
    else:
        items.append(f"No pH between {PH_MIN:g} and {PH_MAX:g} keeps every site >99 % one form: the pKa ± {th.window:g} windows cover the whole range; "
                     "retention will depend on pH however it is set.")
    near = [s for s in sites if abs(ph - s.pka) < th.window]
    if near:
        items.append(f"pH {ph:g} is within ± {th.window:g} of the pKa of " + ", ".join(f"{s.label} ({s.pka:.1f})" for s in near)
                     + ": a small pH shift changes the charge and the retention, so buffer well (or move away from these pKa values).")
    if not permanent:
        neutral = neutral_ph_ranges(f_neutral, ph_grid, th)
        if neutral:
            top = logp + float(np.log10(f_neutral.max()))
            items.append(f"The uncharged form is at least {100 * th.neutral_fraction:.0f} % of the molecules at pH {_range_text(neutral)} (logD up to {logd_text(top)}).")
        else:
            items.append(f"The fully uncharged form never reaches {100 * th.neutral_fraction:.0f} % at any pH (zwitterion or several sites close together): logD stays below logP.")
    return items
