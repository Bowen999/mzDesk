"""Known methods: published chromatography methods for a compound, plus the methods of the most similar compounds.

Sources: Japanese Pharmacopoeia 18 (bundled, jp18.py) and RepoRT (downloaded once, report.py / report_setup.py).
1. Exact: same connectivity key (first 14 characters of the InChIKey after largest fragment + neutralisation).
2. Similar (always computed, listed after the exact hits): score = w_structure * Tanimoto(Morgan r2, 2048 bits) + w_properties * mean(property terms);
   the property terms are exp(-|dlogP|/s), exp(-|dTPSA|/s), exp(-|dMW|/s) and the agreement of the charge class at pH 2.7 and 7.4
   (compound_index.Config holds every weight, scale and cut-off).
``search`` returns the structured result; ``page_outputs`` turns it into the blocks of the page (summary cards, exact-match table, grouped table of
the similar compounds; everything long sits behind "?" popovers, builders in mzdesk/api/blocks.py).
"""
import re
from dataclasses import dataclass, field

from mzdesk.api import blocks as bl
from mzdesk.methods import compound_index as ci
from mzdesk.methods import jp18, report

SOURCE_JP = "Japanese Pharmacopoeia"
SOURCE_REPORT = "RepoRT"
SOURCES = (SOURCE_JP, SOURCE_REPORT)
MISSING = "—"


class QueryError(ValueError):
    """The query structure cannot be used."""


class IndexProblem(RuntimeError):
    """An index is stale or unreadable (not the user's input)."""


@dataclass
class Method:
    """One method (a row of the exact-match table or of a similar compound's group)."""
    source: str                    # "JP18" | "RepoRT"
    compound: str
    match: str                     # "exact · same stereo" | "exact · stereo differs / not specified" | "similar"
    similarity: float | None       # None for exact hits shown as 1.00
    test: str
    column: str
    mobile: str
    detector: str
    rt: str | None
    reference: tuple               # (text, href or None)
    block_id: str | None = None    # JP only
    condition: str | None = None   # JP only: printed conditions
    label: str = ""                # text of the "Test / dataset" cell: JP "Purity (4) · p. 623", RepoRT "0382 AA_LC1 (RP)"
    monograph: str | None = None   # JP only: shown before the label when it is not the compound's own monograph
    help: list = field(default_factory=list)       # test, reference and conditions / gradient lines for the row's "?"


@dataclass
class SimilarCompound:
    name: str
    smiles: str
    key14: str
    score: float
    tanimoto: float
    d_logp: float
    d_tpsa: float
    d_mw: float
    q_low: int
    q_neutral: int
    sources: list
    methods: list = field(default_factory=list)


@dataclass
class Result:
    query: ci.Query
    exact: list
    similar: list
    searched: list                 # sources searched
    status: dict                   # source -> short description
    warnings: list
    cfg: ci.Config
    report_installed: bool | None = None       # None when RepoRT was not asked for
    report_commit: str = ""
    n_datasets: int = 0
    n_jp: int = 0

    @property
    def match(self):
        return "exact" if self.exact else "similar" if self.similar else "none"

    @property
    def methods(self):
        return self.exact + [m for s in self.similar for m in s.methods]


# ---------------------------------------------------------------- method rows
def _jp_method(jp, group, name, match, similarity, stereo_block=None, ik2_of=None, role_of=None):
    """One JP18 row for a group of records with identical conditions. Match: "exact · same stereo", "similar" ...; when the compound is only the
    monograph's substance and the test measures something else, " · monograph substance (test measures: ...)" is added."""
    r = jp.records[group[0]]
    text = jp18.test_text(r)
    if len(group) > 1:
        text += f" (+{len(group) - 1} more with identical conditions)"
    if stereo_block is not None:
        match = f"exact · {'same stereo' if ik2_of(group[0]) == stereo_block else 'stereo differs / not specified'}"
    if role_of and role_of(group[0]):
        head, _, tail = match.partition(" · ")
        match = f"{head} · monograph substance (test measures: {jp18.measures_text(r)})" + (f" · {tail}" if tail else "")
    reference = jp18.reference(r)
    more = f" (+{len(group) - 1})" if len(group) > 1 else ""
    page = r.get("page_printed") or r.get("page")
    label = _short_test(r) + (f" · p. {page}" if page else "") + more
    block_id = r["block_id"] + more
    condition = jp18.condition_text(r, jp.by_block)
    return Method(source="JP18", compound=name, match=match, similarity=similarity, test=text, column=jp18.column_text(r),
                  mobile=jp18.mobile_phase_text(r), detector=jp18.detector_text(r), rt=None, reference=reference, block_id=block_id, condition=condition,
                  label=label, monograph=r.get("monograph") or None,
                  help=[f"Test: {text}", f"Reference: {reference[0]}, block {block_id}", f"Printed conditions: {condition}"])


def _report_method(rep, ds_idx, rts, name, match, similarity, stereo_block=None, ik2s=()):
    ds = rep.datasets[ds_idx]
    if stereo_block is not None:
        match = f"exact · {'same stereo' if stereo_block in ik2s else 'stereo differs / not specified'}"
    reference = report.reference(ds)
    gradient = _gradient_points(ds)
    return Method(source="RepoRT", compound=name, match=match, similarity=similarity, test=report.dataset_title(ds),
                  column=report.column_text(ds), mobile=report.mobile_phase_text(ds), detector="MS", rt=ci.rt_text(rts),
                  reference=reference, label=report.dataset_title(ds),
                  help=[f"Reference: {reference[0]}"] + ([gradient] if gradient else []))


def _short_test(r):
    """"Purity (4)" from "Purity (4) 1,1′-[3,3′-dithiobis...]-L-diproline"; a short label is kept, a long one is cut."""
    t = r.get("test_label") or r.get("entry") or "Operating conditions"
    m = re.match(r"^(.*?\(\d+\))\s", t)
    if m:
        return m.group(1)
    return t if len(t) <= 40 else t[:39].rstrip() + "…"


def _gradient_points(ds):
    """"Gradient, % B: 0 min 5, 2 min 5, 13 min 95" from the dataset's gradient table (the varying eluent; B alone when there are two); "" when there is none."""
    g = ds.get("gradient") or {}
    t = g.get("t") or []
    letters = [L for L in "ABCD" if g.get(L)]
    if len(t) < 2 or not letters:
        return ""
    varying = [L for L in letters if max(g[L]) - min(g[L]) > 1e-9]
    show = "B" if varying in (["A", "B"], ["B"]) or not varying else varying[0]
    pts = [f"{x:g} min {round(v, 1):g}" for x, v in zip(t, g[show])]
    return f"Gradient, % {show}: " + ", ".join(pts[:12]) + (" …" if len(pts) > 12 else "")


def _jp_groups(jp, comp):
    """(groups of records with identical conditions, record -> stereo block, record -> role) of a compound."""
    refs, _rts, ik2 = jp.index.entries(comp)
    ids = [int(x) for x in refs]
    ik2_of = {int(r): bytes(k).decode() for r, k in zip(refs, ik2)}
    role = {int(r): int(x) for r, x in zip(refs, jp.index.roles(comp))}
    return jp18.groups(jp.records, ids, role), (lambda i: ik2_of[i]), (lambda i: role[i])


def _report_datasets(rep, comp):
    """{dataset index: ([rt, ...], {stereo blocks})} of a compound, best documented dataset first."""
    refs, rts, ik2 = rep.index.entries(comp)
    out = {}
    for r, t, k in zip(refs, rts, ik2):
        rec = out.setdefault(int(r), ([], set()))
        rec[0].append(float(t))
        rec[1].add(bytes(k).decode())
    return dict(sorted(out.items(), key=lambda kv: (report.metadata_gaps(rep.datasets[kv[0]]), kv[0])))


# ---------------------------------------------------------------- search
def search(smiles, sources=SOURCES, cfg=ci.CONFIG, report_path=None, jp_dir=None):
    """Methods for the compound (exact) and for the most similar compounds. Raises QueryError for a structure that cannot be used."""
    try:
        q = ci.make_query(smiles, cfg)
    except ValueError as exc:
        raise QueryError(str(exc)) from None
    warnings = [f"Structure: {n}." for n in q.norm.notes]
    if q.method != ci.METHOD_MOLGPKA:
        warnings.append("MolGpKa weights are not installed (python -m mzdesk.platform_fit.pka_setup): the charge class of the query comes from "
                        "functional groups (carboxylic / sulfonic / phosphoric acid: negative, aliphatic amine / amidine / guanidine / quaternary N: positive).")
    sources = [s for s in SOURCES if s in set(sources)]
    res = Result(query=q, exact=[], similar=[], searched=[], status={}, warnings=warnings, cfg=cfg)
    jp = rep = None
    if SOURCE_JP in sources:
        jp = jp18.load(cache_dir=jp_dir)
        _check(jp.index, cfg, "JP18")
        res.searched.append(SOURCE_JP)
        res.n_jp = len(jp.records)
        res.status[SOURCE_JP] = f"{len(jp.records)} LC / GC methods"
    if SOURCE_REPORT in sources:
        rep = report.load(report_path)
        res.report_installed = rep is not None
        if rep is None:
            res.status[SOURCE_REPORT] = f"not installed: run {report.SETUP_COMMAND}"
            warnings.append(f"RepoRT is not installed, so it was not searched: run {report.SETUP_COMMAND} (one-time download of about 100 MB).")
        else:
            _check(rep.index, cfg, "RepoRT")
            res.searched.append(SOURCE_REPORT)
            res.n_datasets, res.report_commit = len(rep.datasets), rep.commit
            res.status[SOURCE_REPORT] = f"{len(rep.datasets)} datasets, commit {rep.commit[:7]}"
    for lib in (jp, rep):
        if lib is not None and lib.index.charge_method != q.method:
            warnings.append(f"Charge classes: the query used {q.method}, the {lib.index.meta.get('source')} index {lib.index.charge_method}"
                            + (f" (rebuild with {report.SETUP_COMMAND} --force after installing MolGpKa)" if lib is rep else "") + ".")
    _exact(res, q, jp, rep)
    _similar(res, q, jp, rep)
    return res


def _check(index, cfg, label):
    try:
        index.check(cfg)
    except ValueError as exc:
        raise IndexProblem(f"{label}: {exc}") from None


def _exact(res, q, jp, rep):
    cfg = res.cfg
    stereo = q.norm.stereo_block
    if jp is not None:
        for comp in jp.index.exact(q.norm.key14):
            groups, ik2_of, role_of = _jp_groups(jp, comp)
            name = jp.index.names[comp] or q.norm.smiles
            res.exact += [_jp_method(jp, g, name, "exact", None, stereo, ik2_of, role_of) for g in groups[:cfg.max_exact_per_source]]
            if len(groups) > cfg.max_exact_per_source:
                res.warnings.append(f"JP18: {len(groups)} distinct conditions for this compound; the first {cfg.max_exact_per_source} are shown.")
    if rep is not None:
        for comp in rep.index.exact(q.norm.key14):
            name = rep.index.names[comp] or q.norm.smiles
            datasets = list(_report_datasets(rep, comp).items())
            for ds_idx, (rts, ik2s) in datasets[:cfg.max_exact_per_source]:
                res.exact.append(_report_method(rep, ds_idx, rts, name, "exact", None, stereo, ik2s))
            if len(datasets) > cfg.max_exact_per_source:
                res.warnings.append(f"RepoRT: this compound is in {len(datasets)} datasets; the {cfg.max_exact_per_source} best documented are shown.")


def _ranked(libs, q, cfg):
    """Similar compounds over all sources, best first, one at a time: yields (key14, {source: compound row}, stats).

    Each source is scored with one vectorised pass; the sources are merged lazily by score, and a compound found in both is yielded once
    (the compound itself is left out).
    """
    per = []
    for name, lib in libs:
        scores = lib.index.score(q, cfg)
        per.append((name, lib, scores, lib.index.candidates(scores, q.norm.key14, cfg)))
    ptr = [0] * len(per)
    seen = set()
    while True:
        best = None
        for i, (_name, lib, scores, rows) in enumerate(per):
            if ptr[i] < len(rows):
                c = rows[ptr[i]]
                order = (-float(scores.total[c]), -float(scores.structure[c]), bytes(lib.index.ik14[c]))
                if best is None or order < best[0]:
                    best = (order, i, int(c))
        if best is None:
            return
        order, i, c = best
        ptr[i] += 1
        key = order[2].decode()
        if key in seen:
            continue
        seen.add(key)
        name, lib, scores, _rows = per[i]
        parts = {n: int(hit[0]) for n, other in libs if len(hit := other.index.exact(key))}
        yield key, parts, {"score": -order[0], "tanimoto": float(scores.structure[c]), "d_logp": float(lib.index.logp[c] - q.logp),
                           "d_tpsa": float(lib.index.tpsa[c] - q.tpsa), "d_mw": float(lib.index.mw[c] - q.mw),
                           "q_low": int(lib.index.q_low[c]), "q_neutral": int(lib.index.q_neutral[c])}


def _similar(res, q, jp, rep):
    cfg = res.cfg
    libs = [(n, lib) for n, lib in ((SOURCE_JP, jp), (SOURCE_REPORT, rep)) if lib is not None]
    used = {}                                                    # RepoRT dataset index -> similar hits already listed from it
    for key, parts, c in _ranked(libs, q, cfg):
        methods, names = [], []
        label = None
        if SOURCE_JP in parts:
            comp = parts[SOURCE_JP]
            names.append(jp.index.names[comp])
            label = label or jp.index.smiles[comp]
            groups, _ik2_of, role_of = _jp_groups(jp, comp)
            methods += [_jp_method(jp, g, "", "similar", c["score"], role_of=role_of) for g in groups]
        report_hits = []
        if SOURCE_REPORT in parts:
            comp = parts[SOURCE_REPORT]
            names.append(rep.index.names[comp])
            label = label or rep.index.smiles[comp]
            report_hits = sorted(((ds, rts) for ds, (rts, _k) in _report_datasets(rep, comp).items()),
                                 key=lambda h: (used.get(h[0], 0), report.metadata_gaps(rep.datasets[h[0]]), h[0]))
        name = next((n for n in names if n and not ci.is_id_like(n)), None) or next((n for n in names if n), None) or label
        for m in methods:
            m.compound = name
        methods = methods[:cfg.max_methods_per_similar]
        for ds, rts in report_hits:
            if len(methods) >= cfg.max_methods_per_similar:
                break
            if used.get(ds, 0) >= cfg.max_similar_per_dataset:
                continue
            used[ds] = used.get(ds, 0) + 1
            methods.append(_report_method(rep, ds, rts, name, "similar", c["score"]))
        if not methods:
            continue
        res.similar.append(SimilarCompound(name=name, smiles=label, key14=key, score=c["score"], tanimoto=c["tanimoto"], d_logp=c["d_logp"],
                                           d_tpsa=c["d_tpsa"], d_mw=c["d_mw"], q_low=c["q_low"], q_neutral=c["q_neutral"],
                                           sources=[s for s in SOURCES if s in parts], methods=methods))
        if len(res.similar) >= cfg.max_similar:
            break


# ---------------------------------------------------------------- page blocks
def _q(c):
    return ci.CHARGE_NAMES[int(c)]


def _signed(v, nd):
    r = round(float(v), nd)
    return f"{'+' if r >= 0 else '−'}{abs(r):.{nd}f}"


def _where(sources):
    return " + ".join("JP18" if x == SOURCE_JP else x for x in sources)


def _match_notes(m):
    """The parts of the match string that matter in the row's "?": monograph substance and a stereo difference (same stereo is the default and not said)."""
    notes = []
    for part in m.match.split(" · ")[1:]:
        if part == "same stereo":
            continue
        measures = re.search(r"\(test measures: (.*)\)$", part)
        notes.append(f"Monograph substance: the test measures {measures.group(1)}." if part.startswith("monograph substance") and measures
                     else part[:1].upper() + part[1:] + ".")
    return notes


_ELUENT = re.compile(r"^[A-D]: ")


def _phase_cell(mobile):
    """Mobile phase for the table: one line per eluent ("A: …", "B: …"), then the rest (gradient · flow, or JP18's
    "isocratic · methanol · pH 3.0") on one line. The page shows each line separately."""
    parts = [p.strip() for p in mobile.split("; ") if p.strip()]
    eluents = [p for p in parts if _ELUENT.match(p)]
    rest = " · ".join(p for p in parts if not _ELUENT.match(p))
    return "\n".join(eluents + ([rest] if rest else [])) or mobile


def _method_row(m, compound):
    """Source badge, Test / dataset (linked to the JP page or the DOI), Column, Mobile phase, Detector, RT; the row's "?" has the rest."""
    text = m.label
    if m.monograph and m.monograph.strip().lower() != (compound or "").strip().lower():
        text = f"{m.monograph}: {text}"
    badge = bl.badge("JP18", "jp") if m.source == "JP18" else bl.badge(m.source, "report")
    return bl.row([badge, bl.link(text, m.reference[1]), m.column, _phase_cell(m.mobile), m.detector, m.rt or MISSING], help=_match_notes(m) + m.help)


def _title(name):
    return name if len(name) <= 60 else name[:59] + "…"


def page_outputs(res):
    """Result -> {block id: data} for the Known methods page: summary cards, the exact matches and the similar compounds (one collapsible group each)."""
    cfg = res.cfg
    methods = res.methods
    n_exact, n_sim = len(res.exact), len(methods) - len(res.exact)
    if res.exact:
        same = any("same stereo" in m.match for m in res.exact)
        note = "same stereo" if same else "stereo differs / not specified"
        if all("monograph substance" in m.match for m in res.exact):
            note = "monograph substance only; " + note
        match = {"value": "exact", "note": note}
    elif res.similar:
        best = res.similar[0]
        match = {"value": "similar", "note": f"best: {best.name} ({best.score:.2f})"}
    else:
        match = {"value": "none", "note": f"no compound with score ≥ {cfg.min_score:g} and Tanimoto ≥ {cfg.min_structure:g}"}
    searched = ", ".join(("JP18" if s == SOURCE_JP else s) for s in res.searched) or MISSING
    sources_note = f"RepoRT not installed: run {report.SETUP_COMMAND}" if res.report_installed is False else ""
    groups = [bl.group(_title(s.name), [_method_row(m, s.name) for m in s.methods], meta=f"{s.score:.2f} · {len(s.methods)} method{'' if len(s.methods) == 1 else 's'}",
                       help=_similar_help(s, res)) for s in res.similar]
    return {
        "summary": {"values": {
            "Match": match,
            "Methods found": {"value": str(len(methods)), "note": f"{n_exact} exact · {n_sim} for {len(res.similar)} similar compounds"},
            "Sources": {"value": searched, **({"note": sources_note} if sources_note else {}), "help": _sources_help(res)}}},
        "exact": bl.table_data([_method_row(m, m.compound) for m in res.exact], count=n_exact),
        "similar": bl.grouptable_data(groups, count=len(res.similar), help=_score_help(res)),
    }


def _score_help(res):
    """How the similarity score is made, the query's own properties and the cut-offs (the "?" of the Similar compounds block)."""
    q, cfg = res.query, res.cfg
    items = [
        f"Score = {cfg.w_structure:g} × Tanimoto (Morgan radius {cfg.fp_radius}, {cfg.fp_bits} bits) + {cfg.w_properties:g} × mean of five property terms: "
        f"exp(−|ΔlogP| / {cfg.scale_logp:g}), exp(−|ΔTPSA| / {cfg.scale_tpsa:g}), exp(−|ΔMW| / {cfg.scale_mw:g}), and 1 / 0 for the same / a different charge class "
        f"(sign of the dominant net charge) at pH {cfg.ph_acidic:g} and pH {cfg.ph_neutral:g}. Δ = library compound − query. "
        f"Listed: score ≥ {cfg.min_score:g} and Tanimoto ≥ {cfg.min_structure:g} (whatever the property score), at most {cfg.max_similar} compounds, {cfg.max_methods_per_similar} methods per compound, "
        f"{cfg.max_similar_per_dataset} per RepoRT dataset.",
        f"Query: logP {q.logp:.2f} (Crippen), TPSA {q.tpsa:.1f} Å², MW {q.mw:.1f}; charge class {_q(q.q_low)} at pH {cfg.ph_acidic:g}, {_q(q.q_neutral)} at pH {cfg.ph_neutral:g} "
        f"({q.method}; molecules are reduced to the largest fragment and neutralised first).",
    ]
    if not res.similar:
        items.append(f"No compound reached the score cut-off of {cfg.min_score:g} (with Tanimoto ≥ {cfg.min_structure:g}) in the sources searched.")
    return items


def _similar_help(s, res):
    """Why one compound is similar: Tanimoto, the property differences and the charge classes (the "?" of its group)."""
    q, cfg = res.query, res.cfg
    items = [f"Tanimoto {s.tanimoto:.2f}; ΔlogP {_signed(s.d_logp, 1)}, ΔTPSA {_signed(s.d_tpsa, 0)} Å², ΔMW {_signed(s.d_mw, 0)}.",
             f"Charge class at pH {cfg.ph_acidic:g}: {_q(q.q_low)} (query) vs {_q(s.q_low)}; at pH {cfg.ph_neutral:g}: {_q(q.q_neutral)} vs {_q(s.q_neutral)}.",
             f"In {_where(s.sources)}. SMILES: {s.smiles}"]
    return items


def _sources_help(res):
    """What was searched (counts, commit) and the licence / citation of each source (the "?" of the Sources card)."""
    items = [f"{'JP18' if s == SOURCE_JP else s}: {t}" for s, t in res.status.items()]
    if SOURCE_JP in res.searched:
        items.append(f"{jp18.LICENCE_NOTE} Page links open the MHLW PDFs.")
    if res.report_installed:
        items.append(f"RepoRT: {report.CITATION} Data {report.LICENCE} ({report.REPO_URL}, commit {res.report_commit[:10]}, {res.n_datasets} datasets); "
                     "each dataset keeps its own origin (the link of its row). Column, eluent and gradient texts are summaries of the dataset "
                     "metadata; RepoRT is LC-MS, retention times are in minutes. Reuse of RepoRT-derived values requires attribution and the same licence.")
    elif res.report_installed is False:
        items.append(f"RepoRT ({report.CITATION}; data {report.LICENCE}) is not installed: run {report.SETUP_COMMAND} (one-time download of about 100 MB).")
    return items
