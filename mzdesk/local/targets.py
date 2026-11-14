"""Which page or workflow needs which component: needs("workflow:esi_ms2_annotation", form values) -> what is missing, blocked, unsupported.

Decided from the code (runners in mzdesk/api/runners, workflows in mzdesk/annotation and mzdesk/workflows), not guessed:
  module:ms2_library      libraries (spectrum type EI or MS/MS, the ticked libraries)
  module:ms2_structure    msfinder           module:ms2_denovo  frigid        module:ms2_esi  cfmid (Docker)      module:ms2_ei  neims
  module:ccs              ccs                module:known_methods  repo + molgpka (the RepoRT index is classified with MolGpKa)
  module:separation / ionization / adducts   molgpka (pka_engine).  Not platform_fit: its runner uses the lipid / derivatization rules only.
  workflow:esi_ms2_annotation   libraries + msfinder + molgpka without asking, ccs when an observed CCS is given, frigid (asks) when Run FRIGID is on
  workflow:ei_ms_annotation     libraries + neims (both ask)
  workflow:targeted_method      molgpka + repo
  module:gnps_redu        redu               module:gnps_fastmasst  redu when Join Pan-ReDU is ticked
Everything else (the ionization page's EI half, Targeted Method's GC route) has a rule-of-thumb fallback when NEIMS is missing and does not ask for it.
"""
from mzdesk.local import components, prereq

MODULES = {
    "ms2_library": ("libraries",),
    "ms2_structure": ("msfinder",),
    "ms2_denovo": ("frigid",),
    "ms2_esi": ("cfmid",),
    "ms2_ei": ("neims",),
    "ccs": ("ccs",),
    "known_methods": ("repo", "molgpka"),
    "separation": ("molgpka",),
    "ionization": ("molgpka",),
    "adducts": ("molgpka",),
    "gnps_redu": ("redu",),
}
WHY = {
    "libraries": "Library search needs the MassBank and GNPS spectral libraries.",
    "msfinder": "MS-FINDER predicts the structure candidates.",
    "frigid": "FRIGID generates de novo structures.",
    "cfmid": "CFM-ID predicts the ESI MS/MS spectrum.",
    "neims": "NEIMS predicts the EI spectrum.",
    "molgpka": "MolGpKa predicts the pKa values.",
    "ccs": "The observed CCS is compared with the CCS reference data.",
    "repo": "RepoRT lists the published chromatography methods.",
    "redu": "Dataset Finder searches the Pan-ReDU sample metadata of GNPS2.",
}
DOCKER_MISSING = "CFM-ID runs in Docker. Install Docker Desktop, start it, and run again."
DOCKER_DOWN = "Docker is not running. Start Docker Desktop and run again."


class UnknownTarget(ValueError):
    pass


def _truthy(value):
    return value is True or str(value).strip().lower() in ("true", "1", "yes", "on")


def _filled(value):
    return value is not None and str(value).strip() != ""


def _libraries_missing(selected):
    """True when one of the ticked libraries has no index (``selected`` None = both)."""
    names = components.LIBRARY_NAMES if selected is None else [n for n in selected if n in components.LIBRARY_NAMES]
    return any(not components.library_installed(n) for n in names)


def requirements(target, inputs):
    """-> [(component id, ask override or None, why or None)] for the target. UnknownTarget when the target is not a page of the app."""
    from mzdesk.spec import MODULES as SPEC_MODULES, PIPELINES

    inputs = inputs if isinstance(inputs, dict) else {}
    kind, _, name = str(target).partition(":")
    if kind == "module" and name in {m["id"] for m in SPEC_MODULES}:
        if name == "ms2_library":
            chosen = inputs.get("libraries")
            return [("libraries", None, None)] if not isinstance(chosen, list) or _libraries_missing(chosen) else []
        if name == "gnps_fastmasst":
            return [("redu", None, "Join Pan-ReDU is ticked.")] if _truthy(inputs.get("redu_join")) else []
        return [(cid, None, None) for cid in MODULES.get(name, ())]
    if kind == "workflow" and name in {p["id"] for p in PIPELINES}:
        if name == "esi_ms2_annotation":
            chosen = inputs.get("libraries")
            chosen = chosen if isinstance(chosen, list) else []
            out = [("libraries", False, None)] if chosen and _libraries_missing(chosen) else []
            out += [("msfinder", False, None), ("molgpka", False, None)]
            if _filled(inputs.get("ccs_observed")):
                out.append(("ccs", False, "The observed CCS is compared with the CCS reference data."))
            if _truthy(inputs.get("run_frigid")):
                out.append(("frigid", True, "Run FRIGID is ticked."))
            return out
        if name == "ei_ms_annotation":
            return [("libraries", True, None), ("neims", True, None)]
        if name == "targeted_method":
            return [("molgpka", None, None), ("repo", None, None)]
        return []
    raise UnknownTarget(f"Unknown target {target!r}.")


def needs(target, inputs=None):
    """-> {"missing": [{id, title, size, ask, why}], "blocked": [{id, title, message, link}], "unsupported": [{id, title, message}]}.

    ask = the user must confirm the download (the component's own flag, except where the workflow says otherwise: MS/MS Annotation installs the
    libraries and MS-FINDER without asking). A component this platform cannot have is "unsupported" and never blocks; one that cannot be fixed
    by downloading (Docker missing or not running) is "blocked". RepoRT brings MolGpKa with it when that is missing too."""
    missing, blocked, unsupported = [], [], []
    seen = set()
    wanted = list(requirements(target, inputs))
    if any(cid == "repo" for cid, _, _ in wanted) and not any(cid == "molgpka" for cid, _, _ in wanted):
        wanted.append(("molgpka", None, None))
    for cid, ask, why in wanted:
        if cid in seen:
            continue
        seen.add(cid)
        comp = components.REGISTRY[cid]
        if not comp.supported():
            unsupported.append({"id": cid, "title": comp.title, "message": comp.unsupported_message()})
            continue
        if cid == "cfmid":
            from mzdesk.ms2 import cfmid_engine as cf
            state = prereq.docker_state(cf.IMAGE)
            if state == "no-docker":
                blocked.append({"id": cid, "title": comp.title, "message": DOCKER_MISSING, "link": prereq.DOCKER_LINK})
                continue
            if state == "down":
                blocked.append({"id": cid, "title": comp.title, "message": DOCKER_DOWN, "link": None})
                continue
            if state == "ready":
                continue
        elif comp.installed():
            continue
        missing.append({"id": cid, "title": comp.title, "size": comp.size, "ask": comp.ask if ask is None else ask, "why": why or WHY[cid]})
    return {"missing": missing, "blocked": blocked, "unsupported": unsupported}
