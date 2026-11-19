"""Step 6a: merge the layers' candidates into one list with one Candidate per structure.

Candidates are the same structure when their InChIKey first blocks (ik14: formula + connectivity + H layers, no stereo, no
isotopes) are equal, so stereoisomers merge (a spectrum does not tell them apart) and constitutional isomers do not. The key is
computed from the SMILES (largest fragment); a candidate without a valid SMILES falls back to the InChIKey the layer gave.

A merged candidate keeps every contributing layer in `sources` (the best contribution per layer). Its own fields come from the
layer with the strongest evidence (library > MS-FINDER > FRIGID); fields that layer lacks are filled from the next one. A candidate
that has no valid SMILES after that is kept, flagged FLAG_NO_STRUCTURE, and is never scored by detectability.
"""
from __future__ import annotations

import copy

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.rdBase import BlockLogs

from mzdesk.annotation.levels import assign_level
from mzdesk.annotation.schema import FLAG_NO_STRUCTURE, LAYER_FORMULA, LAYERS, Candidate, evidence_of


# ---------------------------------------------------------------- structures
def parse_smiles(smiles):
    """RDKit Mol of a SMILES, or None when it is empty, unparsable or has no heavy atom."""
    if not smiles or not str(smiles).strip():
        return None
    with BlockLogs():
        mol = Chem.MolFromSmiles(str(smiles).strip())
    return mol if mol is not None and mol.GetNumHeavyAtoms() else None


def has_structure(smiles):
    """True when `smiles` is a usable structure (detectability can only be scored for these)."""
    return parse_smiles(smiles) is not None


def structure_keys(smiles):
    """(InChIKey, ik14) of the largest fragment of a SMILES; (None, None) when it is not a usable structure."""
    mol = parse_smiles(smiles)
    if mol is None:
        return None, None
    parts = Chem.GetMolFrags(mol, asMols=True)
    if len(parts) > 1:                                  # salt or mixture: the largest part, as the other pages do
        mol = max(parts, key=lambda p: (p.GetNumHeavyAtoms(), Descriptors.ExactMolWt(p)))
    with BlockLogs():
        key = Chem.MolToInchiKey(mol)
    return (key, key[:14]) if key else (None, None)


def _key14(candidate):
    """ik14 a candidate is merged by: from its SMILES, else the InChIKey the layer gave; None when it has neither."""
    ik14 = structure_keys(candidate.smiles)[1]
    if ik14:
        return ik14
    given = (candidate.ik14 or (candidate.inchikey or "")[:14]).strip().upper()
    return given if len(given) == 14 else None


# ---------------------------------------------------------------- merge
def merge_candidates(layer_results):
    """LayerResults (in layer order) -> unique Candidates, one per ik14, with id, sources, flags and level filled.

    Ids are "c1", "c2", ... in order of first appearance. The input candidates are not changed. A candidate without any key is
    never merged. Several candidates of one layer with the same ik14 count once (the best layer_score_norm).
    """
    groups = {}                                         # ik14 (or a unique tuple) -> [(rank_in_layer, candidate)], first seen first
    for result in layer_results:
        for rank, cand in enumerate(result.candidates, start=1):
            if cand.layer not in LAYERS:
                raise ValueError(f"A layer result can only hold candidates of {', '.join(LAYERS)}, not {cand.layer!r}.")
            groups.setdefault(_key14(cand) or ("unkeyed", len(groups)), []).append((rank, cand))
    return [_merge_group(f"c{n}", members) for n, members in enumerate(groups.values(), start=1)]


def _merge_group(new_id, members):
    best = {}                                           # layer -> (rank_in_layer, candidate): the best contribution of each layer
    for rank, cand in members:
        if cand.layer not in best or cand.layer_score_norm > best[cand.layer][1].layer_score_norm:
            best[cand.layer] = (rank, cand)
    ranked = [best[layer] for layer in LAYERS if layer in best]          # strongest evidence first
    contributors = [cand for _, cand in ranked]
    merged = copy.deepcopy(contributors[0])

    for name in ("name", "formula", "adduct", "mz_error_ppm"):
        if getattr(merged, name) in (None, ""):
            setattr(merged, name, next((getattr(c, name) for c in contributors if getattr(c, name) not in (None, "")), getattr(merged, name)))
    smiles = next((c.smiles for c in contributors if has_structure(c.smiles)), None)
    if smiles:
        merged.smiles = smiles
        merged.inchikey, merged.ik14 = structure_keys(smiles)
    else:                                               # no usable structure: keep what the layers gave
        merged.smiles = next((c.smiles for c in contributors if c.smiles), None)
        merged.inchikey = next((c.inchikey for c in contributors if c.inchikey), None)
        merged.ik14 = _key14(merged)

    merged.id = new_id
    merged.evidence = {c.layer: copy.deepcopy(evidence_of(c, c.layer)) for c in contributors}
    merged.hypothesis = all(c.hypothesis for c in contributors)
    merged.sources = [{"layer": c.layer, "layer_score": c.layer_score, "layer_score_norm": c.layer_score_norm, "rank_in_layer": rank}
                      for rank, c in ranked]
    merged.flags = [f for f in dict.fromkeys(f for c in contributors for f in c.flags) if f != FLAG_NO_STRUCTURE]
    if not smiles:
        merged.flags.append(FLAG_NO_STRUCTURE)
    merged.detect = merged.base_score = merged.final_score = merged.rank = merged.isomer = None
    assign_level(merged)
    return merged


# ---------------------------------------------------------------- formula-only rows
def formula_only_candidates(formulas):
    """Rows for the case that no structure was proposed: one Candidate per molecular formula, Level 4, ranked in the given order.

    `formulas`: the shared formula list of step 2, as formula strings or dicts with "formula" and optionally "adduct", "prob"
    (BUDDY probability: becomes layer_score and layer_score_norm), "ppm" (mz_error_ppm) and "isotope_sim". These rows are not
    scored by detectability: base_score and final_score stay None, rank is the position in the list.
    """
    rows = []
    for rank, item in enumerate(formulas or [], start=1):
        item = {"formula": item} if isinstance(item, str) else item
        prob = item.get("prob")
        row = Candidate(id=f"f{rank}", layer=LAYER_FORMULA, formula=item["formula"], adduct=item.get("adduct"),
                        mz_error_ppm=item.get("ppm"), layer_score=prob, layer_score_norm=float(prob or 0.0),
                        evidence={LAYER_FORMULA: {k: item[k] for k in ("prob", "ppm", "isotope_sim") if item.get(k) is not None}},
                        flags=[FLAG_NO_STRUCTURE], rank=rank)
        assign_level(row)
        rows.append(row)
    return rows
