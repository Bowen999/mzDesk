"""Ionization page rules: ESI polarity and response for LC-MS, EI usability for GC-MS.

Only ESI (+ / -) and EI. No APCI, APPI or CI, and no adduct lists (that is the Adducts page). All thresholds are in
``Thresholds`` (one place). The rules are coarse and qualitative: a starting point, not a measurement.

ESI+ response (first rule that fires wins)
  strong    permanent positive charge (quaternary N+ ...), or a basic site with pKa >= pH + 1 (protonated in solution)
  moderate  a basic site with pKa in [pH - 1, pH + 1) (partly protonated); or a basic site below pH - 1 (not protonated in
            solution, can still take H+ in the droplet); or no basic site but >= 3 O/N atoms (Na+ / NH4+ adducts) or a weakly
            basic group (amide, carbonyl, ether)
  weak      only C / H / halogens; or a single O/N atom in a large non-polar molecule; or O/N atoms without any weakly basic group
  Refinement A  the "no basic site" route is capped at weak when the molecule has an acidic site with pKa <= 10: an anion
                former responds far better in ESI-, so its weak [M+H]+ / adduct signal does not count.

ESI- response
  strong    permanent negative charge, or an acidic site with pKa <= 6 (carboxylic, sulfonic / sulfate, phosphate, tetrazole ...:
            these deprotonate in negative ESI even from an acidic mobile phase), or pKa <= pH - 1 (Refinement B)
  moderate  strongest acid pKa in (6, 10] (phenols, imides, sulfonamides); or >= 3 alcohol OH (polyols, sugars: formate /
            chloride adducts)
  weak      otherwise
  Refinement B  an acid that is deprotonated in solution at the chosen pH (pKa <= pH - 1) is strong whatever its pKa
                (a phenol at pH 12), the mirror image of the ESI+ rule.
  Refinement C  the polyol route is capped at weak when the molecule has a basic site (an amino sugar prefers ESI+).
  Refinement D  with a permanent positive charge, the strongest acid sites (as many as the charge) are the counterions of that
                charge (betaine, carnitine, phosphatidylcholine): there is no free acidic H, so they do not count for [M-H]-.

pKa sites: MolGpKa (the Separation page's prediction, base sites on non-basic N already dropped); without its weights the textbook
values of ``FG_BASES`` / ``FG_ACIDS`` (a note says so).
  Refinement E  amine N-oxides (R3N+-O-, pyridine N-oxide) are not matched by the MolGpKa patterns, but the O takes a proton (pKa about
                4.7 aliphatic, 0.8 aromatic): added as base sites with these textbook values (marked as not predicted).

Polarity
  both      ESI+ and ESI- both strong, or both at least moderate (with three levels the one-step condition always holds)
  poor      both weak
  + or -    otherwise, the stronger one
  For "both" a note says which polarity the chosen pH favours (ions already formed in solution).

EI
  GC-amenable  derivatization.gc_assess(), the verdict the Platform page shows (largest fragment, neutralized)
  M+•          NEIMS relative intensity at the nominal molecular mass: >= 5 % visible, 1-5 % weak, < 1 % not expected.
               NEIMS predicts the underivatized molecule. Without NEIMS weights: textbook rule of thumb instead.
"""
from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from mzdesk.platform_fit import derivatization as deriv
from mzdesk.platform_fit import lipid_metabolite_classifier as clf
from mzdesk.platform_fit import pka_engine as pe

PH_MIN, PH_MAX = 0.0, 14.0
MINUS = "−"
POS, NEG, BOTH, POOR = "+", MINUS, "both", "poor"
STRONG, MODERATE, WEAK = "strong", "moderate", "weak"
RANK = {WEAK: 0, MODERATE: 1, STRONG: 2}
VISIBLE, WEAK_ION, ABSENT = "visible", "weak", "not expected"
NEIMS_MISSING = "NEIMS not installed"
NEIMS_REPEATS = ("The structure has several components", "The input charge")
METHOD_MOLGPKA, METHOD_GROUPS = "MolGpKa", "functional groups"


@dataclass(frozen=True)
class Thresholds:
    # ESI+
    base_protonated: float = 1.0      # strong: basic site pKa >= pH + this
    base_partial: float = 1.0         # moderate: basic site pKa >= pH - this
    n_het_adducts: int = 3            # no basic site: this many O/N atoms -> moderate (adducts)
    large_heavy_atoms: int = 15       # "large" for the single-heteroatom rule (heavy atoms) ...
    nonpolar_logp: float = 3.0        # ... and "non-polar" (Crippen logP)
    # ESI-
    acid_strong_pka: float = 6.0      # strong: acidic site pKa <= this
    acid_deprotonated: float = 1.0    # strong: acidic site pKa <= pH - this
    acid_moderate_pka: float = 10.0   # moderate: strongest acid pKa <= this (also: sites above this are not shown as ionizable)
    polyol_min_oh: int = 3            # moderate without an acid: at least this many alcohol OH
    # EI
    mol_ion_visible: float = 5.0      # M+• relative intensity (%) at or above this: visible
    mol_ion_weak: float = 1.0         # ... at or above this (and below the previous): weak; below: not expected


THRESHOLDS = Thresholds()


# ---------------------------------------------------------------- pKa profile (MolGpKa, or functional groups when the weights are missing)
# Textbook pKa by group, used only when the MolGpKa weights are not installed. First pattern that matches an atom wins.
# (name, SMARTS, index of the atom that holds / takes the proton, pKa)
_NOT_AMIDE = "!$(N-[#6]=[#7,#8,#16])"
FG_BASES = (
    ("amidine / guanidine", "[NX3;!$(N-C=O);!$(N-S(=O)=O)]-[CX3;!$(C-[OX2,SX2])]=[NX2;!$(N-C=O);!$(N-[!#6;!#1])]", 2, 12.0),
    ("aliphatic amine", f"[NX3;+0;{_NOT_AMIDE};!$(N-C#N);!$(N-a);!$(N-[!#6;!#1]);!$(N-C=C)]", 0, 10.0),
    ("aniline-type N", f"[NX3;+0;{_NOT_AMIDE};!$(N-S(=O)=O);!$(N-[!#6;!#1])]-a", 0, 4.6),
    ("diazine N", "[nX2;r6;+0;$(n:c:n),$(n:n),$(n:c:c:n)]", 0, 1.5),
    ("pyridine-type N", "[nX2;r6;+0]", 0, 5.0),
    ("imidazole-type N", "[nX2;r5;+0]:c:[nX3]", 0, 6.5),
    ("azole N", "[nX2;r5;+0]", 0, 1.5),
)
FG_ACIDS = (
    ("carboxylic acid", "[CX3](=O)[OX2H1]", 2, 4.5),
    ("sulfonic / sulfate acid", "[SX4](=O)(=O)[OX2H1]", 3, 0.0),
    ("phosphate / phosphonate", "[PX4](=O)[OX2H1]", 2, 2.0),
    ("tetrazole NH", "[nH;$(n1cnnn1),$(n1nncn1)]", 0, 4.9),
    ("acyl sulfonamide", "C(=O)[NX3;H1]S(=O)=O", 2, 4.5),
    ("phenol", "[OX2H1]c", 0, 10.0),
    ("sulfonamide NH", "[NX3;H1,H2]S(=O)(=O)", 0, 9.5),
    ("imide NH", "[$([NX3H1](C=O)C=O),$([nH](c=O)c=O)]", 0, 9.5),
    ("hydroxamic acid", "C(=O)[NX3][OX2H1]", 3, 9.0),
)
N_OXIDES = (                       # (name, SMARTS, index of the O, textbook pKa of the protonated N-oxide)
    ("amine N-oxide", "[NX4+]-[OX1-]", 1, 4.7),
    ("aromatic N-oxide", "[n+]-[OX1-]", 1, 0.8),
)
_fg_compiled = []
_noxide_compiled = []


def _fg_patterns():
    if not _fg_compiled:
        for kind, table in (("base", FG_BASES), ("acid", FG_ACIDS)):
            _fg_compiled.extend((kind, name, Chem.MolFromSmarts(smarts), idx, pka) for name, smarts, idx, pka in table)
    return _fg_compiled


def group_sites(mol):
    """Fallback pKa sites from functional groups -> [pe.Site] (one per atom and type)."""
    seen, sites = set(), []
    for kind, _name, patt, idx, pka in _fg_patterns():
        for match in mol.GetSubstructMatches(patt):
            atom = match[idx]
            if (kind, atom) not in seen:
                seen.add((kind, atom))
                sites.append(pe.Site(atom=atom, element=mol.GetAtomWithIdx(atom).GetSymbol(), kind=kind, pka=pka))
    return sorted(sites + n_oxide_sites(mol), key=lambda s: (s.pka, s.atom))


def n_oxide_sites(mol):
    """Base sites on the O of amine N-oxides (Refinement E) -> [pe.Site]."""
    if not _noxide_compiled:
        _noxide_compiled.extend((Chem.MolFromSmarts(smarts), idx, pka) for _name, smarts, idx, pka in N_OXIDES)
    sites = {}
    for patt, idx, pka in _noxide_compiled:
        for match in mol.GetSubstructMatches(patt):
            sites.setdefault(match[idx], pe.Site(atom=match[idx], element="O", kind="base", pka=pka))
    return sorted(sites.values(), key=lambda s: s.atom)


@dataclass
class Profile:
    mol: Chem.Mol                 # neutral form of the largest fragment (atom order of the input when it has one fragment)
    notes: list                   # preparation notes (neutralized, salt ...)
    sites: list                   # [pe.Site], base sites on non-basic N already dropped
    ignored: list                 # [(Site, reason)]
    fixed_charge: int
    permanent: bool
    logp: float
    method: str                   # METHOD_MOLGPKA | METHOD_GROUPS
    fallback: str = ""            # why the functional-group fallback was used
    textbook: list = field(default_factory=list)      # sites whose pKa is a textbook value, not a prediction (N-oxides)


def profile(smiles):
    """SMILES -> Profile. Raises ValueError for a structure that cannot be used (unparsable, too large)."""
    try:
        a = pe.analyze(smiles)
        extra = n_oxide_sites(a.mol)
        return Profile(a.mol, list(a.notes), sorted(a.sites + extra, key=lambda s: (s.pka, s.atom)), list(a.ignored), a.fixed_charge, a.permanent,
                       a.logp, METHOD_MOLGPKA, textbook=extra)
    except pe.PkaError as exc:
        prep = pe.prepare(smiles)
        return Profile(prep.mol, list(prep.notes), group_sites(prep.mol), [], prep.fixed_charge, prep.permanent,
                       pe.crippen_logp(prep.mol), METHOD_GROUPS, fallback=str(exc))


# ---------------------------------------------------------------- structural counts for the "no basic site" route
_ACID_GROUPS = [Chem.MolFromSmarts(s) for s in ("[CX3](=[OX1])[OX2H1,OX1-]", "[SX4](=[OX1])(=[OX1])[OX2H1,OX1-]", "[PX4](=[OX1])[OX2H1,OX1-]")]
_CHARGE_SEPARATED = Chem.MolFromSmarts("[#7+]~[#8]")                 # nitro, N-oxide: not basic
_WEAK_GROUPS = (
    ("amide", Chem.MolFromSmarts("[NX3][CX3]=[OX1,SX1]")),
    ("carbonyl", Chem.MolFromSmarts("[CX3;!$(C-[OX2H1,OX1-]);!$(C-[NX3])]=[OX1]")),     # ketone, aldehyde, ester; not acid, not amide
    ("ether", Chem.MolFromSmarts("[OX2;H0]([#6])[#6]")),
)


def count_heteroatoms(mol):
    """O and N atoms that can bind H+ or Na+: neutral, not in an acid group (carboxyl, sulfonic, phosphate), not nitro / N-oxide."""
    skip = {i for patt in _ACID_GROUPS for m in mol.GetSubstructMatches(patt) for i in m}
    skip |= {i for m in mol.GetSubstructMatches(_CHARGE_SEPARATED) for i in m}
    return sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() in (7, 8) and a.GetFormalCharge() == 0 and a.GetIdx() not in skip)


def weak_base_groups(mol):
    return [name for name, patt in _WEAK_GROUPS if mol.HasSubstructMatch(patt)]


def count_alcohols(mol):
    return len(deriv.find_groups(mol).get("alcohol", {"atoms": []})["atoms"])


# ---------------------------------------------------------------- ESI+ and ESI-
@dataclass(frozen=True)
class Response:
    level: str            # strong / moderate / weak
    reason: str           # the deciding reason, one short sentence fragment (card note and table)


def _acids(prof):
    return [s for s in prof.sites if s.kind == "acid"]


def _free_acids(prof):
    """Acid sites without the counterions of a permanent positive charge (Refinement D), strongest first."""
    acids = sorted(_acids(prof), key=lambda s: (s.pka, s.atom))
    n_counter = prof.fixed_charge if prof.fixed_charge > 0 else 0
    return acids[n_counter:], acids[:n_counter]


def _plural(n, word):
    return f"{n} {word}" + ("" if n == 1 else "s")


def esi_positive(prof, ph, th=THRESHOLDS):
    bases = sorted((s for s in prof.sites if s.kind == "base"), key=lambda s: (-s.pka, s.atom))
    if prof.fixed_charge > 0:
        return Response(STRONG, f"permanent positive charge ({prof.fixed_charge:+d}): always an ion")
    if bases:
        b = bases[0]
        if b.pka >= ph + th.base_protonated:
            return Response(STRONG, f"basic site {b.label} (pKa {b.pka:.1f}) is protonated in solution at pH {ph:g} (pKa ≥ pH + {th.base_protonated:g})")
        if b.pka >= ph - th.base_partial:
            return Response(MODERATE, f"basic site {b.label} (pKa {b.pka:.1f}) is within ± {th.base_partial:g} of pH {ph:g}: partly protonated")
        return Response(MODERATE, f"basic site {b.label} (pKa {b.pka:.1f}) is not protonated in solution at pH {ph:g}, but can still take H+ in the droplet")
    # no basic site: heteroatoms and weakly basic groups
    n_het, groups = count_heteroatoms(prof.mol), weak_base_groups(prof.mol)
    heavy = prof.mol.GetNumHeavyAtoms()
    large_nonpolar = heavy >= th.large_heavy_atoms and prof.logp >= th.nonpolar_logp
    if n_het == 0:
        return Response(WEAK, "no basic site and no O/N atom that can bind H+ or Na+ (only C, H, halogens, or nitro / N-oxide / acid groups)")
    if n_het == 1 and large_nonpolar:
        return Response(WEAK, f"no basic site; a single O/N atom in a large non-polar molecule ({heavy} heavy atoms, logP {prof.logp:.1f})")
    if n_het >= th.n_het_adducts or groups:
        why = f"{n_het} O/N atoms" + (f" and weakly basic {' / '.join(groups)}" if groups else "")
        level, reason = MODERATE, f"no basic site; {why}: weak [M+H]+, mostly Na+ / NH4+ adducts"
        free, _ = _free_acids(prof)
        acid = next((s for s in free if s.pka <= th.acid_moderate_pka), None)
        if acid:                                                   # Refinement A
            return Response(WEAK, f"no basic site; {why} would give adducts, but the acidic site {acid.label} (pKa {acid.pka:.1f}) makes it an anion former: ESI− is the better mode")
        return Response(level, reason)
    return Response(WEAK, f"no basic site; {_plural(n_het, 'O/N atom')} and no carbonyl, amide or ether group")


def esi_negative(prof, ph, th=THRESHOLDS):
    if prof.fixed_charge < 0:
        return Response(STRONG, f"permanent negative charge ({prof.fixed_charge:+d}): always an ion")
    free, counter = _free_acids(prof)
    if free:
        a = free[0]
        if a.pka <= th.acid_strong_pka:
            return Response(STRONG, f"acidic site {a.label} (pKa {a.pka:.1f} ≤ {th.acid_strong_pka:g}): deprotonates in negative ESI even from an acidic mobile phase")
        if a.pka <= ph - th.acid_deprotonated:                    # Refinement B
            return Response(STRONG, f"acidic site {a.label} (pKa {a.pka:.1f}) is deprotonated in solution at pH {ph:g}")
        if a.pka <= th.acid_moderate_pka:
            return Response(MODERATE, f"weak acid {a.label} (pKa {a.pka:.1f}, between {th.acid_strong_pka:g} and {th.acid_moderate_pka:g}): deprotonates only partly")
    n_oh = count_alcohols(prof.mol)
    has_base = prof.fixed_charge > 0 or any(s.kind == "base" for s in prof.sites)
    if n_oh >= th.polyol_min_oh:
        if has_base:                                               # Refinement C
            return Response(WEAK, f"no acidic site; {n_oh} alcohol OH would give formate / chloride adducts, but the basic site favours ESI+")
        return Response(MODERATE, f"no acidic site below pKa {th.acid_moderate_pka:g}; {n_oh} alcohol OH (polyol): formate / chloride adducts")
    if counter:
        c = counter[0]
        return Response(WEAK, f"the acid {c.label} (pKa {c.pka:.1f}) is the counterion of the permanent positive charge: no free acidic H, only anion adducts")
    if free:
        return Response(WEAK, f"strongest acidic site {free[0].label} has pKa {free[0].pka:.1f} (above {th.acid_moderate_pka:g}): not deprotonated in ESI−")
    return Response(WEAK, "no acidic site")


def solution_favour(prof, ph, th=THRESHOLDS):
    """Which charge is already formed in solution at ``ph``: POS, NEG or None (both or neither)."""
    n_pos = sum(s.kind == "base" and s.pka >= ph + th.base_protonated for s in prof.sites)
    free, _ = _free_acids(prof)
    n_neg = sum(s.pka <= ph - th.acid_deprotonated for s in free)
    return POS if n_pos > n_neg else NEG if n_neg > n_pos else None


def polarity(pos, neg, prof=None, ph=None, th=THRESHOLDS):
    """(+ / − / both / poor, note)."""
    rp, rn = RANK[pos.level], RANK[neg.level]
    summary = f"ESI+ {pos.level}, ESI− {neg.level}"
    if rp == 0 and rn == 0:
        return POOR, summary
    if rp >= 1 and rn >= 1:
        if prof is None:
            return BOTH, summary
        fav = solution_favour(prof, ph, th)
        if fav:
            charge = "positive" if fav == POS else "negative"
            return BOTH, f"{summary}; at pH {ph:g} {charge} ions are already formed in solution, so ESI{fav} is usually the more sensitive mode"
        return BOTH, f"{summary}; try both modes at pH {ph:g}"
    return (POS if rp > rn else NEG), summary


# ---------------------------------------------------------------- EI: GC amenability and the molecular ion
@dataclass
class MolIon:
    value: str                    # card value: visible / weak / not expected / "NEIMS not installed" / "not available"
    level: str = ""               # VISIBLE / WEAK_ION / ABSENT, "" when unknown
    rel: float = None             # NEIMS relative intensity (%) at the nominal mass
    mz: int = None
    note: str = ""
    warnings: list = field(default_factory=list)


_AROMATIC = Chem.MolFromSmarts("a")
_CONJUGATED = Chem.MolFromSmarts("[CX3]=[CX3]-[CX3]=[CX3,OX1,NX2]")
_ALIPHATIC_ALCOHOL = Chem.MolFromSmarts("[OX2H1][CX4]")
_BRANCHED = Chem.MolFromSmarts("[CX4;H0;!R;$(C(-[#6])(-[#6])(-[#6])-[#6])]")        # acyclic quaternary carbon


def ion_level(rel, th=THRESHOLDS):
    return VISIBLE if rel >= th.mol_ion_visible else WEAK_ION if rel >= th.mol_ion_weak else ABSENT


def rule_of_thumb(mol):
    """Textbook M+• trend at 70 eV (used only when NEIMS is not installed) -> (level, reason)."""
    if mol.HasSubstructMatch(_AROMATIC) or mol.HasSubstructMatch(_CONJUGATED):
        return VISIBLE, "aromatic or conjugated systems stabilize M+•"
    if mol.HasSubstructMatch(_ALIPHATIC_ALCOHOL):
        return ABSENT, "aliphatic alcohols lose H2O or an alkyl group so fast that M+• is usually absent"
    if mol.HasSubstructMatch(_BRANCHED):
        return ABSENT, "highly branched alkyl: M+• cleaves at the branch points and is usually absent"
    if mol.GetRingInfo().NumRings():
        return VISIBLE, "ring systems need two cleavages to fall apart, so M+• usually survives"
    return WEAK_ION, "acyclic saturated molecules: M+• is weak and fades with chain length"


def molecular_ion(smiles, mol, th=THRESHOLDS):
    """M+• of the underivatized molecule: NEIMS relative intensity at the nominal mass, or the rule of thumb without NEIMS weights."""
    from mzdesk.ms2 import neims_engine as ne

    try:
        res = ne.predict(smiles)
    except ne.NEIMSError:
        level, why = rule_of_thumb(mol)
        return MolIon(NEIMS_MISSING, level, note=f"rule of thumb: M+• {level} ({why})")
    except ValueError as exc:
        level, why = rule_of_thumb(mol)
        return MolIon("not available", level, note=f"NEIMS cannot predict this structure ({exc}); rule of thumb: M+• {level} ({why})")
    mz = res["molecular_ion_mz"]
    rel = dict(res["peaks"]).get(mz, 0.0)                        # peaks below NEIMS's display threshold count as 0
    level = ion_level(rel, th)
    warnings = [w for w in res["warnings"] if not w.startswith(NEIMS_REPEATS)]       # salt / neutralization notes: already given by the pKa step
    return MolIon(level, level, rel, mz, f"NEIMS: {rel:.1f} % of the base peak at m/z {mz} (nominal)", warnings)


def gc_verdict(smiles):
    """The Platform page's GC-MS verdict for the neutral form of the largest fragment (salts and mixtures: GC and EI see the main component)
    -> (gc_assess dict, neutral molecule)."""
    mol = Chem.MolFromSmiles(smiles)
    if len(Chem.GetMolFrags(mol)) > 1:
        mol = rdMolStandardize.LargestFragmentChooser(preferOrganic=True).choose(mol)
    mol, _ = deriv.neutralize(mol)
    groups = deriv.find_groups(mol)
    return deriv.gc_assess(mol, groups, clf.detect_glycosyl(mol)[0]), mol


EI_VERDICT = {"yes": "usable", "after derivatization": "usable after derivatization", "no": "not usable"}


def ion_caveat(gc):
    """When the M+• does not describe the molecule that reaches the detector ("" otherwise)."""
    if gc["gc"] == "after derivatization":
        return "derivatization is needed and NEIMS predicts the underivatized molecule, so the derivative has another M+•"
    if gc["gc"] == "no":
        return "hypothetical: the molecule is not GC-amenable"
    return ""


def ion_card_note(gc, ion):
    """Where the M+• value comes from, and when it does not describe the molecule that reaches the detector."""
    caveat = ion_caveat(gc)
    return f"{ion.note}; {caveat}" if caveat else ion.note


def ei_reason(gc, ion):
    """Reason text of the EI row: the GC rule, then the molecular ion."""
    parts = [r.removeprefix("Rule: ") for r in gc["reasons"]]
    if ion.rel is not None:
        parts.append(f"M+• (m/z {ion.mz}) {ion.value}: {ion.rel:.1f} % of the base peak by NEIMS")
    else:
        parts.append(f"{ion.value}; {ion.note}" if ion.value == NEIMS_MISSING else ion.note)
    if gc["gc"] == "after derivatization":
        parts.append("the M+• is for the underivatized molecule")
    return "; ".join(parts)


# ---------------------------------------------------------------- everything for one structure
@dataclass
class Assessment:
    ph: float
    profile: Profile
    pos: Response
    neg: Response
    polarity: str
    polarity_note: str
    gc: dict
    ion: MolIon
    warnings: list
    hits: list


def assess(smiles, ph=2.7, th=THRESHOLDS):
    """SMILES -> Assessment. Raises ValueError for an unusable structure or a pH outside 0-14."""
    if not PH_MIN <= ph <= PH_MAX:
        raise ValueError(f"pH must be between {PH_MIN:g} and {PH_MAX:g}.")
    prof = profile(smiles)
    pos, neg = esi_positive(prof, ph, th), esi_negative(prof, ph, th)
    pol, pol_note = polarity(pos, neg, prof, ph, th)
    gc, gc_mol = gc_verdict(smiles)
    ion = molecular_ion(smiles, gc_mol, th)
    warnings = [n for n in prof.notes if not n.startswith("Permanent charge")]
    warnings += ion.warnings
    return Assessment(ph, prof, pos, neg, pol, pol_note, gc, ion, warnings, rule_hits(prof, ph, pos, neg, gc, ion, th))


def rule_rows(a):
    """[[mode, verdict, reason], ...] for ESI+, ESI−, EI."""
    return [["ESI+", a.pos.level, a.pos.reason], [f"ESI{MINUS}", a.neg.level, a.neg.reason],
            ["EI", EI_VERDICT[a.gc["gc"]], ei_reason(a.gc, a.ion)]]


def rule_hits(prof, ph, pos, neg, gc, ion, th=THRESHOLDS):
    hits = []
    if prof.method == METHOD_MOLGPKA:
        hits.append("pKa sites from MolGpKa (the same prediction as the Separation page); the pKa of a site carries an error of about one unit.")
    else:
        hits.append(f"MolGpKa weights not available ({prof.fallback}): pKa values below are textbook values of the functional groups, "
                    "not a prediction; the rules still run.")
    sites = sorted(prof.sites, key=lambda s: (s.pka, s.atom))
    shown = [s for s in sites if s.kind == "base" or s.pka <= th.acid_moderate_pka]
    if shown:
        hits.append("Ionizable sites used: " + ", ".join(f"{s.label} {s.kind} {s.pka:.1f}" for s in shown) + ".")
    if prof.textbook:
        hits.append("N-oxide site" + ("s" if len(prof.textbook) > 1 else "") + " " + ", ".join(s.label for s in prof.textbook)
                    + ": the MolGpKa patterns do not cover N-oxides, so the pKa is a textbook value (not predicted).")
    if len(shown) < len(sites):
        n = len(sites) - len(shown)
        hits.append(f"{_plural(n, 'acid site')} with pKa above {th.acid_moderate_pka:g} (alcohol OH, amide NH ...) {'is' if n == 1 else 'are'} too weak to deprotonate in ESI− and "
                    f"{'is' if n == 1 else 'are'} not listed.")
    if not sites and not prof.permanent:
        hits.append("No ionizable site found.")
    hits += [f"{s.label} {why}: not basic, base site ignored (MolGpKa predicted pKa {s.pka:.1f})." for s, why in prof.ignored]
    if prof.fixed_charge > 0:
        hits.append(f"Permanent positive charge ({prof.fixed_charge:+d}; quaternary N+, N-alkylpyridinium ...): an ion at any pH, so ESI+ is strong; "
                    "an acid in the same molecule is its counterion.")
    hits.append(f"ESI+ rests on the solution pKa of the basic site against pH {ph:g} (strong at pKa ≥ {ph + th.base_protonated:.1f}); ESI− on the acid pKa "
                f"(strong at ≤ {th.acid_strong_pka:g}, moderate to {th.acid_moderate_pka:g}). Solution pKa is not gas-phase behaviour: gas-phase proton affinity "
                "and the droplet chemistry can move a compound up or down, mostly for weak bases (pyridine-type N, carbonyls) and for acids at low pH.")
    hits.append("Response strength is qualitative (strong / moderate / weak, relative to other small molecules); it says nothing about absolute sensitivity, "
                "matrix effects or the ion source settings.")
    if gc["gc"] == "after derivatization":
        hits.append("Derivatization is needed for GC: NEIMS predicts the underivatized molecule, so the M+• and the spectrum shown are not those of the "
                    "TMS or MeOx-TMS derivative that is actually measured (masses: Derivatization page).")
    if ion.rel is not None:
        hits.append(f"M+• rule: NEIMS relative intensity at the nominal mass ≥ {th.mol_ion_visible:g} % visible, {th.mol_ion_weak:g}-{th.mol_ion_visible:g} % weak, "
                    f"< {th.mol_ion_weak:g} % not expected. A neural-network prediction: the intensity of a single peak is far less certain than the peak pattern.")
    else:
        hits.append("NEIMS weights not installed (python -m mzdesk.ms2.neims_setup), or the structure is outside its range: the M+• line is a textbook rule of thumb.")
    return hits
