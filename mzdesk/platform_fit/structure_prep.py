"""结构预处理的共用小函数：取最大片段、质子化不挨着阳离子的阴离子。pka_engine（pKa / 电离 / 加合物 / 分离页）和 derivatization（Platform / Derivatization 页）共用。

只依赖 RDKit，不带副作用（不改 RDKit 日志），所以 derivatization 可以直接导入。
"""
from rdkit import Chem
from rdkit.Chem import Descriptors


def neutralize_anions(mol):
    """Protonate O- / S- / N- that are not next to a cation: carboxylate, phosphate, phenolate, sulfonate ... -> (mol, changed).

    RDKit's Uncharger leaves such an anion alone when a quaternary N+ is in the molecule (betaine, phosphatidylcholine), which
    would hide an ionizable site. Charge-separated groups (nitro, N-oxide, azide: anion next to a cation) are left as they are.
    """
    rw = Chem.RWMol(mol)
    changed = False
    for atom in rw.GetAtoms():
        if atom.GetFormalCharge() == -1 and atom.GetSymbol() in ("O", "S", "N") and not any(n.GetFormalCharge() > 0 for n in atom.GetNeighbors()):
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
            changed = True
    if not changed:
        return mol, False
    try:
        Chem.SanitizeMol(rw)
    except Exception:  # noqa: BLE001 — keep the molecule as it was
        return mol, False
    return rw.GetMol(), True


def largest_fragment(mol):
    """-> (the main component of a salt or mixture, number of components): the one with carbon, then the most heavy atoms, then the highest MW."""
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) == 1:
        return mol, 1
    def key(f):
        return (any(a.GetSymbol() == "C" for a in f.GetAtoms()), f.GetNumHeavyAtoms(), Descriptors.MolWt(f))
    return max(frags, key=key), len(frags)


def parts_note(n_parts):
    """The warning every page gives when it analyses only the largest of ``n_parts`` (> 1) disconnected parts."""
    return f"{n_parts} disconnected parts (salt or mixture): only the largest is analysed."
