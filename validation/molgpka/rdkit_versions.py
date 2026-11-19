"""Does the RDKit version change the numpy engine's predictions?  Run with any Python that has RDKit and numpy:

    <python> -I validation/molgpka/rdkit_versions.py

Compares the engine, on the RDKit of that Python, with the per-site pKa the ORIGINAL code produced (original_output.json,
same_graph level, made under the RDKit that run_original.sh used). Molecules this RDKit cannot parse are skipped and counted.
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))
from rdkit import Chem, __version__          # noqa: E402

from mzdesk.platform_fit import pka_engine as E          # noqa: E402


def main():
    data = json.loads((HERE / "original_output.json").read_text())
    n = equal = sites = skipped = 0
    worst, big = 0.0, 0
    for rec in data["molecules"]:
        mol = Chem.MolFromSmiles(rec["smiles"])
        if mol is None:
            skipped += 1
            continue
        mol_h = Chem.AddHs(mol)
        acid, base = E.find_site_atoms(mol_h)
        n += 1
        equal += (sorted(int(a) for a in rec["same_graph"]["acid"]) == acid and sorted(int(a) for a in rec["same_graph"]["base"]) == base)
        for kind, aids in (("acid", acid), ("base", base)):
            got = E.predict_aids(mol_h, kind, aids)
            for a in aids:
                ref = rec["same_graph"][kind].get(str(a))
                if ref is not None:
                    sites += 1
                    worst = max(worst, abs(got[a] - ref))
                    big += abs(got[a] - ref) > 0.05
    print(f"RDKit {__version__}: {n} molecules ({skipped} unparsable), site sets equal {equal}, {sites} sites, max |diff| {worst:.2g} pKa, {big} sites differ by > 0.05")


if __name__ == "__main__":
    main()
