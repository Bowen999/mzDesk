#!/usr/bin/env python3
"""Package the corrected CASMI 2022 key and extracted spectra as JSON Lines.

The generated ``inputs`` object is a blind mzDesk ESI MS/MS Annotation form;
the official answer is kept separately in ``ground_truth``.  Requires the 500
locally extracted SIRIUS .ms files, which are not tracked in git.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


BASE = Path(__file__).resolve().parent
OUTPUT = BASE / "casmi2022_ground_truth.jsonl"
CANONICAL_ADDUCTS = {
    "[M+H]+": "[M+H]+", "[M-H]-": "[M-H]-", "[M+NH4]+": "[M+NH4]+",
    "[M+Na]+": "[M+Na]+", "[M+Cl]-": "[M+Cl]-", "[M+FA-H]-": "[M+HCOO]-",
}


def table(path, id_column, delimiter):
    with path.open(encoding="utf-8", newline="") as handle:
        result = {}
        for row in csv.DictReader(handle, delimiter=delimiter):
            key = int(row[id_column])
            if key in result:
                raise ValueError(f"duplicate challenge {key} in {path}")
            result[key] = row
    return result


def parse_ms(path):
    """Read each peak block without combining collision energies."""
    meta, blocks, block = {}, [], None
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            key, _, value = line[1:].partition(" ")
            if key in ("ms1", "ms2", "collision"):
                block = {"kind": key, "collision_energy": float(value) if key == "collision" else None,
                         "peaks": []}
                blocks.append(block)
            else:
                if key in meta or block is not None:
                    raise ValueError(f"unexpected header in {path}: {line}")
                meta[key] = value
        else:
            if block is None:
                raise ValueError(f"peak outside spectrum block: {path}")
            parts = line.split()
            if len(parts) != 2:
                raise ValueError(f"invalid peak in {path}: {line}")
            mz, intensity = map(float, parts)
            if not all(math.isfinite(x) and x > 0 for x in (mz, intensity)):
                raise ValueError(f"non-positive or non-finite peak in {path}: {line}")
            block["peaks"].append([mz, intensity])
    return meta, blocks


def merge_peaks(blocks, ppm=10.0):
    """Match benchmark/tools/jm/prep.py: fixed-start 10 ppm clusters, max intensity."""
    peaks = sorted((p for block in blocks for p in block["peaks"]), key=lambda p: p[0])
    merged, start = [], 0
    for i in range(1, len(peaks) + 1):
        if i == len(peaks) or peaks[i][0] - peaks[start][0] > peaks[start][0] * ppm * 1e-6:
            merged.append(max(peaks[start:i], key=lambda p: p[1]))
            start = i
    return merged


def peak_text(peaks):
    return "\n".join(f"{mz:.6f} {intensity:.2f}" for mz, intensity in peaks)


def check_equal(left, right, fields, cid):
    for field in fields:
        if left[field].strip() != right[field].strip():
            raise ValueError(f"challenge {cid}: {field} differs between metadata and key")


def build_records():
    meta = table(BASE / "metadata/casmi2022_challenges.csv", "challenge_id", ",")
    truth = table(BASE / "sirius_run/results/truth.tsv", "challenge_id", "\t")
    report = table(BASE / "sirius_run/ms/extraction_report.tsv", "id", "\t")
    expected_ids = set(range(1, 501))
    if any(set(t) != expected_ids for t in (meta, truth, report)):
        raise ValueError("metadata, truth and extraction report must each contain IDs 1–500")

    records = []
    for cid in sorted(expected_ids):
        m, t, r = meta[cid], truth[cid], report[cid]
        check_equal(m, t, ("priority", "polarity", "adduct", "formula", "inchikey", "smiles"), cid)
        check_equal(m, r, ("priority", "polarity", "file"), cid)
        for field in ("precursor_mz", "monoisotopic_mass"):
            if abs(float(m[field]) - float(t[field])) > 1e-7:
                raise ValueError(f"challenge {cid}: inconsistent {field}")
        for source, reported in (("rt_min", "rt_given"), ("precursor_mz", "precursor_given")):
            if abs(float(m[source]) - float(r[reported])) > 1e-7:
                raise ValueError(f"challenge {cid}: inconsistent {source}")

        path = BASE / "sirius_run/ms" / f"challenge_{cid:03d}.ms"
        header, blocks = parse_ms(path)
        ms1 = [b for b in blocks if b["kind"] == "ms1"]
        ms2 = [b for b in blocks if b["kind"] in ("ms2", "collision")]
        if (header.get("compound") != f"challenge_{cid:03d}"
                or header.get("ionization") != ("[M+?]+" if m["polarity"] == "positive" else "[M+?]-")
                or abs(float(header["parentmass"]) - float(r["parentmass_used"])) > 1e-7
                or len(ms1) != 1 or len(ms1[0]["peaks"]) != int(r["n_ms1_isotope_peaks"])
                or len(ms2) != int(r["n_ms2_scans_used"])
                or sum(len(b["peaks"]) for b in ms2) != int(r["n_ms2_peaks"])
                or not ms2 or not all(b["peaks"] for b in ms2)):
            raise ValueError(f"challenge {cid}: spectrum and extraction report disagree")

        merged = merge_peaks(ms2)
        if not merged or len(merged) > 5000:
            raise ValueError(f"challenge {cid}: invalid merged MS2 peak count")
        adduct = m["adduct"].strip()
        if adduct not in CANONICAL_ADDUCTS:
            raise ValueError(f"challenge {cid}: unknown adduct {adduct}")

        inputs = {
            "precursor": {"precursor_mz": float(header["parentmass"]),
                          "ion_mode": m["polarity"].capitalize(), "adduct": "Auto"},
            "ms2": {"source": "Paste", "value": peak_text(merged)},
            "ms1": peak_text(ms1[0]["peaks"]),
            "libraries": ["MassBank", "GNPS"], "tol_ppm": 10.0, "tol_da": 0.02,
            "omics": "Either / unknown", "halogen": False, "local_isomers": False,
            "additive": "Unknown", "alkali": "Unknown", "ph": 2.7,
        }
        records.append({
            "challenge_id": cid,
            "inputs": inputs,
            "ground_truth": {
                "formula": t["formula"], "adduct_official": adduct,
                "adduct_mzdesk": CANONICAL_ADDUCTS[adduct],
                "inchikey": t["inchikey"], "inchikey14": t["inchikey14"],
                "smiles_official": t["smiles"], "monoisotopic_mass": float(t["monoisotopic_mass"]),
                "expected_precursor_mz": float(t["expected_mz"]),
                "given_precursor_error_ppm": float(t["ppm_err"]),
                "structure_notes": t["notes"].split(";") if t["notes"] else [],
            },
            "provenance": {
                "priority": m["priority"], "file": m["file"], "mzml_path": m["mzml_path"],
                "rt_min_given": float(m["rt_min"]), "precursor_mz_given": float(m["precursor_mz"]),
                "collection": m["collection"], "rack_id": m["rack_id"], "well_id": m["well_id"],
                "mix_or_set": m["mix_or_set"],
                "extraction_status": r["status"], "extraction_notes": r["notes"],
                "ms1_rt_min": float(r["ms1_rt_used"]) if r["ms1_rt_used"] else None,
                "ms2_rt_min": [float(x) for x in r["ms2_rt_used"].split(";") if x],
                "ms2_precursor_error_ppm": float(r["ms2_precursor_ppm_error"]) if r["ms2_precursor_ppm_error"] else None,
                "ms2_rt_offset_min": float(r["ms2_rt_offset_min"]) if r["ms2_rt_offset_min"] else None,
                "parentmass_source": r["parentmass_source"],
                "extracted_ms_path": str(path.relative_to(BASE)),
                "extracted_ms_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            },
            "ms2_by_collision_energy": [{"collision_energy": b["collision_energy"], "peaks": b["peaks"]} for b in ms2],
        })
    return records


def main():
    records = build_records()  # validate everything before touching the output
    with OUTPUT.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=False, separators=(",", ":"), allow_nan=False) + "\n")
    print(f"Wrote {len(records)} challenges to {OUTPUT}")


if __name__ == "__main__":
    main()
