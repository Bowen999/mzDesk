#!/usr/bin/env python3
"""Export the official CASMI 2022 final key and extracted spectra as CSV.

Each challenge has its official fields, an MS1 isotope peak list and a merged
MS2 peak list. Analysis settings are not included.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


BASE = Path(__file__).resolve().parent
OUTPUT = BASE / "casmi2022_ground_truth.csv"
CANONICAL_ADDUCTS = {
    "[M+H]+": "[M+H]+", "[M-H]-": "[M-H]-", "[M+NH4]+": "[M+NH4]+",
    "[M+Na]+": "[M+Na]+", "[M+Cl]-": "[M+Cl]-", "[M+FA-H]-": "[M+HCOO]-",
}

# Official fields map directly to columns in the final key.
OFFICIAL_COLUMNS: tuple[tuple[str, str], ...] = (
    ("challenge_id", "int"),
    ("formula", "str"),
    ("adduct_official", "str"),
    ("inchikey", "str"),
    ("smiles_official", "str"),
    ("monoisotopic_mass", "float"),
    ("priority", "str"),
    ("file", "str"),
    ("rt_min_given", "float"),
    ("precursor_mz_given", "float"),
    ("collection", "str"),
    ("rack_id", "str"),
    ("well_id", "str"),
    ("mix_or_set", "str"),
)
# Peak lists use one newline-separated ``m/z intensity`` pair per line.
COLUMNS = OFFICIAL_COLUMNS + (("ms1_isotope", "str"), ("ms2", "str"))
KIND = dict(COLUMNS)
NAMES = [name for name, _ in COLUMNS]
OFFICIAL_HEADERS = {
    "challenge_id": "Compound Number", "formula": "Formula", "adduct_official": "Adduct",
    "inchikey": "InChIKey", "smiles_official": "SMILES", "monoisotopic_mass": "Monoisotopic Mass",
    "priority": "Priority/Bonus", "file": "File", "rt_min_given": "RT [min]",
    "precursor_mz_given": "Precursor m/z (Da)", "collection": "Collection",
    "rack_id": "Rack ID or UCD Barcode", "well_id": "Well ID", "mix_or_set": "Mix or Set Number",
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


def encode_cell(value, kind, name):
    """Serialise one value; ``str`` of a float round-trips through ``float``."""
    if value is None:
        if kind.endswith("?"):
            return ""
        raise ValueError(f"column {name}: missing value")
    if kind in ("int", "float", "float?"):
        return str(value)
    if kind == "bool":
        return "true" if value else "false"
    if kind in ("list", "floatlist"):
        return ";".join(str(item) for item in value)
    if kind == "json":
        return json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    return value


def decode_cell(text, kind, name):
    if text is None:
        raise ValueError(f"column {name}: row has too few fields")
    if kind == "int":
        return int(text)
    if kind == "bool":
        if text not in ("true", "false"):
            raise ValueError(f"column {name}: {text!r} is not a boolean")
        return text == "true"
    if kind in ("float", "float?"):
        if not text:
            if kind == "float?":
                return None
            raise ValueError(f"column {name}: empty numeric cell")
        return float(text)
    if kind == "list":
        return text.split(";") if text else []
    if kind == "floatlist":
        return [float(item) for item in text.split(";")] if text else []
    if kind == "json":
        return json.loads(text)
    return text


def flatten(record):
    """Serialise one record in ``NAMES`` order."""
    cid = record.get("challenge_id")
    if set(record) != set(NAMES):
        raise ValueError(f"challenge {cid}: fields differ from the columns: {sorted(set(record) ^ set(NAMES))}")
    return {name: encode_cell(record[name], KIND[name], name) for name in NAMES}


def unflatten(row):
    """Decode one CSV row."""
    return {name: decode_cell(row[name], kind, name) for name, kind in COLUMNS}


def write_records(records, path=OUTPUT):
    """Write the records as CSV."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(NAMES)
        for record in records:
            writer.writerow(flatten(record).values())


def load_records(path=OUTPUT):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != NAMES:
            raise ValueError(f"{path}: header differs from the column list")
        rows = list(reader)
    if any(row.get(None) for row in rows):
        raise ValueError(f"{path}: a row has the wrong number of fields")
    return [unflatten(row) for row in rows]


def build_extracted_records():
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
                or len(ms1) != 1 or not ms1[0]["peaks"]
                or len(ms1[0]["peaks"]) != int(r["n_ms1_isotope_peaks"])
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


def build_records():
    """Combine the official final answer workbook and extracted peak lists."""
    from openpyxl import load_workbook

    source = BASE / "metadata/source/casmi2022_key_final.xlsx"
    book = load_workbook(source, read_only=True, data_only=True)
    try:
        rows = book["CHALLENGES & KEY ALL"].iter_rows(values_only=True)
        next(rows)  # Share / DO NOT SHARE banner
        headers = next(rows)
        if len(set(headers)) != len(headers) or set(OFFICIAL_HEADERS.values()) != set(headers):
            raise ValueError(f"{source}: unexpected official key columns")
        indices = {name: headers.index(header) for name, header in OFFICIAL_HEADERS.items()}
        records = []
        for row in rows:
            if row[indices["challenge_id"]] is None:
                continue
            record = {}
            for name, kind in OFFICIAL_COLUMNS:
                value = row[indices[name]]
                if value is None:
                    raise ValueError(f"{source}: missing {name} for challenge {row[indices['challenge_id']]}")
                if kind == "int":
                    if int(value) != value:
                        raise ValueError(f"{source}: invalid challenge ID {value!r}")
                    value = int(value)
                elif kind == "float":
                    value = float(value)
                    if not math.isfinite(value):
                        raise ValueError(f"{source}: non-finite {name}")
                else:
                    # Some official identifiers are stored as integral Excel numbers.
                    if name in ("rack_id", "mix_or_set") and isinstance(value, (int, float)) and value == int(value):
                        value = int(value)
                    value = str(value).strip()
                    if not value:
                        raise ValueError(f"{source}: empty {name}")
                record[name] = value
            records.append(record)
        records.sort(key=lambda record: record["challenge_id"])
        if [record["challenge_id"] for record in records] != list(range(1, 501)):
            raise ValueError(f"{source}: expected each challenge ID from 1 through 500 exactly once")
    finally:
        book.close()

    extracted = {record["challenge_id"]: record for record in build_extracted_records()}
    for record in records:
        cid = record["challenge_id"]
        spectrum = extracted[cid]
        truth, provenance = spectrum["ground_truth"], spectrum["provenance"]
        for name in ("formula", "adduct_official", "inchikey", "smiles_official"):
            if record[name] != truth[name]:
                raise ValueError(f"challenge {cid}: {name} differs between key and spectrum metadata")
        for name in ("priority", "file"):
            if record[name] != provenance[name]:
                raise ValueError(f"challenge {cid}: {name} differs between key and spectrum metadata")
        for name in ("rt_min_given", "precursor_mz_given"):
            if abs(record[name] - provenance[name]) > 1e-7:
                raise ValueError(f"challenge {cid}: {name} differs between key and spectrum metadata")
        record["ms1_isotope"] = spectrum["inputs"]["ms1"]
        record["ms2"] = spectrum["inputs"]["ms2"]["value"]
    return records


def main():
    records = build_records()  # validate everything before touching the output
    write_records(records)
    print(f"Wrote {len(records)} challenges to {OUTPUT}")


if __name__ == "__main__":
    main()
