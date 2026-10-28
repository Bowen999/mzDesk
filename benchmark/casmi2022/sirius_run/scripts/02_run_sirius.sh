#!/usr/bin/env bash
# Blind SIRIUS 6.3.5 run on the CASMI 2022 .ms files (CSI:FingerID + PubChem, MSNovelist).
#
# Usage:
#   scripts/02_run_sirius.sh                 # all ms/challenge_*.ms
#   scripts/02_run_sirius.sh 1 2 3 450       # only these challenge IDs
#   OUT_TAG=pilot scripts/02_run_sirius.sh 1 2 3 4 18 253 300 450   # pilot (separate project/logs/results)
#
# Env overrides: SIRIUS (binary), CORES (default 8), CHUNK (compounds per SIRIUS call, default 25),
#                OUT_TAG (empty = full run), SKIP_EXPORT=1 (do not write summaries).
#
# Outputs (RUN = sirius_run/):
#   full : RUN/project/casmi2022.sirius, RUN/logs/, RUN/results/sirius_summaries/
#   tag  : RUN/project_<tag>/casmi2022.sirius, RUN/logs/<tag>/, RUN/results/<tag>/sirius_summaries/
#
# Resume / re-run: the .sirius project is kept. A ledger (imported_ids.txt next to the project)
# lists the IDs already imported, so they are never imported twice (SIRIUS itself does NOT dedupe:
# importing the same .ms again creates a duplicate compound). Each call first imports+computes new
# chunks, then a final "-i project -o project" pass completes any instance whose tool results are
# missing (SIRIUS skips tools whose results already exist unless --recompute is given).
# To redo a compound from scratch, delete the project dir and the ledger (or the whole project).
# Keep the settings below constant across resumes; changing them on an existing project makes
# SIRIUS recompute/invalidate results.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN="$(cd "$HERE/.." && pwd)"
SIRIUS="${SIRIUS:-~/tools/sirius.app/Contents/MacOS/sirius}"
CORES="${CORES:-8}"          # machine: 10 logical CPUs (5 physical); leave 2 free
CHUNK="${CHUNK:-25}"
OUT_TAG="${OUT_TAG:-}"

if [[ -n "$OUT_TAG" ]]; then
  PROJ_DIR="$RUN/project_$OUT_TAG"; LOG_DIR="$RUN/logs/$OUT_TAG"; SUM_DIR="$RUN/results/$OUT_TAG/sirius_summaries"
else
  PROJ_DIR="$RUN/project"; LOG_DIR="$RUN/logs"; SUM_DIR="$RUN/results/sirius_summaries"
fi
PROJECT="$PROJ_DIR/casmi2022.sirius"
LEDGER="$PROJ_DIR/imported_ids.txt"
mkdir -p "$PROJ_DIR" "$LOG_DIR"
touch "$LEDGER"

# ------------------------------------------------------------------ settings
# Generic for every compound; nothing here depends on a compound or its answer.
# Adducts: enforced = always considered (polarity is taken from '>ionization [M+?]+/-' in the .ms,
# so SIRIUS only uses the matching-polarity entries). Names are SIRIUS 6 spellings:
#   [M+NH4]+ is reported as [M+H3N+H]+, [M+FA-H]- as [M+CH2O2-H]-, [M+HAc-H]- as [M+C2H4O2-H]-.
ADDUCTS_POS='[M+H]+,[M+Na]+,[M+NH4]+,[M+K]+,[M+H-H2O]+'
ADDUCTS_NEG='[M-H]-,[M+Cl]-,[M+FA-H]-,[M+HAc-H]-,[M-H2O-H]-'
ADDUCTS="$ADDUCTS_POS,$ADDUCTS_NEG"
CFG=(
  --AlgorithmProfile=orbitrap                      # Thermo Q Exactive HF Orbitrap
  --MS1MassDeviation.allowedMassDeviation=10ppm    # measured MS1 m/z is used as parentmass, up to ~8 ppm off
  --MS2MassDeviation.allowedMassDeviation=10ppm    # MS2 at 15k resolution
  "--AdductSettings.enforced=$ADDUCTS"             # always consider these
  "--AdductSettings.fallback=$ADDUCTS"             # (same list, if nothing is detected)
  --AdductSettings.ignoreDetectedAdducts=false
  --FormulaSettings.enforced=C,H,N,O,P,F,I         # = SIRIUS "extended organic" (adds F, I)
  --FormulaSettings.detectable=S,Br,Cl,B,Se,Fe,Zn,Mg  # added when MS1 isotope pattern indicates them
  --FormulaSettings.fallback=S
  --StructureSearchDB=PUBCHEM                      # CSI:FingerID candidate DB
  --NumberOfCandidates=10                          # formula candidates kept (default)
  --NumberOfCandidatesPerIonization=1              # >=1 formula per adduct kept (default)
  --NumberOfStructureCandidates=10000              # default; per-compound structure list cap
  --NumberOfMsNovelistCandidates=128               # MSNovelist (default)
  --Timeout.secondsPerInstance=300                 # default per-compound timeout
)
# Tool chain. CANOPUS ('classes') is part of the CLI chain only because 'structures' must follow it.
CHAIN=( formulas -c 10 fingerprints classes structures -d PUBCHEM denovo-structures -c 128 )

# ------------------------------------------------------------------ select IDs
if [[ $# -gt 0 ]]; then
  IDS=(); for i in "$@"; do IDS+=("$(printf '%03d' "$((10#$i))")"); done
else
  IDS=(); for f in "$RUN"/ms/challenge_*.ms; do b="$(basename "$f" .ms)"; IDS+=("${b#challenge_}"); done
fi
NEW=()
for id in "${IDS[@]}"; do
  [[ -f "$RUN/ms/challenge_$id.ms" ]] || { echo "missing ms/challenge_$id.ms" >&2; exit 2; }
  grep -qx "$id" "$LEDGER" || NEW+=("$id")
done
echo "[$(date +%T)] project=$PROJECT  requested=${#IDS[@]}  new=${#NEW[@]}  cores=$CORES"
"$SIRIUS" --version 2>/dev/null | grep -E '^(SIRIUS|CSI)' | tee "$LOG_DIR/versions.txt" || true

run_sirius() {  # run_sirius <logfile> <sirius args...>
  local log="$1"; shift
  local t0=$SECONDS
  "$SIRIUS" --cores "$CORES" "$@" >"$log" 2>&1
  echo "  exit 0, $((SECONDS - t0)) s, log: ${log#$RUN/}"
}

# ------------------------------------------------------------------ import + compute, chunked
n=0; k=0
while [[ $n -lt ${#NEW[@]} ]]; do
  chunk=("${NEW[@]:n:CHUNK}"); n=$((n + CHUNK)); k=$((k + 1))
  stage="$LOG_DIR/stage_$k"; rm -rf "$stage"; mkdir -p "$stage"
  for id in "${chunk[@]}"; do cp "$RUN/ms/challenge_$id.ms" "$stage/"; done
  echo "[$(date +%T)] chunk $k: ${chunk[0]}..${chunk[$((${#chunk[@]}-1))]} (${#chunk[@]} compounds)"
  run_sirius "$LOG_DIR/chunk_$(printf '%03d' $k).log" \
    -i "$stage" -o "$PROJECT" config "${CFG[@]}" "${CHAIN[@]}"
  printf '%s\n' "${chunk[@]}" >> "$LEDGER"     # only after SIRIUS exited 0
  rm -rf "$stage"
done

# ------------------------------------------------------------------ complete unfinished instances
if [[ -s "$LEDGER" ]]; then
  echo "[$(date +%T)] completion pass over the project (skips finished compounds)"
  run_sirius "$LOG_DIR/complete_pass.log" -i "$PROJECT" -o "$PROJECT" config "${CFG[@]}" "${CHAIN[@]}"
fi

# ------------------------------------------------------------------ export
if [[ "${SKIP_EXPORT:-0}" != 1 ]]; then
  tmp="$SUM_DIR.tmp"; rm -rf "$tmp"; mkdir -p "$(dirname "$SUM_DIR")"
  echo "[$(date +%T)] exporting summaries"
  # --full-summary: ALL candidates (structure_identifications_all.tsv, denovo_structure_identifications_all.tsv,
  #                 formula_identifications_all.tsv, ...). --top-hit-summary: rank-1 per compound incl. confidence.
  # --top-hit-adduct-summary: best formula per adduct. TSV.
  run_sirius "$LOG_DIR/export.log" -i "$PROJECT" -o "$PROJECT" \
    summaries --full-summary --top-hit-summary --top-hit-adduct-summary --format=TSV --output "$tmp"
  rm -rf "$SUM_DIR"; mv "$tmp" "$SUM_DIR"
  echo "[$(date +%T)] summaries in $SUM_DIR"
  ls -l "$SUM_DIR" | awk 'NR>1{print "   ",$5,$9}'
  # QC: which requested compounds have formula / PubChem-structure / MSNovelist results (missing ones are
  # flagged 0 in qc_results_present.tsv; e.g. "No suitable fragmentation tree" -> no structure, or timeout, or MSNovelist mass limit).
  python3 - "$SUM_DIR" "${IDS[@]}" <<'PY'
import sys, re, csv, pandas as pd
d, ids = sys.argv[1], ["challenge_" + i for i in sys.argv[2:]]
def have(fn):
    # quoting=QUOTE_NONE is essential: compound names contain a stray '"' that would otherwise make
    # pandas swallow following lines as one quoted field and silently drop whole challenges (v1 bug).
    out = set()
    try:
        for ch in pd.read_csv(f"{d}/{fn}", sep="\t", usecols=["mappingFeatureId"], quoting=csv.QUOTE_NONE,
                              dtype=str, chunksize=500000):
            out |= set(ch.mappingFeatureId.dropna().unique())
    except Exception: return set()
    return {m.group(1) for m in (re.search(r"_(challenge_\d+)_challenge", x) for x in out) if m}
sets = {"formula": have("formula_identifications_all.tsv"),
        "pubchem_structure": have("structure_identifications_all.tsv"),
        "msnovelist": have("denovo_structure_identifications_all.tsv")}
rows = [(c, *[int(c in v) for v in sets.values()]) for c in ids]
q = pd.DataFrame(rows, columns=["challenge", *sets])
q.to_csv(f"{d}/qc_results_present.tsv", sep="\t", index=False)
print("   QC (of %d requested): " % len(ids) + ", ".join("%s=%d" % (k, q[k].sum()) for k in sets))
PY
fi
echo "[$(date +%T)] done"
