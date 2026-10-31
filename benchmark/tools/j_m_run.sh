#!/usr/bin/env bash
# JESTR + MIST on all CASMI 2022 challenges (blind, same candidates as CSI:FingerID), with progress bars.
#
# Usage (from anywhere):
#   bash j_m_run.sh                  # all 500 challenges
#   bash j_m_run.sh 3 4 26 163       # only these challenge IDs (quick test)
#
# Env overrides:
#   JOBS=6        parallel worker processes (machine: 10 logical CPUs; each worker uses 1 thread)
#   OUT=...       output dir (default: casmi2022/jm_run)
#   STEPS="prep jestr mist score"   run only some steps
#   FORCE=1       delete previous outputs in OUT and start from scratch
#   MIST_CHUNK=10 challenges per MIST CLI call
#
# Resume: every challenge is written atomically; re-running the script skips finished challenges.
#
# Settings (see tools/jm/*.py headers for details)
#   spectra    : casmi2022/sirius_run/ms (blind extraction, measured precursor m/z); the 3 collision
#                energies are merged (10 ppm, max intensity) for JESTR; MIST merges them itself.
#   candidates : SIRIUS 6.3.5 PubChem candidates (structure_identifications_all.tsv), dedup by
#                InChIKey 2D block -> exactly the list CSI:FingerID ranked. Formula/adduct from SIRIUS.
#   JESTR      : NPLIB1 weights, 74-dim atom features, intensity -> 999, m/z <= 1000, cosine; all 500.
#   MIST       : canopus pretrained FP + contrastive models, subformulae 20 ppm / 50 formulae
#                (MIST training settings), Morgan4096 r2, cosine; positive mode only (304),
#                because the released weights only know positive adducts.
#   scoring    : InChIKey 2D block; optimistic rank under ties (pessimistic reported too);
#                no output = not found; subsets all / Priority / Bonus / pos / neg.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"            # benchmark/
TOOLS="$ROOT/tools/jm"
SIRIUS_RUN="$ROOT/casmi2022/sirius_run"
MS_DIR="$SIRIUS_RUN/ms"
SUMMARIES="$SIRIUS_RUN/results/sirius_summaries"
TRUTH="$SIRIUS_RUN/results/truth.tsv"
JESTR_REPO="$ROOT/JESTR"
MIST_REPO="$ROOT/MIST"
JPY="$JESTR_REPO/.venv/bin/python"
MPY="$MIST_REPO/.venv/bin/python"

JOBS="${JOBS:-6}"
OUT="${OUT:-$ROOT/casmi2022/jm_run}"
STEPS="${STEPS:-prep jestr mist score}"
MIST_CHUNK="${MIST_CHUNK:-10}"
PROG="$OUT/.progress"
LOGS="$OUT/logs"

IDS=""
for i in "$@"; do IDS="${IDS:+$IDS,}$((10#$i))"; done

# ------------------------------------------------------------------ checks
fail() { echo "ERROR: $*" >&2; exit 2; }
for f in "$JPY" "$MPY" "$TRUTH" "$SUMMARIES/structure_identifications_all.tsv" \
         "$JESTR_REPO/data/NPLIB1/pretrained_mol_enc_model_1707829192911_best.pt" \
         "$JESTR_REPO/data/NPLIB1/pretrained_spec_enc_model_1707829192911_best.pt" \
         "$MIST_REPO/pretrained_models/mist_fp_canopus_pretrain.ckpt" \
         "$MIST_REPO/pretrained_models/mist_contrastive_canopus_pretrain.ckpt"; do
  [[ -e "$f" ]] || fail "missing $f"
done
ls "$MS_DIR"/challenge_*.ms >/dev/null 2>&1 || fail "no .ms files in $MS_DIR (run sirius_run/scripts/01_extract_ms.py)"
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || fail "JOBS must be a positive integer"

if [[ "${FORCE:-0}" == 1 ]]; then rm -rf "$OUT"; fi
mkdir -p "$OUT" "$PROG" "$LOGS"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 DGLBACKEND=pytorch PYTHONWARNINGS=ignore

# ------------------------------------------------------------------ progress bar
BAR_W=30
IS_TTY=0; [[ -t 2 ]] && IS_TTY=1
ALL_PIDS=""

cleanup() {
  local p
  for p in $ALL_PIDS; do pkill -TERM -P "$p" 2>/dev/null || true; kill "$p" 2>/dev/null || true; done
}
trap 'echo; echo "interrupted: stopping workers (finished challenges are kept)" >&2; cleanup; exit 130' INT TERM

fmt_t() { printf '%02d:%02d:%02d' $(($1 / 3600)) $(($1 % 3600 / 60)) $(($1 % 60)); }

# draw <step> <label> <unit> <t0>: sums "<units done> <units total> <items done> <items total>" of all shards
draw() {
  local step=$1 label=$2 unit=$3 t0=$4 el=$((SECONDS - $4)) line
  line=$( { cat "$PROG/$step".* 2>/dev/null || true; } | awk -v w="$BAR_W" -v el="$el" -v lab="$label" -v unit="$unit" '
    {d+=$1; t+=$2; i+=$3; n+=$4}
    END {
      f = (t > 0) ? d / t : 0; if (f > 1) f = 1
      k = int(f * w + 0.5); bar = ""
      for (j = 0; j < k; j++) bar = bar "█"
      for (; j < w; j++) bar = bar "░"
      if (f > 0.002 && f < 1) { e = int(el * (1 - f) / f); eta = sprintf("%02d:%02d:%02d", e/3600, (e%3600)/60, e%60) }
      else if (f >= 1) eta = "00:00:00"; else eta = "--:--:--"
      printf "%-6s %s %5.1f%%  %d/%d %s  ETA %s", lab, bar, f * 100, i, n, unit, eta
    }')
  line="$line  elapsed $(fmt_t "$el")"
  if [[ $IS_TTY == 1 ]]; then printf '\r\033[K%s' "$line" >&2; else echo "$line" >&2; fi
}

# monitor <step> <label> <unit> <pid...>: redraws until all workers end; stops early if one fails
monitor() {
  local step=$1 label=$2 unit=$3; shift 3
  local left="$*" t0=$SECONDS last=0 p rc still
  while [[ -n "$left" ]]; do
    still=""
    for p in $left; do
      if kill -0 "$p" 2>/dev/null; then still="$still $p"; continue; fi
      rc=0; wait "$p" || rc=$?
      if [[ $rc -ne 0 ]]; then
        draw "$step" "$label" "$unit" "$t0"; echo >&2
        echo "ERROR: $label worker (pid $p) exited with $rc. Last log lines:" >&2
        tail -n 25 "$(cat "$PROG/pid.$p")" >&2 || true
        cleanup; exit 1
      fi
    done
    left="$still"
    ls "$PROG/$step".* >/dev/null 2>&1 || { sleep 1; continue; }   # workers still starting up
    if [[ $IS_TTY == 1 || $((SECONDS - last)) -ge 30 || -z "$left" ]]; then
      draw "$step" "$label" "$unit" "$t0"; last=$SECONDS
    fi
    [[ -n "$left" ]] && sleep 1
  done
  [[ $IS_TTY == 1 ]] && echo >&2
  return 0
}

# launch <log> <cmd...>: start in background, remember pid -> log
launch() {
  local log=$1; shift
  "$@" >"$log" 2>&1 &
  local p=$!
  echo "$log" >"$PROG/pid.$p"
  ALL_PIDS="$ALL_PIDS $p"
  LAUNCHED="${LAUNCHED:-} $p"
}

has_step() { [[ " $STEPS " == *" $1 "* ]]; }
NSTEP=0; for s in $STEPS; do NSTEP=$((NSTEP + 1)); done
K=0
header() { K=$((K + 1)); echo "[$(date +%T)] [$K/$NSTEP] $*" >&2; }

T_ALL=$SECONDS
echo "JESTR + MIST on CASMI 2022 | out=$OUT | jobs=$JOBS | ids=${IDS:-all}" >&2

# ------------------------------------------------------------------ 1. inputs
if has_step prep; then
  header "prepare inputs (spectra, SIRIUS PubChem candidates, MIST hypotheses)"
  STAMP="$OUT/inputs/.done"
  if [[ -f "$STAMP" && "$(cat "$STAMP")" == "ids=${IDS:-all}" ]]; then
    echo "       inputs up to date (ids=${IDS:-all}), skipped" >&2
  else
    rm -rf "$OUT/inputs" "$PROG"/prep.*
    LAUNCHED=""
    launch "$LOGS/prep.log" "$JPY" "$TOOLS/prep.py" --ms-dir "$MS_DIR" --summaries "$SUMMARIES" \
      --out "$OUT/inputs" --ids "$IDS" --progress-file "$PROG/prep.0"
    monitor prep "prep" "challenges" $LAUNCHED
    echo "ids=${IDS:-all}" >"$STAMP"
    echo "       $(tail -n 1 "$LOGS/prep.log")" >&2
  fi
fi
[[ -f "$OUT/inputs/.done" ]] || fail "inputs missing: run with STEPS containing 'prep'"

# ------------------------------------------------------------------ 2. JESTR
if has_step jestr; then
  header "JESTR ranking (all polarities, $JOBS workers)"
  rm -f "$PROG"/jestr.*
  LAUNCHED=""
  for ((s = 0; s < JOBS; s++)); do
    launch "$LOGS/jestr_shard$s.log" "$JPY" "$TOOLS/run_jestr.py" --inputs "$OUT/inputs" --out "$OUT/jestr" \
      --repo "$JESTR_REPO" --shard "$s" --nshards "$JOBS" --threads 1 --progress-file "$PROG/jestr.$s"
  done
  monitor jestr "JESTR" "challenges" $LAUNCHED
fi

# ------------------------------------------------------------------ 3. MIST
if has_step mist; then
  header "MIST fingerprint + contrastive ranking (positive mode, $JOBS workers)"
  rm -f "$PROG"/mist.*
  LAUNCHED=""
  for ((s = 0; s < JOBS; s++)); do
    launch "$LOGS/mist_shard$s.log" "$MPY" "$TOOLS/run_mist.py" --inputs "$OUT/inputs" --out "$OUT/mist" \
      --repo "$MIST_REPO" --shard "$s" --nshards "$JOBS" --chunk "$MIST_CHUNK" --progress-file "$PROG/mist.$s"
  done
  monitor mist "MIST" "challenges" $LAUNCHED
fi

# ------------------------------------------------------------------ 4. score
if has_step score; then
  header "score vs CASMI truth (InChIKey 2D block)"
  "$JPY" "$TOOLS/score.py" --run-dir "$OUT" --truth "$TRUTH" --out "$OUT/score" 2>&1 | tee "$LOGS/score.log" >&2
fi

rm -f "$PROG"/pid.*
echo "[$(date +%T)] done in $(fmt_t $((SECONDS - T_ALL))) | results: $OUT/score/{summary,per_challenge}.tsv | logs: $LOGS" >&2
