#!/usr/bin/env python
"""FRIGID-base driver: positive-mode MS/MS + known molecular formula -> ranked SMILES candidates.

Run with the dedicated env:   ~/tools/frigid/env/bin/python frigid_predict.py --input job.json --output result.json

Pipeline (all FRIGID's own code, nothing re-implemented except the subformula assignment):
  peaks --(subformula_assign.py, MIST-style, 20 ppm / top-50 / 0.001)--> subformula JSON
        --(mist.data.featurizers.PeakFormula, max_peaks=10, cls=ms1)--> batch
        --(mist_msg_aug.ckpt, ScaledMIST)--> 4096-bit Morgan fingerprint probabilities --(>= 0.187)--> binary FP
        --(dlm.sampler.Sampler + dlm.utils.benchmark_utils.generate_with_formula_filter, DLM.ckpt)--> SMILES
  score = Tanimoto(Morgan r=2 4096 of candidate, MIST-predicted binary FP)   [FRIGID's own similarity]
  order = FRIGID's own: formula-matched (unique 2D InChIKey) by score desc, then padded non-matched by score desc.

Input JSON (all keys except the optional ones are required):
  {"formula": "C8H10N4O2",        # neutral molecular formula, elements limited to C H N O P S F Cl Br I B Si Se As Fe Co Na K
   "precursor_mz": 195.0877,
   "adduct": "[M+H]+",            # [M+H]+ [M+Na]+ [M+K]+ [M+NH4]+ [M-H2O+H]+ [M-2H2O+H]+ [M]+ ; negative mode is rejected
   "ion_mode": "positive",        # optional, "negative" is rejected
   "peaks": [[mz, intensity], ...],
   "mode": "quick" | "full",      # quick: one batch of 16 samples; full: until 10 unique formula-matched or 100 samples
   "device": "auto" | "mps" | "cpu",  # optional, default "auto": MPS only if available memory >= 12 GB AND swap free > 1 GB at
                                  # job start (numbers + decision are written to the output as `device_decision`), else CPU
   "seed": 42,                    # optional (default 42)
   "length_model": "auto" | "none" | "/path/model.joblib",   # SAFE-token-length model (FRIGID's NGBoost, features JSON next to it).
                                  # "auto" (default) = $FRIGID_TOKEN_MODEL or FRIGID/token_models/models/best_ngboost_MSG.joblib if it
                                  # exists, else FRIGID's built-in fallback (uniform random length 20..99). "none" forces the fallback.
   "sigma_lambda": 3.0,           # optional, FRIGID benchmark default: length ~ N(mean, variance = sigma * sigma_lambda)
   "mps_empty_cache": false,      # optional: torch.mps.empty_cache() after each generation batch (lower memory, MPS only)
   "name": "caffeine",            # optional label
   "truth_smiles": "..."}         # optional, ONLY used to report `true_rank` afterwards; never given to the model

Output JSON: see the dict built at the bottom of `predict()`; candidates[] = {rank, smiles, formula, score,
formula_match, inchikey_2d, isotope_labeled, radical, charged}. Progress/log text goes to stderr; stdout is silent.
Exit code: 0 ok, 2 invalid input (incl. negative mode), 1 runtime failure, 3 memory guard tripped (result JSON is always
written; ok=false for 1/2/3).

Memory guard (driver/memguard.py, runs inside this process): CPU -> process-tree RSS cap 6144 MB; MPS -> process-tree physical
footprint cap 12288 MB (RSS does not see MPS unified memory). Override with env FRIGID_RSS_CAP_MB / FRIGID_FOOTPRINT_CAP_MB.
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import random
try:
    import resource
except ImportError:          # Windows has no resource module
    resource = None
import sys
import tempfile
import time
import traceback
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
# FRIGID's get_tokenizer() falls back to downloading the (46 KB) datamol-io/safe-gpt tokenizer from the HF Hub. Keep that cache
# inside ~/tools/frigid and stay offline once it is there (cache was populated by the very first smoke run).
os.environ.setdefault("HF_HOME", os.path.join(os.path.dirname(HERE), "hf_cache"))
os.environ.setdefault("HF_HUB_OFFLINE", "1")
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
FRIGID_ROOT = os.environ.get("FRIGID_ROOT", os.path.join(os.path.dirname(HERE), "FRIGID"))
MIST_CKPT = os.environ.get("FRIGID_MIST_CKPT", os.path.join(os.path.dirname(HERE), "checkpoints", "mist_msg_aug.ckpt"))
DLM_CKPT = os.environ.get("FRIGID_DLM_CKPT", os.path.join(os.path.dirname(HERE), "checkpoints", "DLM.ckpt"))
DEFAULT_LENGTH_MODEL = os.environ.get("FRIGID_TOKEN_MODEL", os.path.join(FRIGID_ROOT, "token_models", "models", "best_ngboost_MSG.joblib"))
sys.path.insert(0, os.path.join(FRIGID_ROOT, "src"))
sys.path.insert(0, HERE)

# FRIGID settings for the MassSpecGym MIST encoder (configs/spec2mol_benchmark_msg.yaml)
FP_BITS, FP_RADIUS, FP_THRESHOLD = 4096, 2, 0.187
SOFTMAX_TEMP, RANDOMNESS, BATCH_SIZE = 1.0, 10.0, 16
MODES = {"quick": dict(n_required=16, max_attempts=16), "full": dict(n_required=10, max_attempts=100)}
ENCODER_KW = dict(cls_type="ms1", inten_transform="float", magma_modulo=2048, max_peaks=10)


class InputError(ValueError):
    pass


def log(*a):
    print("[frigid_predict]", *a, file=sys.stderr, flush=True)


def validate(job: dict) -> dict:
    import subformula_assign as sa

    j = dict(job)
    for k in ("formula", "precursor_mz", "adduct", "peaks", "mode"):
        if k not in j:
            raise InputError(f"missing key {k!r}")
    adduct = str(j["adduct"]).replace(" ", "")
    mode_txt = str(j.get("ion_mode", "positive")).lower()
    if mode_txt.startswith("neg") or adduct.endswith("-") or adduct.endswith("-]") or adduct.endswith("]-"):
        raise InputError("negative ion mode is not supported (FRIGID/MIST were trained on positive-mode adducts only)")
    if not mode_txt.startswith("pos"):
        raise InputError(f"ion_mode must be 'positive', got {j.get('ion_mode')!r}")
    try:
        j["adduct"] = sa.standardize_adduct(adduct)
        sa.formula_to_dense(str(j["formula"]))
    except ValueError as e:
        raise InputError(str(e)) from None
    j["formula"] = str(j["formula"]).strip()
    j["precursor_mz"] = float(j["precursor_mz"])
    if not j["precursor_mz"] > 0:
        raise InputError("precursor_mz must be > 0")
    peaks = [(float(p[0]), float(p[1])) for p in j["peaks"]]
    peaks = [p for p in peaks if p[0] > 0 and p[1] > 0]
    if not peaks:
        raise InputError("no usable peaks (need [[mz, intensity], ...] with positive values)")
    j["peaks"] = peaks
    if j["mode"] not in MODES:
        raise InputError(f"mode must be one of {list(MODES)}")
    j["device"] = str(j.get("device", "auto")).lower()
    if j["device"] not in ("mps", "cpu", "auto"):
        raise InputError("device must be 'mps', 'cpu' or 'auto'")
    j["seed"] = int(j.get("seed", 42))
    j["mps_empty_cache"] = bool(j.get("mps_empty_cache", False))
    j["sigma_lambda"] = float(j.get("sigma_lambda", 3.0))
    lm = str(j.get("length_model", "auto"))
    if lm not in ("auto", "none") and not os.path.isfile(lm):
        raise InputError(f"length_model file not found: {lm}")
    j["length_model"] = lm
    return j


def load_models(device: str) -> dict:
    """Load MIST encoder + DLM sampler once (reusable across jobs in a long-lived worker)."""
    import torch
    from dlm.sampler import Sampler
    from mist.models.checkpoint import load_mist_checkpoint

    os.environ["FRIGID_DEVICE"] = device  # read by the patched dlm.sampler._default_device()
    t = time.perf_counter()
    mist = load_mist_checkpoint(MIST_CKPT, device=torch.device(device))
    t_mist = time.perf_counter() - t
    t = time.perf_counter()
    sampler = Sampler(DLM_CKPT)
    t_dlm = time.perf_counter() - t
    return dict(mist=mist, sampler=sampler, device=device, t_mist=t_mist, t_dlm=t_dlm)


def predict(job: dict, models: dict) -> dict:
    import numpy as np
    import torch
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import AllChem, rdMolDescriptors

    from dlm.utils.benchmark_utils import (binarize_fingerprint, generate_with_formula_filter, load_token_model,
                                           normalize_formula, predict_token_count)
    from mist.data.data import Spectra
    from mist.data.featurizers import PeakFormula
    import subformula_assign as sa

    RDLogger.DisableLog("rdApp.*")
    device, mist, sampler = models["device"], models["mist"], models["sampler"]
    dev = torch.device(device)
    out_warn: list[str] = []

    # ---- seeds -------------------------------------------------------------------------------------------------
    seed = job["seed"]
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if device == "mps":
        torch.mps.manual_seed(seed)

    # ---- subformula assignment + featurisation (CPU) ------------------------------------------------------------
    t0 = time.perf_counter()
    sub = sa.assign_subformulae(job["formula"], job["adduct"], job["peaks"], job["precursor_mz"],
                                ppm=20.0, max_peaks=50, inten_thresh=0.001)
    tbl = sub["output_tbl"]
    n_assigned = 0 if tbl is None else len(tbl["formula"])
    if n_assigned == 0:
        out_warn.append("no peak could be assigned a sub-formula; the encoder only sees the precursor token")
    elif n_assigned < 3:
        out_warn.append(f"only {n_assigned} peaks assigned a sub-formula; fingerprint prediction will be weak")
    with tempfile.TemporaryDirectory(prefix="frigid_subform_") as tmp:
        with open(os.path.join(tmp, "query.json"), "w") as f:
            json.dump(sub, f)
        feat = PeakFormula(subform_folder=tmp, **ENCODER_KW)
        spec = Spectra(spectra_name="query", instrument="")
        batch = PeakFormula.collate_fn([feat.featurize(spec)])
    batch = {k: v.to(dev) if isinstance(v, torch.Tensor) else v for k, v in batch.items()}
    t_subform = time.perf_counter() - t0

    # ---- MIST encoder -> predicted fingerprint --------------------------------------------------------------------
    t0 = time.perf_counter()
    with torch.no_grad():
        probs, _ = mist(batch)
    probs = probs.float().cpu().numpy()[0]
    if device == "mps":
        torch.mps.synchronize()
    fp_bin = binarize_fingerprint(probs, FP_THRESHOLD)
    t_enc = time.perf_counter() - t0
    if fp_bin.sum() == 0:
        out_warn.append("predicted fingerprint has no bit above threshold; generation is effectively formula-only")

    # ---- DLM generation with FRIGID's own formula filter ------------------------------------------------------------
    target_formula = normalize_formula(job["formula"])
    mps_peak = {"driver": 0, "n_batches": 0}
    orig_generate = sampler.generate
    if device == "mps":
        def _generate_tracked(*a, **k):  # per-batch hook only; FRIGID's generate() itself is untouched
            try:
                return orig_generate(*a, **k)
            finally:
                mps_peak["n_batches"] += 1
                mps_peak["driver"] = max(mps_peak["driver"], torch.mps.driver_allocated_memory())
                if job.get("mps_empty_cache"):  # optional: release the MPS allocator pool after every batch
                    torch.mps.empty_cache()
        sampler.generate = _generate_tracked
    cfg = MODES[job["mode"]]
    # SAFE-token-length model: passed explicitly to FRIGID's generate_with_formula_filter (same as benchmark_spec2mol.py --token-model)
    lm_path = {"none": None, "auto": DEFAULT_LENGTH_MODEL}.get(job["length_model"], job["length_model"])
    token_model = token_features = None
    is_ngb = False
    length_info = dict(used=False, path=None)
    if lm_path and os.path.isfile(lm_path):
        token_model, token_features, is_ngb = load_token_model(lm_path)
        mean_len, sigma = predict_token_count(token_model, token_features, target_formula, is_ngb, job["sigma_lambda"])
        if mean_len is not None:
            length_info = dict(used=True, path=lm_path, type="NGBoost" if is_ngb else "sklearn", features=list(token_features),
                               predicted_total_tokens=mean_len, sigma=None if sigma is None else round(float(sigma), 4),
                               sigma_lambda=job["sigma_lambda"],
                               sampling_std=None if sigma is None else round(float((sigma * job["sigma_lambda"]) ** 0.5), 3))
        else:
            out_warn.append("length model could not predict for this formula; FRIGID's fallback lengths are used")
    elif job["length_model"] == "auto":
        out_warn.append(f"no length model at {lm_path}; using FRIGID's fallback lengths (uniform random 20..99 tokens)")
    t0 = time.perf_counter()
    (smiles, sims, total_gen, total_valid, total_matched, _counter, _last, gen_time) = generate_with_formula_filter(
        sampler=sampler, fingerprint_array=fp_bin, target_formula=target_formula, target_smiles="",
        n_required=cfg["n_required"], max_attempts=cfg["max_attempts"], batch_size=BATCH_SIZE,
        softmax_temp=SOFTMAX_TEMP, randomness=RANDOMNESS, fp_bits=FP_BITS, fp_radius=FP_RADIUS,
        token_model=token_model, token_features=token_features, is_ngboost=is_ngb, sigma_lambda=job["sigma_lambda"])
    if device == "mps":
        torch.mps.synchronize()
        sampler.generate = orig_generate
    t_gen_wall = time.perf_counter() - t0

    # ---- build candidate list -----------------------------------------------------------------------------------------
    truth_key = None
    if job.get("truth_smiles"):
        m = Chem.MolFromSmiles(job["truth_smiles"])
        truth_key = Chem.MolToInchiKey(m).split("-")[0] if m is not None else None
    cands, true_rank = [], None
    for i, (smi, sim) in enumerate(zip(smiles, sims), start=1):
        mol = Chem.MolFromSmiles(smi)
        formula = normalize_formula(rdMolDescriptors.CalcMolFormula(mol)) if mol is not None else None
        key = Chem.MolToInchiKey(mol).split("-")[0] if mol is not None else None
        atoms = list(mol.GetAtoms()) if mol is not None else []
        cands.append(dict(rank=i, smiles=smi, formula=formula, score=round(float(sim), 6),
                          formula_match=(formula == target_formula), inchikey_2d=key,
                          # FRIGID's formula check (RDKit CalcMolFormula) ignores isotope labels, so e.g. a [3H]- or
                          # [11C]-labelled molecule can count as a formula match; flag these for the caller to filter.
                          isotope_labeled=any(a.GetIsotope() != 0 for a in atoms),
                          radical=any(a.GetNumRadicalElectrons() != 0 for a in atoms),
                          charged=any(a.GetFormalCharge() != 0 for a in atoms)))
        if truth_key and key == truth_key and true_rank is None:
            true_rank = i
    n_match = sum(c["formula_match"] for c in cands)
    if n_match == 0:
        out_warn.append("no generated molecule matches the requested formula; candidates are formula-mismatched padding")
    elif n_match < cfg["n_required"] and job["mode"] == "full":
        out_warn.append(f"only {n_match} formula-matched unique molecules found in {total_gen} samples")

    if resource is not None:
        mem = dict(ru_maxrss_mb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576, 1))  # bytes on macOS
    else:                                                                    # Windows: the kernel's peak working set of this process
        import memguard as _mg
        mem = dict(ru_maxrss_mb=round((_mg.proc_memory(os.getpid()) or {}).get("lifetime_max_footprint_mb", 0.0), 1))
    if device == "mps":
        mem["mps_current_allocated_mb"] = round(torch.mps.current_allocated_memory() / 1048576, 1)
        mem["mps_driver_allocated_mb"] = round(torch.mps.driver_allocated_memory() / 1048576, 1)
        mem["mps_driver_allocated_peak_after_batch_mb"] = round(mps_peak["driver"] / 1048576, 1)
        mem["mps_empty_cache_each_batch"] = bool(job.get("mps_empty_cache"))

    if not length_info["used"]:
        length_info["fallback"] = ("FRIGID built-in: uniform random total token length in [20, 99] (data/len.pk missing -> "
                                   "range(20,100), no length model passed)")
    length_model = length_info
    return dict(
        ok=True, name=job.get("name"), formula=job["formula"], adduct=job["adduct"], precursor_mz=job["precursor_mz"],
        mode=job["mode"], seed=seed, device_used=device, device_decision=job.get("_device_decision"),
        n_samples_generated=total_gen, n_valid=total_valid, n_formula_matched_with_duplicates=total_matched,
        n_candidates=len(cands), n_candidates_formula_match=n_match,
        true_rank=true_rank if truth_key else None, truth_evaluated=bool(truth_key),
        subformula=dict(n_peaks_in=len(job["peaks"]), n_peaks_assigned=n_assigned,
                        peaks_used_by_encoder=min(n_assigned, ENCODER_KW["max_peaks"])),
        fingerprint=dict(bits=FP_BITS, threshold=FP_THRESHOLD, n_bits_on=int(fp_bin.sum()),
                         max_prob=round(float(probs.max()), 4)),
        timings_s=dict(model_load_mist=round(models["t_mist"], 3), model_load_dlm=round(models["t_dlm"], 3),
                       model_load_total=round(models["t_mist"] + models["t_dlm"], 3),
                       subformula_and_featurize=round(t_subform, 3), encoder=round(t_enc, 3),
                       generation=round(t_gen_wall, 3), generation_sampler_only=round(gen_time, 3)),
        memory=mem, length_model=length_model,
        settings=dict(fp_threshold=FP_THRESHOLD, softmax_temp=SOFTMAX_TEMP, randomness=RANDOMNESS,
                      batch_size=BATCH_SIZE, **cfg, encoder=ENCODER_KW, mist_ckpt=os.path.basename(MIST_CKPT),
                      dlm_ckpt=os.path.basename(DLM_CKPT), torch=torch.__version__),
        warnings=out_warn, candidates=cands)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    t_start = time.perf_counter()
    result: dict
    code = 0
    guard = None

    def write(res):
        with open(args.output, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=1)

    def on_exceed(reason):  # called from the guard thread: leave a structured failure behind, then die hard
        write(dict(ok=False, error_type="memory_guard", error=reason, memory_guard=guard.snapshot() if guard else None))
        log("MEMORY GUARD TRIPPED:", reason)
        os._exit(3)

    try:
        job = validate(json.load(open(args.input, encoding="utf-8")))
        warnings.filterwarnings("ignore")
        import torch
        import memguard

        decision = memguard.choose_device(job["device"], mps_available=torch.backends.mps.is_available())
        device = decision["chosen"]
        job["_device_decision"] = decision
        log(f"job={job.get('name')} formula={job['formula']} mode={job['mode']} seed={job['seed']} | device {job['device']} -> {device}: {decision['reason']}")
        caps = {"cpu": int(os.environ.get("FRIGID_RSS_CAP_MB", memguard.DEFAULT_CAPS_MB["cpu"])),
                "mps": int(os.environ.get("FRIGID_FOOTPRINT_CAP_MB", memguard.DEFAULT_CAPS_MB["mps"]))}
        guard = memguard.MemoryGuard(device, cap_mb=caps[device], on_exceed=on_exceed).start()
        models = load_models(device)
        log(f"models loaded in {models['t_mist'] + models['t_dlm']:.1f}s")
        result = predict(job, models)
        result["timings_s"]["process_total_incl_imports"] = round(time.perf_counter() - t_start, 3)
        result["memory_guard"] = guard.stop()
        result["memory"]["peak_phys_footprint_mb"] = result["memory_guard"]["kernel_lifetime_peak_footprint_mb"]
        del models
        gc.collect()
    except InputError as e:
        result, code = dict(ok=False, error_type="invalid_input", error=str(e)), 2
    except Exception as e:  # noqa: BLE001
        result, code = dict(ok=False, error_type=type(e).__name__, error=str(e), traceback=traceback.format_exc()), 1
        log("FAILED:", type(e).__name__, e)
    write(result)
    log(f"wrote {args.output} ok={result.get('ok')}")
    return code


if __name__ == "__main__":
    sys.exit(main())
