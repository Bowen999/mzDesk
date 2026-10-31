#!/usr/bin/env python3
"""Reproducible B10 Simvastatin vs DMSO feature analysis and standalone HTML.

Input: 45-file MZmine area, per-cell detection state, and ISA sample mapping.
No extra imputation; Welch tests use biological (not injection) replicates.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import scipy
from scipy import stats
from scipy.special import stdtr
from plotly.offline.offline import get_plotlyjs

HERE = Path(__file__).resolve().parent
OUTPUT = Path(os.environ.get('SIMVA_B10_OUTPUT_DIR', str(HERE))).resolve()
SOURCE = Path('/Volumes/G-TECH/MetabolomicsWorkbench/MTBLS12332/RP18_pos/mzmine/result')
AREA = SOURCE / 'MTBLS12332_RP18pos_B10_alignment.csv'
STATE = SOURCE / 'MTBLS12332_RP18pos_B10_feature_state.csv'
META = SOURCE / 'MTBLS12332_RP18pos_B10_sample_map.csv'
TEMPLATE = HERE / 'simvastatin_b10_template.html'
PREFIX = 'simvastatin_b10_differential'


def sha256(path: Path) -> str:
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def number(value):
    if value is None or str(value).strip().lower() in ('', 'na', 'null', 'none'):
        return None
    try:
        n = float(value)
        return n if math.isfinite(n) else None
    except (ValueError, TypeError):
        return None


def positive(value):
    n = number(value)
    return n if n is not None and n > 0 else None


def classify_ms2(flag, spectrum, precursor, charge):
    """Apply the user-defined strong non-precursor MS/MS criterion."""
    if flag not in ('true', 'false'):
        raise ValueError(f'Unexpected MS/MS flag: {flag}')
    text = spectrum.strip()
    present = text.lower() not in ('', 'na', 'null', 'none')
    if (flag == 'true') != present:
        raise ValueError('MS/MS flag and exported spectrum disagree')
    if not present:
        return False, 0, '无导出谱', False, 0, False
    peaks = []
    for token in text.split():
        parts = token.split(':')
        if len(parts) != 2:
            raise ValueError(f'Malformed MS/MS peak: {token}')
        mz, intensity = number(parts[0]), number(parts[1])
        if mz is None or mz <= 0 or intensity is None or intensity < 0:
            raise ValueError(f'Invalid MS/MS peak: {token}')
        peaks.append((mz, intensity))
    if not peaks or precursor is None or precursor <= 0:
        raise ValueError('MS/MS spectrum or precursor is missing')
    tolerance = max(0.001, precursor * 10e-6)
    precursor_only = len(peaks) == 1 and abs(peaks[0][0] - precursor) <= tolerance
    max_intensity = max(intensity for _, intensity in peaks)
    strong_nonprecursor = 0
    for mz, intensity in peaks:
        if intensity <= 0.05 * max_intensity or abs(mz - precursor) <= tolerance:
            continue
        if charge is not None:
            isotope_positions = (precursor + 1.003354835 / charge, precursor + 2.00670967 / charge)
            if any(abs(mz - isotope) <= max(0.001, isotope * 10e-6) for isotope in isotope_positions):
                continue
        strong_nonprecursor += 1
    high_quality = strong_nonprecursor >= 3
    charge_unknown = charge is None
    if precursor_only:
        status = '仅前体峰（不计为有MS/MS）'
    elif high_quality and charge_unknown:
        status = '符合强峰数量；电荷未知，同位素待复核'
    elif high_quality:
        status = '符合强峰数量'
    elif len(peaks) == 1 and peaks[0][0] > precursor + tolerance:
        status = '单峰高于前体，需复核'
    else:
        status = f'非前体强峰不足3个（{strong_nonprecursor}个）'
    return not precursor_only, len(peaks), status, high_quality, strong_nonprecursor, charge_unknown


def evaluate_filters(treated, control, blank):
    """Presence includes DETECTED and ESTIMATED; an all-missing blank is retained."""
    t = [v for v in treated if v is not None]
    c = [v for v in control if v is not None]
    b = [v for v in blank if v is not None]
    pass_presence = max(len(t), len(c)) >= 9
    if not pass_presence:
        return False, False, False
    no_blank = not b
    pass_blank = no_blank or max(t + c) >= 3 * max(b)
    return True, pass_blank, no_blank


def independent_welch(a, b):
    va = float(np.var(a, ddof=1)) / len(a)
    vb = float(np.var(b, ddof=1)) / len(b)
    den = va + vb
    if den <= 0:
        return None
    df_den = va * va / (len(a) - 1) + vb * vb / (len(b) - 1)
    if df_den <= 0:
        return None
    t = (float(np.mean(a)) - float(np.mean(b))) / math.sqrt(den)
    df = den * den / df_den
    return 2 * float(stdtr(df, -abs(t)))


def bh_adjust(pvalues):
    order = np.argsort(pvalues, kind='stable')
    sorted_p = np.asarray(pvalues, dtype=float)[order]
    m = len(sorted_p)
    adjusted = np.minimum.accumulate((sorted_p * m / np.arange(1, m + 1))[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    out = np.empty(m, dtype=float)
    out[order] = adjusted
    return out


def selftest():
    """Focused synthetic edge cases for presence, blank policy and BH."""
    assert evaluate_filters([1.] * 9 + [None] * 3, [None] * 12, [None] * 12) == (True, True, True)
    assert evaluate_filters([3.] * 9, [None] * 12, [1.]) == (True, True, False)
    assert evaluate_filters([2.99] * 9, [None] * 12, [1.]) == (True, False, False)
    assert evaluate_filters([1.] * 8, [1.] * 8, [None]) == (False, False, False)
    assert np.allclose(bh_adjust([.01, .04, .03]), [.03, .04, .04], rtol=0, atol=1e-14)
    assert classify_ms2('true', '200.0005:120', 200.0, 1)[0] is False
    assert classify_ms2('true', '150.0000:120', 200.0, 1)[3] is False
    assert classify_ms2('true', '100:100 120:5 130:6 140:7 200:100', 200.0, 1)[4] == 3
    assert classify_ms2('true', '100:100 120:5 130:6 200:100', 200.0, 1)[3] is False
    assert classify_ms2('true', '100:100 120:20 201.00335:20 202.00671:20', 200.0, 1)[3] is False
    assert classify_ms2('true', '100:100 120:20 201.00335:20 202.00671:20', 200.0, None)[3] is True
    assert classify_ms2('false', 'null', 200.0, None)[0] is False


def load_inputs():
    with META.open(newline='', encoding='utf-8') as f:
        metas = list(csv.DictReader(f))
    by_name = {m['sample_name']: m for m in metas}
    if len(by_name) != 45:
        raise ValueError('Sample metadata must contain 45 unique names')
    groups = {g: [m['sample_name'] for m in metas if m['sample_class'] == g]
              for g in ('Simvastatin', 'DMSO', 'PBS', 'QC')}
    if {k: len(v) for k, v in groups.items()} != {'Simvastatin': 12, 'DMSO': 12, 'PBS': 12, 'QC': 9}:
        raise ValueError(f'Unexpected groups: {groups}')
    with AREA.open(newline='', encoding='utf-8') as af, STATE.open(newline='', encoding='utf-8') as sf:
        ar = csv.reader(af)
        metadata_rows = [next(ar) for _ in range(4)]
        header = next(ar)
        state_reader = csv.DictReader(sf)
        if state_reader.fieldnames is None or state_reader.fieldnames[0] != 'Alignment ID':
            raise ValueError('Unexpected state-file header')
        sample_order = state_reader.fieldnames[1:]
        if set(sample_order) != set(by_name) or len(sample_order) != 45:
            raise ValueError('Sample columns in state table do not match metadata')
        if header[35:80] != sample_order:
            raise ValueError('Sample order in area and state tables differs')
        for i, (label, field) in enumerate((('Class', 'sample_class'), ('File type', 'file_type'),
                                             ('Injection order', 'acquisition_id'), ('Batch ID', 'batch'))):
            if metadata_rows[i][34] != label:
                raise ValueError(f'Missing sample metadata row: {label}')
            if metadata_rows[i][35:80] != [by_name[name][field] for name in sample_order]:
                raise ValueError(f'Sample metadata mismatch: {label}')
        col = {name: i for i, name in enumerate(header)}
        for field in ('Alignment ID', 'Average Mz', 'Average Rt(min)', 'Metabolite name', 'Adduct type',
                      'MS/MS assigned', 'MS/MS spectrum',
                      'Formula', 'CCS (A2)', 'Average 1/K0 (Vs/cm2)', 'Weighted dot product'):
            if field not in col:
                raise ValueError(f'Missing area-table column: {field}')
        filtered = []
        count = Counter()
        ids = set()
        for row, state_row in zip(ar, state_reader, strict=True):
            count['input_features'] += 1
            fid = row[col['Alignment ID']]
            if fid in ids or state_row['Alignment ID'] != fid:
                raise ValueError(f'Duplicate or misaligned feature ID: {fid}')
            ids.add(fid)
            if len(row) != len(header):
                raise ValueError(f'Invalid area row width: {fid}')
            ms2_flag = row[col['MS/MS assigned']].strip().lower()
            charge_value = positive(row[col['Charge']])
            if charge_value is not None and (not charge_value.is_integer() or charge_value > 10):
                raise ValueError(f'Unexpected charge: {fid}, {charge_value}')
            charge = int(charge_value) if charge_value is not None else None
            has_ms2, ms2_peak_count, ms2_status, has_hq_ms2, ms2_strong_peaks, ms2_charge_unknown = classify_ms2(
                ms2_flag, row[col['MS/MS spectrum']], positive(row[col['Average Mz']]), charge)
            if ms2_flag == 'true':
                count['input_features_with_exported_ms2'] += 1
            if ms2_status == '仅前体峰（不计为有MS/MS）':
                count['input_features_precursor_only_ms2'] += 1
            if has_ms2:
                count['input_features_with_ms2'] += 1
            if has_hq_ms2:
                count['input_features_with_hq_ms2'] += 1
            vals = {}
            for name in sample_order:
                value = positive(row[col[name]])
                state = state_row[name].upper()
                if state not in ('DETECTED', 'ESTIMATED', 'UNKNOWN'):
                    raise ValueError(f'Unexpected feature state: {fid}, {name}, {state}')
                if (state == 'UNKNOWN') != (value is None):
                    raise ValueError(f'Value/state mismatch: {fid}, {name}')
                vals[name] = value
            present, blankpass, blank_missing = evaluate_filters(
                [vals[n] for n in groups['Simvastatin']], [vals[n] for n in groups['DMSO']],
                [vals[n] for n in groups['PBS']])
            if not present:
                count['failed_presence'] += 1
                continue
            count['passed_presence'] += 1
            if not blankpass:
                count['failed_blank'] += 1
                continue
            count['passed_blank'] += 1
            b = [vals[n] for n in groups['PBS'] if vals[n] is not None]
            bs = [state_row[n] for n in groups['PBS'] if vals[n] is not None]
            if blank_missing:
                blank_status = '空白全部缺失'
                count['blank_all_missing'] += 1
            elif 'DETECTED' in bs:
                blank_status = '空白有直接检出'
                count['blank_direct_support'] += 1
            else:
                blank_status = '空白仅补峰'
                count['blank_estimated_only'] += 1
            name = row[col['Metabolite name']].strip()
            annotated = name.lower() not in ('', 'unknown', 'na', 'n/a', 'null', 'none')
            if annotated and not has_ms2:
                count['annotated_without_qualifying_representative_ms2'] += 1
            count['annotated_after_filters' if annotated else 'unknown_after_filters'] += 1
            if has_ms2:
                count['ms2_after_filters'] += 1
                if not annotated:
                    count['unannotated_ms2_after_filters'] += 1
            if has_hq_ms2:
                count['hq_ms2_after_filters'] += 1
                if not annotated:
                    count['unannotated_hq_ms2_after_filters'] += 1
            qc_direct = [vals[n] for n in groups['QC'] if state_row[n] == 'DETECTED' and vals[n] is not None]
            qc_cv = (float(np.std(qc_direct, ddof=1) / np.mean(qc_direct))
                     if len(qc_direct) >= 3 and np.mean(qc_direct) > 0 else None)
            filtered.append({
                'id': fid, 'mz': positive(row[col['Average Mz']]), 'rt': number(row[col['Average Rt(min)']]),
                'name': name if annotated else '未注释', 'annotated': annotated, 'has_ms2': has_ms2,
                'ms2_peak_count': ms2_peak_count, 'ms2_status': ms2_status,
                'has_hq_ms2': has_hq_ms2, 'ms2_strong_peaks': ms2_strong_peaks,
                'ms2_charge_unknown': ms2_charge_unknown,
                'is_treatment_compound': annotated and name.casefold() == 'simvastatin',
                'adduct': row[col['Adduct type']], 'formula': row[col['Formula']],
                'mobility': number(row[col['Average 1/K0 (Vs/cm2)']]), 'ccs': number(row[col['CCS (A2)']]),
                'library_score': number(row[col['Weighted dot product']]), 'blank_status': blank_status,
                'blank_max': max(b) if b else None, 'qc_cv_direct': qc_cv, 'qc_direct_n': len(qc_direct),
                'treated_present': sum(vals[n] is not None for n in groups['Simvastatin']),
                'control_present': sum(vals[n] is not None for n in groups['DMSO']),
                'treated_direct': sum(state_row[n] == 'DETECTED' for n in groups['Simvastatin']),
                'control_direct': sum(state_row[n] == 'DETECTED' for n in groups['DMSO']),
                'treated_estimated': sum(state_row[n] == 'ESTIMATED' for n in groups['Simvastatin']),
                'control_estimated': sum(state_row[n] == 'ESTIMATED' for n in groups['DMSO']),
                'values': vals,
            })
    if count['input_features'] != count['failed_presence'] + count['passed_presence']:
        raise AssertionError('Presence accounting mismatch')
    if count['passed_presence'] != count['failed_blank'] + count['passed_blank']:
        raise AssertionError('Blank accounting mismatch')
    return filtered, groups, by_name, count


def calculate(filtered, groups, metas, count):
    samples = groups['Simvastatin'] + groups['DMSO']
    matrix = np.array([[row['values'][n] if row['values'][n] is not None else np.nan for n in samples]
                       for row in filtered], dtype=float)
    sample_medians = np.nanmedian(matrix, axis=0)
    if not np.all(np.isfinite(sample_medians) & (sample_medians > 0)):
        raise ValueError('A sample has no usable positive median')
    reference_median = float(np.median(sample_medians))
    norm = matrix / sample_medians * reference_median
    indexes = {n: i for i, n in enumerate(samples)}
    rep_samples = {}
    for group in ('Simvastatin', 'DMSO'):
        grouped = defaultdict(list)
        for name in groups[group]:
            grouped[metas[name]['biological_replicate']].append(indexes[name])
        if len(grouped) != 4 or any(len(v) != 3 for v in grouped.values()):
            raise ValueError(f'Unexpected biological/technical structure: {group}')
        rep_samples[group] = [grouped[k] for k in sorted(grouped, key=int)]
    pvalues = []
    tested = []
    manual_checks = 0
    for i, row in enumerate(filtered):
        biological = {}
        tech_counts = {}
        for group in ('Simvastatin', 'DMSO'):
            biological[group] = []
            tech_counts[group] = []
            for ids in rep_samples[group]:
                valid = norm[i, ids]
                valid = valid[np.isfinite(valid) & (valid > 0)]
                biological[group].append(float(np.mean(valid)) if len(valid) else None)
                tech_counts[group].append(int(len(valid)))
        row['replicates_treated'] = biological['Simvastatin']
        row['replicates_control'] = biological['DMSO']
        row['technical_counts_treated'] = tech_counts['Simvastatin']
        row['technical_counts_control'] = tech_counts['DMSO']
        a = np.log2([v for v in biological['Simvastatin'] if v is not None])
        b = np.log2([v for v in biological['DMSO'] if v is not None])
        row['n_treated'] = len(a)
        row['n_control'] = len(b)
        row['log2fc'] = None
        row['fc'] = None
        row['p'] = None
        row['q'] = None
        if len(a) < 3 or len(b) < 3:
            row['test_status'] = '生物学重复不足（每组至少3个）'
            count['excluded_insufficient_bio'] += 1
            continue
        lfc = float(np.mean(a) - np.mean(b))
        check_p = independent_welch(a, b)
        if check_p is None:
            row['test_status'] = '两组方差均为零，无法可靠检验'
            count['excluded_zero_variance'] += 1
            continue
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            result = stats.ttest_ind(a, b, equal_var=False)
        p = float(result.pvalue)
        if not math.isfinite(p) or p < 0 or p > 1:
            row['test_status'] = 'Welch 检验数值无效'
            count['excluded_nonfinite_p'] += 1
            continue
        if not math.isclose(p, check_p, rel_tol=1e-9, abs_tol=1e-14):
            raise AssertionError(f'Independent Welch check failed: {row["id"]}')
        manual_checks += 1
        row.update(log2fc=lfc, fc=math.exp2(lfc), p=p, test_status='已检验')
        tested.append(row)
        pvalues.append(p)
    if not tested:
        raise ValueError('No statistically testable features')
    adjusted = bh_adjust(pvalues)
    ordered = sorted(pvalues)
    for j, original_p in enumerate(ordered):
        brute = min(1., min(ordered[k] * len(ordered) / (k + 1) for k in range(j, len(ordered))))
        actual = float(np.sort(adjusted)[j])
        if not math.isclose(brute, actual, rel_tol=1e-10, abs_tol=1e-14):
            raise AssertionError(f'Independent BH check failed at rank {j}')
        if actual + 1e-14 < original_p:
            raise AssertionError('Adjusted q is smaller than raw p')
    for row, q in zip(tested, adjusted):
        row['q'] = float(q)
    count['tested'] = len(tested)
    count['welch_independent_checks'] = manual_checks
    count['bh_independent_checks'] = len(tested)
    if count['passed_blank'] != (count['tested'] + count['excluded_insufficient_bio'] +
                                  count['excluded_zero_variance'] + count['excluded_nonfinite_p']):
        raise AssertionError('Test eligibility accounting mismatch')
    return tested, {name: float(sample_medians[i]) for i, name in enumerate(samples)}, reference_median


def write_csv(path, records):
    fields = ['id', 'mz', 'rt', 'name', 'annotated', 'has_ms2', 'has_hq_ms2', 'ms2_strong_peaks', 'ms2_charge_unknown', 'ms2_peak_count', 'ms2_status', 'is_treatment_compound', 'adduct', 'formula', 'mobility', 'ccs', 'library_score',
              'blank_status', 'blank_max', 'qc_cv_direct', 'qc_direct_n', 'treated_present', 'control_present',
              'treated_direct', 'control_direct', 'treated_estimated', 'control_estimated', 'n_treated',
              'n_control', 'log2fc', 'fc', 'p', 'q', 'test_status', 'replicates_treated', 'replicates_control',
              'technical_counts_treated', 'technical_counts_control']
    with path.open('x', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        for row in records:
            obj = {k: (json.dumps(row[k], ensure_ascii=False) if isinstance(row.get(k), list) else row.get(k))
                   for k in fields}
            w.writerow(obj)


def main():
    selftest()
    paths = {name: OUTPUT / f'{PREFIX}{suffix}' for name, suffix in
             [('html', '.html'), ('csv', '.csv'), ('audit', '_audit.json')]}
    collisions = [str(p) for p in paths.values() if p.exists()]
    if collisions:
        raise FileExistsError(f'Refusing to overwrite existing outputs: {collisions}')
    filtered, groups, metas, count = load_inputs()
    tested, medians, refmedian = calculate(filtered, groups, metas, count)
    display_keys = ('id', 'mz', 'rt', 'name', 'annotated', 'has_ms2', 'has_hq_ms2', 'ms2_strong_peaks', 'ms2_charge_unknown', 'ms2_peak_count', 'ms2_status', 'is_treatment_compound', 'adduct', 'formula', 'mobility', 'ccs',
                    'library_score', 'blank_status', 'blank_max', 'qc_cv_direct', 'qc_direct_n',
                    'treated_present', 'control_present', 'treated_direct', 'control_direct',
                    'treated_estimated', 'control_estimated', 'n_treated', 'n_control', 'log2fc',
                    'fc', 'p', 'q', 'replicates_treated', 'replicates_control',
                    'technical_counts_treated', 'technical_counts_control')
    points = [{key: row[key] for key in display_keys} for row in tested]
    initial = [r for r in tested if r['q'] <= .05 and abs(r['log2fc']) >= math.log2(1.5)]
    initial_counts = {'total': len(initial), 'up': sum(r['log2fc'] > 0 for r in initial),
                      'down': sum(r['log2fc'] < 0 for r in initial),
                      'annotated': sum(r['annotated'] for r in initial),
                      'unannotated': sum(not r['annotated'] for r in initial),
                      'unannotated_with_ms2': sum(not r['annotated'] and r['has_ms2'] for r in initial),
                      'features_with_ms2': sum(r['has_ms2'] for r in initial),
                      'unannotated_with_hq_ms2': sum(not r['annotated'] and r['has_hq_ms2'] for r in initial),
                      'unannotated_hq_ms2_charge_unknown': sum(not r['annotated'] and r['has_hq_ms2'] and r['ms2_charge_unknown'] for r in initial),
                      'features_with_hq_ms2': sum(r['has_hq_ms2'] for r in initial),
                      'blank_missing': sum(r['blank_status']=='空白全部缺失' for r in initial),
                      'blank_estimated_only': sum(r['blank_status']=='空白仅补峰' for r in initial),
                      'treatment_compound': sum(r['is_treatment_compound'] for r in initial)}
    if initial_counts['up'] + initial_counts['down'] != initial_counts['total']:
        raise AssertionError('Initial direction counts do not sum to total')
    audit = {
        'source_sha256': {str(p): sha256(p) for p in (AREA, STATE, META)},
        'parameters': {'presence_per_group': 'at least 9/12 injections in either group; detected and MZmine-estimated count',
                       'background': 'max(treated and DMSO) >= 3 * max(PBS); all-missing PBS retained and flagged',
                       'normalization': 'sample positive median after filters, scaled to median of sample medians',
                       'technical_replicates': 'arithmetic mean of available normalized injections for each biological replicate',
                       'extra_imputation': False, 'min_biological_replicates_per_group': 3,
                       'test': 'two-sided Welch t-test of log2 biological means',
                       'fold_change': '2^(mean(log2 treated biological means) - mean(log2 control biological means))',
                       'multiple_testing': 'Benjamini-Hochberg across all finite test p-values',
                       'ms2_definition': 'High-quality criterion: >=3 peaks strictly above 5% of the maximum intensity across the complete representative spectrum, excluding the precursor and, when charge is known, precursor M+1 and M+2 within max(0.001 Da, 10 ppm). Unknown charge is not inferred and remains provisional.',
                       'qc_filter': False, 'default_fc': 1.5, 'default_q': .05},
        'counts': dict(count), 'sample_medians': medians, 'reference_median': refmedian,
        'default_threshold_counts': initial_counts,
        'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__},
    }
    display = {'points': points, 'counts': dict(count), 'initial': initial_counts, 'parameters': audit['parameters'],
               'sample_medians': medians, 'source_sha256': audit['source_sha256']}
    content = TEMPLATE.read_text(encoding='utf-8')
    if content.count('%%PLOTLY_JS%%') != 1 or content.count('%%DATA_JSON%%') != 1:
        raise ValueError('Invalid HTML template placeholders')
    plotly = get_plotlyjs().replace('</script', '<\\/script')
    data = json.dumps(display, ensure_ascii=False, allow_nan=False, separators=(',', ':')).replace('<', '\\u003c')
    content = content.replace('%%PLOTLY_JS%%', plotly).replace('%%DATA_JSON%%', data)
    write_csv(paths['csv'], filtered)
    with paths['audit'].open('x', encoding='utf-8') as f:
        f.write(json.dumps(audit, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    with paths['html'].open('x', encoding='utf-8') as f:
        f.write(content)
    print(json.dumps({'outputs': {k: str(v) for k, v in paths.items()}, 'counts': dict(count),
                      'default_threshold_counts': initial_counts}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
