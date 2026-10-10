#!/usr/bin/env python3
"""Strict, read-only acceptance of factors plus original solver diagnostics."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
import fit_native16_bayesgraph as native


def validate(directory, sample, modality='16um', tolerance=1e-9):
    reasons, factors, solvers = [], [], []
    for rank in (0, 8):
        path = directory / f'{sample}.{modality}.BayesGraph_rank{rank}_all_gene.factorized.h5'
        if not path.is_file():
            reasons.append(f'missing rank{rank} factor')
            continue
        with h5py.File(path, 'r') as h:
            if int(h.attrs.get('complete', 0)) != 1:
                reasons.append(f'incomplete rank{rank} factor')
                continue
            summary = json.loads(h.attrs['summary_json'])
            finite = all(np.isfinite(h['model/' + key][:]).all() for key in
                         ['class_gene_probability', 'program_loadings', 'cell_program_activity'])
            if not finite:
                reasons.append(f'nonfinite rank{rank} model')
            class_rows = []
            for row in summary.get('class_refinement', []):
                errors = {
                    'weak_nullspace': row.get('weak_model', {}).get('nullspace_centering_max_abs_error', 0.),
                    'inverse_nullspace': row.get('inverse', {}).get('contrast_centering_max_abs_error', 0.),
                    'orthogonal_nullspace': row.get('orthogonalization', {}).get('centering_max_abs_error', 0.)}
                statuses = row.get('inverse', {}).get('cg_status', [])
                if any(not np.isfinite(v) or v > tolerance for v in errors.values()):
                    reasons.append(f'rank{rank} {row["class"]}: final centering error exceeds {tolerance}')
                if any(v != 0 for v in statuses):
                    reasons.append(f'rank{rank} {row["class"]}: CG nonconvergence')
                class_rows.append({'class': row['class'], 'gate': row.get('gate'),
                                   'centering_max_abs_error': errors, 'cg_status': statuses})
            factors.append({'rank': rank, 'path': str(path), 'finite': finite,
                            'program_count_by_class': h['model/program_count_by_class'][:].tolist(),
                            'class_refinement': class_rows})
    checkpoint = directory / 'numerical_checkpoints'
    for path in sorted(checkpoint.glob('lsmr.*.h5')):
        with h5py.File(path, 'r') as h:
            if int(h.attrs.get('complete', 0)) != 1:
                reasons.append(f'incomplete numerical checkpoint {path.name}')
                continue
            row = {key: json.loads(h[f'result/{index}'].attrs['value']) for key, index in
                   [('istop', 1), ('iterations', 2), ('norm_residual', 3), ('normal_residual', 4)]}
            row.update(file=str(path), seconds=float(h.attrs['elapsed_seconds']))
            solvers.append(row)
    bad_lsmr = [r for r in solvers if r['istop'] not in [0, 1, 2, 4, 5]]
    if bad_lsmr:
        reasons.append(f'{len(bad_lsmr)} unique LSMR checkpoints reported condition/iteration-limit nonconvergence')
    status = 'validated' if not reasons and len(factors) == 2 else 'numerical_validation_failed'
    return {'sample': sample, 'modality': modality, 'status': status,
            'nullspace_absolute_tolerance': tolerance, 'tolerance_not_relaxed': True,
            'acceptance_requires_converged_LSMR_and_CG': True,
            'reasons': reasons, 'factors': factors,
            'unique_lsmr_checkpoints': solvers, 'unique_lsmr_nonconverged_count': len(bad_lsmr),
            'counts_are_unique_cached_solve_results_not_invocation_counts': True}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, required=True)
    p.add_argument('--sample', required=True)
    p.add_argument('--modality', default='16um')
    p.add_argument('--output-json', type=Path, required=True)
    a = p.parse_args()
    report = validate(a.directory, a.sample, a.modality)
    native.write_json(a.output_json, report)
    print(json.dumps({k: v for k, v in report.items() if k not in ['factors', 'unique_lsmr_checkpoints']}, indent=2))
    raise SystemExit(0 if report['status'] == 'validated' else 2)


if __name__ == '__main__':
    main()
