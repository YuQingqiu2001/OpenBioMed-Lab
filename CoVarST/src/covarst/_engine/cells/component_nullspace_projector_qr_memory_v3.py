"""Stable component row-space projector without squaring the condition number.

Same x-A.T(AA.T)^+Ax target and unchanged actual-call gates from v1. Small
components use SVD; larger bounded components use rank-revealing pivoted QR of
A.T. Numerical rank cutoff is machine precision, not a biological parameter.
No expression reference or model predictions are used to determine the basis.
"""
from collections import OrderedDict
import json
from pathlib import Path
import numpy as np
from scipy import linalg, sparse
from scipy.sparse.csgraph import connected_components
import component_nullspace_projector as audited

AUDIT_PATH = None
DETAIL_PATH = None
CACHE = OrderedDict()
# Resource ceiling only; numerical rank and acceptance tolerances are unchanged.
MAX_DENSE_BYTES = 40 * 2**30


class Projector:
    def __init__(self, a):
        self.parts = []
        self.counts = {'dense_SVD': 0, 'rank_revealing_QR_AT': 0}
        _, labels = connected_components((a @ a.T).tocsr(), directed=False)
        for label in np.unique(labels):
            rows = np.flatnonzero(labels == label)
            sub = a[rows]
            columns = np.unique(sub.indices)
            sub = sub[:, columns]
            if not len(columns):
                continue
            required = sub.shape[0] * sub.shape[1] * 8
            if required > MAX_DENSE_BYTES:
                raise MemoryError(f'Component {sub.shape} exceeds verified dense QR budget')
            available = int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines()
                                 if line.startswith('MemAvailable:'))) * 1024
            estimated_peak = 8 * required + 4 * 2**30
            if available < estimated_peak + 32 * 2**30:
                raise MemoryError('Insufficient available RAM for bounded QR workspace and 32GiB reserve')
            print('QR memory guard: ' + json.dumps(dict(component=int(label), matrix_bytes=required,
                  estimated_peak_bytes=estimated_peak, available_bytes=available,
                  ceiling_bytes=MAX_DENSE_BYTES)), flush=True)
            if len(rows) <= 256:
                _, singular, vh = linalg.svd(sub.toarray(), full_matrices=False, check_finite=False)
                threshold = np.finfo(float).eps * max(sub.shape) * singular.max()
                basis = vh[singular > threshold].T
                rank = basis.shape[1]
                kind = 'dense_SVD'
            else:
                q, r, permutation = linalg.qr(sub.toarray().T, mode='economic',
                    pivoting=True, check_finite=False, overwrite_a=True)
                diagonal = np.abs(np.diag(r))
                threshold = np.finfo(float).eps * max(sub.shape) * diagonal.max()
                rank = int(np.count_nonzero(diagonal > threshold))
                # Pivoted-QR rank is contiguous; reject anomalous non-monotone cutoff.
                if np.any(diagonal[rank:] > threshold) or np.any(diagonal[:rank] <= threshold):
                    raise ArithmeticError('Non-contiguous numerical QR rank')
                basis = np.ascontiguousarray(q[:, :rank])
                del q, r, permutation
                kind = 'rank_revealing_QR_AT'
            self.parts.append((columns, basis))
            self.counts[kind] += 1
            if DETAIL_PATH is not None and kind != 'dense_SVD':
                row = dict(component=int(label), shape=list(sub.shape), numerical_rank=rank,
                           machine_precision_cutoff=float(threshold), method=kind,
                           matrix_sha256=audited.matrix_key(sub))
                with Path(DETAIL_PATH).open('a', encoding='utf-8') as stream:
                    stream.write(json.dumps(row) + '\n')
                print('Stable QR component: ' + json.dumps(row), flush=True)

    def apply(self, x):
        z = x.copy()
        for columns, q in self.parts:
            val = x[columns]
            z[columns] = val - q @ (q.T @ val)
        return z


def measurement_nullspace_center(measurement, values, tolerance=1e-9):
    # Reuse exactly the previous independent per-call checks, not a weaker gate.
    audited.Projector = Projector
    audited.CACHE = CACHE
    audited.AUDIT_PATH = AUDIT_PATH
    return audited.measurement_nullspace_center(measurement, values, tolerance)


def self_test():
    rng = np.random.default_rng(42)
    for m, n, rank in ((5, 11, 5), (8, 12, 4), (300, 420, 280)):
        left, _ = np.linalg.qr(rng.normal(size=(m, rank)))
        right, _ = np.linalg.qr(rng.normal(size=(n, rank)))
        dense = (left * np.geomspace(1, 1e-6, rank)) @ right.T
        x = rng.normal(size=(n, 3))
        expected = x - right @ (right.T @ x)
        observed, error = measurement_nullspace_center(sparse.csr_matrix(dense), x)
        np.testing.assert_allclose(observed, expected, atol=1e-8, rtol=1e-8)
        assert error <= 1e-9
    print('PASS: independent known-row-space oracle including connected ill-conditioned rank deficiency', flush=True)


if __name__ == '__main__':
    self_test()
