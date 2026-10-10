#!/usr/bin/env python3
"""Numerically verified orthogonal projection onto ker(A), using existing SciPy.

The objective is unchanged: P(x)=x-A.T(AA.T)^+Ax. Decompose disconnected
measurement components; use small dense SVD or structurally matched independent
rows with sparse LU. Row selection is ONLY accepted when ALL original rows
annihilate the result and orthogonality/idempotence checks pass. No ridge,
tolerance relaxation, training labels, expression reference or regularization.
Every actual call is audited; any failed component aborts the new branch.
"""
from collections import OrderedDict
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy import sparse,linalg
from scipy.sparse.csgraph import connected_components,maximum_bipartite_matching
from scipy.sparse.linalg import splu

CACHE=OrderedDict()
AUDIT_PATH=None


def matrix_key(a):
    digest=hashlib.sha256(str(a.shape).encode())
    for value in (a.data,a.indices,a.indptr):digest.update(np.ascontiguousarray(value).tobytes())
    return digest.hexdigest()


class Projector:
    def __init__(self,a):
        self.a=a;self.parts=[];counts={'dense_SVD':0,'structural_subset_Gram_LU':0}
        _,labels=connected_components((a@a.T).tocsr(),directed=False)
        for label in np.unique(labels):
            rows=np.flatnonzero(labels==label);sub=a[rows];columns=np.unique(sub.indices);sub=sub[:,columns]
            if not len(columns):continue
            if len(rows)<=256:
                _,values,vh=linalg.svd(sub.toarray(),full_matrices=False,check_finite=False)
                cut=np.finfo(float).eps*max(sub.shape)*values.max()
                basis=vh[values>cut]
                solver=('dense_SVD',basis)
            else:
                matched=maximum_bipartite_matching(sub,perm_type='column')
                selected=sub[matched>=0]
                factor=splu((selected@selected.T).tocsc())
                solver=('structural_subset_Gram_LU',(selected,factor))
            counts[solver[0]]+=1;self.parts.append((columns,sub,solver))
        self.counts=counts

    def apply(self,x):
        z=x.copy()
        for columns,sub,(kind,solver) in self.parts:
            val=x[columns]
            if kind=='dense_SVD':out=val-solver.T@(solver@val)
            else:
                b,lu=solver
                out=val-b.T@lu.solve(b@val)
                for _ in range(6):
                    residual=b@out
                    if np.max(np.abs(residual),initial=0)<1e-12:break
                    out-=b.T@lu.solve(residual)
            z[columns]=out
        return z


def measurement_nullspace_center(measurement,values,tolerance=1e-9):
    start=time.monotonic();a=sparse.csr_matrix(measurement,dtype=np.float64)
    x=np.asarray(values,dtype=np.float64)
    if x.ndim==1:x=x[:,None]
    if x.ndim!=2 or len(x)!=a.shape[1] or not np.isfinite(x).all():raise ValueError('Projection input invalid')
    if tolerance>1e-9:raise ValueError('Cannot relax original 1e-9 absolute gate')
    key=matrix_key(a)
    if key not in CACHE:
        CACHE[key]=Projector(a)
        while len(CACHE)>4:CACHE.popitem(last=False)
    CACHE.move_to_end(key);projector=CACHE[key]
    z=projector.apply(x)
    if not np.isfinite(z).all():raise ValueError('Nonfinite projection')
    error=float(np.max(np.abs(a@z),initial=0))
    repeated=projector.apply(z)
    idempotence=float(np.max(np.abs(repeated-z),initial=0))
    dot=np.sum(z*(x-z),axis=0)
    normalized_orthogonality=float(np.max(np.abs(dot)/np.maximum(np.sum(x*x,axis=0),1),initial=0))
    norm_excess=float(np.max(np.linalg.norm(z,axis=0)-np.linalg.norm(x,axis=0),initial=0))
    passed=bool(error<=tolerance and idempotence<=1e-8*max(1,float(np.max(np.abs(x),initial=0)))
        and normalized_orthogonality<=1e-9 and norm_excess<=1e-8*max(1,float(np.linalg.norm(x))))
    audit=dict(schema='verified_component_nullspace_projection_v1',matrix_sha256=key,
        shape=list(a.shape),n_rhs=x.shape[1],original_absolute_tolerance=tolerance,
        max_abs_A_times_centered=error,max_abs_idempotence_error=idempotence,
        normalized_null_removed_inner_product=normalized_orthogonality,norm_excess=norm_excess,
        solver_component_counts=projector.counts,passed=passed,elapsed_seconds=time.monotonic()-start,
        input_sha256=hashlib.sha256(x.tobytes()).hexdigest(),output_sha256=hashlib.sha256(z.tobytes()).hexdigest())
    if AUDIT_PATH is not None:
        with Path(AUDIT_PATH).open('a',encoding='utf-8') as stream:stream.write(json.dumps(audit)+'\n')
    if not passed:raise ArithmeticError('Orthogonal projection failed original constraints: '+json.dumps(audit))
    return z,error


def self_test():
    rng=np.random.default_rng(42)
    for m,n,r in [(5,11,5),(7,12,4),(300,600,290)]:
        # Last fixture exercises structural singularity, not a dense generic matrix.
        if m>256:
            dense=np.zeros((m,n));dense[np.arange(r),np.arange(r)]=1
            dense[:r,r:2*r]=rng.uniform(.2,.8,r)*np.eye(r)
            dense[r:]=dense[:m-r]
        else:dense=rng.normal(size=(m,r))@rng.normal(size=(r,n))
        a=sparse.csr_matrix(dense);x=rng.normal(size=(n,3))
        _,sv,vh=linalg.svd(dense,full_matrices=False)
        q=vh[sv>np.finfo(float).eps*max(dense.shape)*sv.max()]
        exact=x-q.T@(q@x);actual,error=measurement_nullspace_center(a,x)
        np.testing.assert_allclose(actual,exact,atol=1e-9,rtol=1e-9)
        assert error<=1e-9
    print('PASS: independent dense-SVD oracle; full-rank, numerical rank-deficient and structurally rank-deficient fixtures.')


if __name__=='__main__':self_test()
