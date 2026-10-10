#!/usr/bin/env python3
"""New numerical-repair branch; frozen upstream model sources are untouched.

Usage: --engine native16|v55 followed by unchanged original fitter arguments.
Only the nullspace linear solver is replaced with a per-call verified orthogonal
projector. The original numerical acceptance tolerance remains 1e-9. Output and
checkpoint directories must be new. Solver provenance enters every factor and
cache namespace; failed old checkpoints cannot be silently reused.
"""
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import sys
import bayesgraph_runtime as runtime
import fit_native16_bayesgraph as native
import component_nullspace_projector as projector


def main():
    p=argparse.ArgumentParser(add_help=False);p.add_argument('--engine',choices=['native16','v55'],required=True)
    a,rest=p.parse_known_args()
    meta=dict(method='verified component orthogonal nullspace projection',
        mathematical_target_unchanged='x-A.T(AA.T)^+Ax',absolute_tolerance=1e-9,
        original_upstream_sources_unchanged=True,all_original_measurement_rows_checked=True,
        component_source_sha256=hashlib.sha256(Path(projector.__file__).read_bytes()).hexdigest(),
        adapter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    original_install=runtime.install_checkpoints
    def install(mv4,directory,sources):
        directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
        projector.AUDIT_PATH=directory.parent/'component_projection_calls.jsonl'
        mv4.measurement_nullspace_center=projector.measurement_nullspace_center
        original_install(mv4,directory,{**sources,**meta})
        (directory.parent/'component_projection_method.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    original_compat=native.install_runtime_compatibility
    def compat(mv4):
        report=original_compat(mv4);report['nullspace_numerical_repair']=meta;return report
    runtime.install_checkpoints=install;native.install_runtime_compatibility=compat
    module=native if a.engine=='native16' else importlib.import_module('fit_generic_v55_bayesgraph')
    sys.argv=[module.__file__]+rest
    module.main()


if __name__=='__main__':main()
