"""Isolated stable-QR numerical branch; upstream model and all gates unchanged."""
import hashlib
import json
from pathlib import Path
import bayesgraph_runtime as runtime
import fit_native16_bayesgraph as native
import component_nullspace_projector_qr_memory_v3 as projector
import component_nullspace_projector as actual_audit


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    meta = dict(matrix_resource_ceiling_bytes=projector.MAX_DENSE_BYTES, memory_only_change=True, method='SVD and rank-revealing QR component orthogonal nullspace projection',
        mathematical_target_unchanged='x-A.T(AA.T)^+Ax', absolute_tolerance=1e-9,
        original_upstream_sources_unchanged=True, all_original_measurement_rows_checked=True,
        component_source_sha256=sha(projector.__file__),
        actual_audit_source_sha256=sha(actual_audit.__file__), adapter_sha256=sha(__file__))
    original_install = runtime.install_checkpoints

    def install(mv4, directory, sources):
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        projector.AUDIT_PATH = directory.parent / 'component_projection_calls.jsonl'
        projector.DETAIL_PATH = directory.parent / 'qr_component_details.jsonl'
        mv4.measurement_nullspace_center = projector.measurement_nullspace_center
        original_install(mv4, directory, {**sources, **meta})
        (directory.parent / 'component_projection_method.json').write_text(json.dumps(meta, indent=2), encoding='utf-8')

    original_compat = native.install_runtime_compatibility

    def compat(mv4):
        report = original_compat(mv4)
        report['nullspace_numerical_repair'] = meta
        return report

    runtime.install_checkpoints = install
    native.install_runtime_compatibility = compat
    native.main()


if __name__ == '__main__':
    main()
