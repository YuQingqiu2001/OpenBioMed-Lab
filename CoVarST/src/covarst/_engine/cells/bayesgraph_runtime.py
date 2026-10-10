#!/usr/bin/env python3
"""Process-local exact-result checkpoints and timing for unchanged solvers.

Keys hash original source, runtime versions, and every numerical argument. No
approximation or alternative solver is used. HDF5 avoids executable pickle data.
"""
import dataclasses
from functools import wraps
import hashlib
import json
import os
from pathlib import Path
import time
import h5py
import numpy as np
import scipy
from scipy import sparse


def update_digest(digest, value):
    if sparse.issparse(value):
        digest.update(str((type(value).__name__, value.shape)).encode())
        for array in (value.data, value.indices, value.indptr):
            update_digest(digest, array)
    elif isinstance(value, np.ndarray):
        array = np.ascontiguousarray(value)
        digest.update(str((array.dtype.str, array.shape)).encode())
        digest.update(memoryview(array).cast('B'))
    elif isinstance(value, (list, tuple)):
        digest.update(type(value).__name__.encode())
        for item in value:
            update_digest(digest, item)
    elif isinstance(value, dict):
        for key in sorted(value):
            digest.update(str(key).encode())
            update_digest(digest, value[key])
    else:
        digest.update(repr(value).encode())


def save_value(group, value):
    if dataclasses.is_dataclass(value):
        group.attrs['kind'] = 'weak_state'
        save_value(group.create_group('value'), vars(value))
    elif isinstance(value, np.ndarray):
        group.attrs['kind'] = 'array'
        group.create_dataset('value', data=value, **({'compression': 'lzf'} if value.ndim else {}))
    elif isinstance(value, dict):
        group.attrs['kind'] = 'dict'
        for key, item in value.items():
            save_value(group.create_group(str(key)), item)
    elif isinstance(value, (list, tuple)):
        group.attrs['kind'] = type(value).__name__
        for index, item in enumerate(value):
            save_value(group.create_group(str(index)), item)
    else:
        group.attrs['kind'] = 'scalar'
        group.attrs['value'] = json.dumps(value.item() if isinstance(value, np.generic) else value)


def load_value(group, mv4):
    kind = group.attrs['kind']
    if kind == 'array':
        return group['value'][()]
    if kind == 'dict':
        return {key: load_value(group[key], mv4) for key in group}
    if kind in ('list', 'tuple'):
        result = [load_value(group[str(i)], mv4) for i in range(len(group))]
        return tuple(result) if kind == 'tuple' else result
    if kind == 'scalar':
        return json.loads(group.attrs['value'])
    if kind == 'weak_state':
        return mv4.WeakStateResult(**load_value(group['value'], mv4))
    raise ValueError('Unknown checkpoint value type')


def install_checkpoints(mv4, directory, source_hashes):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    source_key = json.dumps({'sources': source_hashes, 'numpy': np.__version__,
                             'scipy': scipy.__version__}, sort_keys=True).encode()
    def install(owner, name):
        original = getattr(owner, name)
        @wraps(original)
        def checked(*args, **kwargs):
            digest = hashlib.sha256(source_key + name.encode())
            update_digest(digest, args)
            update_digest(digest, kwargs)
            key = digest.hexdigest()
            target = directory / f'{name}.{key}.h5'
            if target.is_file():
                with h5py.File(target, 'r') as h:
                    if h.attrs.get('complete') != 1 or h.attrs.get('key') != key:
                        raise ValueError(f'Invalid numerical checkpoint: {target}')
                    result = load_value(h['result'], mv4)
                print(f'[checkpoint] reused {name} {key[:12]}', flush=True)
                return result
            started = time.monotonic()
            detail = ''
            if name == 'lsmr':
                matrix = args[0] if args else kwargs['A']
                detail = f' shape={matrix.shape} maxiter={kwargs.get("maxiter")} atol={kwargs.get("atol")}'
            print(f'[numerics] start {name}{detail} key={key[:12]}', flush=True)
            result = original(*args, **kwargs)
            temporary = target.with_name(target.name + '.partial')
            with h5py.File(temporary, 'x') as h:
                h.attrs.update(complete=0, key=key, function=name,
                               elapsed_seconds=time.monotonic() - started)
                save_value(h.create_group('result'), result)
                h.attrs['complete'] = 1
            os.replace(temporary, target)
            detail = f' istop={result[1]} iterations={result[2]}' if name == 'lsmr' else ''
            print(f'[numerics] complete {name} seconds={time.monotonic()-started:.3f}{detail}', flush=True)
            return result
        setattr(owner, name, checked)
    for name in ['fit_slice_internal_weak_state', 'solve_graph_inverse', 'orthogonalize_nullspace_contrasts']:
        install(mv4, name)
    install(mv4.splinalg, 'lsmr')
