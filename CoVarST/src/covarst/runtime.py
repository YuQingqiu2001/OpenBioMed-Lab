"""Resolve frozen mathematical modules; do not import unrelated legacy CLIs."""
from pathlib import Path
import sys, hashlib, json

PACKAGE = Path(__file__).resolve().parent
ENGINE = PACKAGE / '_engine'
PROJECT = PACKAGE.parents[1]

def activate(*names):
    for name in names:
        path = ENGINE / name
        if str(path) not in sys.path:
            sys.path.insert(0, str(path))

def asset_root(value=None):
    root = Path(value) if value else PROJECT
    if not (root / 'references').is_dir():
        raise FileNotFoundError('Use --assets to point to the CoVarST checkout containing models/ and references/')
    return root

def sha256(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8<<20),b''):digest.update(block)
    return digest.hexdigest()

def fresh(path):
    path=Path(path)
    if path.exists():raise FileExistsError(f'Refusing to overwrite {path}')
    path.parent.mkdir(parents=True,exist_ok=True)
    return path

def write_json(path,payload):
    path=fresh(path)
    path.write_text(json.dumps(payload,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
