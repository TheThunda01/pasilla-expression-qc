"""Verify the source files shared by the analysis and result audit."""
import hashlib
import json
from pathlib import Path

COUNTS = "pasilla_gene_counts.tsv"
ANNOTATIONS = "pasilla_sample_annotation.csv"


def fingerprint(path):
    path = Path(path)
    if not path.is_file() or path.is_symlink():
        raise ValueError(f'Expected a regular non-symlink file: {path.name}')
    digest = hashlib.sha256()
    size = 0
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
            size += len(block)
    if size == 0:
        raise ValueError(f'Empty file: {path.name}')
    return {'path': path.name, 'bytes': size, 'sha256': digest.hexdigest()}


def verify_manifest(data_dir: Path) -> dict:
    manifest = json.loads((data_dir / "manifest.json").read_text(encoding="utf-8"))
    records = manifest["files"]
    names = [record["path"] for record in records]
    if any(Path(name).name != name for name in names):
        raise ValueError("Source manifest paths must be plain filenames")
    if len(names) != len(set(names)):
        raise ValueError("Source manifest contains a duplicate filename")
    missing = {COUNTS, ANNOTATIONS} - set(names)
    if missing:
        raise ValueError("Source manifest is missing required inputs: " + ", ".join(sorted(missing)))
    for record in records:
        actual = fingerprint(data_dir / record["path"])
        if actual["bytes"] != record["bytes"] or actual["sha256"] != record["sha256"]:
            raise ValueError(f"Input fingerprint mismatch: {record['path']}")
    return manifest
