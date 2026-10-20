"""Snapshot and verify a complete, internally consistent pasilla result package.

This checks files and sample labels, not biological validity or numeric reproducibility.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

FILES = (
    'additive_design_matrix.csv', 'condition_read_type_counts.csv',
    'correlations_median_ratio.csv', 'correlations_total_count.csv', 'count_totals.png',
    'design_diagnostics.json', 'distances_median_ratio.csv', 'distances_total_count.csv',
    'leave_one_out_stability.csv', 'leave_one_out_stability.png',
    'normalization_comparison.png', 'normalization_distances.png', 'normalization_factors.csv',
    'pca_logcounts.csv', 'pca_logcounts.png', 'pca_logcpm.csv', 'pca_logcpm.png',
    'pca_median_ratio.csv', 'pca_total_count.csv', 'report.md', 'sample_correlations.csv',
    'sample_qc.csv', 'summary.json',
)
SAMPLE_TABLES = ('sample_qc.csv', 'pca_logcounts.csv', 'pca_logcpm.csv',
                 'pca_median_ratio.csv', 'pca_total_count.csv')


def fingerprint(path):
    path = Path(path)
    if not path.is_file() or path.is_symlink():
        raise ValueError(f'Expected a regular non-symlink file: {path.name}')
    digest = hashlib.sha256()
    size = 0
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block); size += len(block)
    if size == 0:
        raise ValueError(f'Empty file: {path.name}')
    return {'path': path.name, 'bytes': size, 'sha256': digest.hexdigest()}


def table(path):
    with Path(path).open(newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f'Invalid CSV header: {Path(path).name}')
        rows = list(reader)
    if not rows or any(None in row or None in row.values() for row in rows):
        raise ValueError(f'Empty or ragged CSV: {Path(path).name}')
    return rows


def snapshot(results, data):
    results, data = Path(results), Path(data)
    source = json.loads((data/'manifest.json').read_text())
    verified_inputs = {}
    for item in source['files']:
        if Path(item['path']).name != item['path']:
            raise ValueError('Source manifest paths must be plain filenames')
        actual = fingerprint(data/item['path'])
        if actual['bytes'] != item['bytes'] or actual['sha256'] != item['sha256']:
            raise ValueError(f"Input fingerprint mismatch: {item['path']}")
        verified_inputs[item['path']] = actual['sha256']
    summary = json.loads((results/'summary.json').read_text())
    if summary['input_sha256'] != verified_inputs or summary['upstream_commit'] != source['upstream_commit']:
        raise ValueError('Result summary does not match verified input provenance')
    with (data/'pasilla_gene_counts.tsv').open(newline='') as stream:
        header = next(csv.reader(stream, delimiter='\t'))
    samples = header[1:]
    if not samples or len(samples) != len(set(samples)) or summary['samples'] != len(samples):
        raise ValueError('Count header and result sample count disagree')
    annotations = table(data/'pasilla_sample_annotation.csv')
    annotation = {r['file'][:-2]: r for r in annotations if r['file'].endswith('fb')}
    if len(annotation) != len(annotations) or set(annotation) != set(samples):
        raise ValueError('Source annotation mapping is incomplete or duplicated')
    for filename in SAMPLE_TABLES:
        rows = table(results/filename)
        if [r.get('sample') for r in rows] != samples:
            raise ValueError(f'Sample order/membership differs in {filename}')
        for row in rows:
            original = annotation[row['sample']]
            if any(row.get(field) != original[field] for field in ('file', 'condition', 'type')):
                raise ValueError(f'Annotation differs in {filename}')
    return {
        'schema_version': '1.0', 'status': 'complete_package_with_matching_sample_labels',
        'scope': 'Byte integrity and sample-label consistency only; not biological validation or a rerun comparison',
        'annotation_status': 'Upstream package labels retained; article/package allocation discrepancy remains unresolved',
        'analysis_version': summary['analysis_version'], 'input_sha256': verified_inputs,
        'upstream_commit': source['upstream_commit'], 'samples_in_order': samples,
        'files': [fingerprint(results/name) for name in FILES],
    }


def create(results, data, manifest):
    record = snapshot(results, data)
    manifest = Path(manifest)
    protected = {Path(results).resolve()/name for name in FILES}
    protected.update(Path(data).resolve()/name for name in record['input_sha256'])
    if manifest.resolve() in protected:
        raise ValueError('Audit manifest cannot overwrite an input or result')
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open('x', encoding='utf-8') as stream:
        json.dump(record, stream, indent=2, sort_keys=True); stream.write('\n')
    return record


def verify(results, data, manifest):
    expected = json.loads(Path(manifest).read_text())
    actual = snapshot(results, data)
    if expected != actual:
        raise ValueError('Result package differs from the audit snapshot; do not mix run outputs')
    return actual


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['create', 'verify'])
    parser.add_argument('--results', type=Path, default=Path('results'))
    parser.add_argument('--data', type=Path, default=Path('data'))
    parser.add_argument('--manifest', type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        result = (create if args.mode == 'create' else verify)(args.results, args.data, args.manifest)
    except (ValueError, OSError, KeyError, TypeError) as error:
        print(f'Artifact audit failed: {error}', file=sys.stderr); return 2
    print(f"{args.mode}: {len(result['files'])} files and {len(result['samples_in_order'])} sample labels verified")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
