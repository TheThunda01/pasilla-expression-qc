"""Compare the standard numeric result tables while keeping identities exact.

This is a reproducibility check against a chosen reference, not scientific validation.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import sys

RTOL, ATOL = 1e-7, 1e-8
IDENTITY = {'sample', 'file', 'condition', 'type', 'omitted_sample', 'method'}
COUNTS = {'intercept', 'condition_treated', 'type_paired_end', 'single-read', 'paired-end',
          'remaining_samples', 'fixed_pca_genes', 'positive_reference_genes',
          'count_table_total', 'genes_with_positive_count'}
SCHEMAS = {
    'additive_design_matrix.csv': ['sample', 'intercept', 'condition_treated', 'type_paired_end'],
    'condition_read_type_counts.csv': ['condition', 'single-read', 'paired-end'],
    'leave_one_out_stability.csv': ['omitted_sample', 'method', 'remaining_samples', 'fixed_pca_genes', 'positive_reference_genes', 'pc1_variance_fraction', 'pc2_variance_fraction', 'procrustes_residual'],
    'normalization_factors.csv': ['sample', 'file', 'condition', 'type', 'total_count_factor', 'median_ratio_factor', 'median_to_total_factor_ratio', 'largest_gene_fraction'],
    'sample_qc.csv': ['sample', 'file', 'condition', 'type', 'count_table_total', 'zero_fraction_all_genes', 'genes_with_positive_count'],
}
for name in ('pca_logcounts', 'pca_logcpm', 'pca_median_ratio', 'pca_total_count'):
    SCHEMAS[name+'.csv'] = ['sample', 'file', 'condition', 'type', 'PC1', 'PC2']
MATRICES = ('correlations_median_ratio.csv', 'correlations_total_count.csv',
            'distances_median_ratio.csv', 'distances_total_count.csv', 'sample_correlations.csv')
JSON_FILES = ('summary.json', 'design_diagnostics.json')
JSON_KEYS = {
    'summary.json': {'analysis', 'analysis_version', 'condition_counts', 'design_diagnostics', 'filter', 'genes_after_filter', 'genes_before_filter', 'input_sha256', 'normalization', 'normalization_sensitivity', 'package_versions', 'pca', 'read_type_counts', 'samples', 'upstream_commit'},
    'design_diagnostics.json': {'all_design_cells_observed', 'all_intended_levels_observed', 'assessment', 'cell_counts', 'column_count', 'columns', 'equal_cell_counts', 'formula', 'full_column_rank', 'minimum_cell_size', 'missing_levels', 'model_fitted', 'perfect_confounding', 'rank', 'reference_levels', 'residual_dimension', 'samples', 'warnings'},
}
FILES = tuple(sorted((*SCHEMAS, *MATRICES, *JSON_FILES)))


def read_csv(path):
    with Path(path).open(newline='', encoding='utf-8') as stream:
        rows = list(csv.reader(stream, strict=True))
    if len(rows) < 2 or not rows[0] or len(set(rows[0])) != len(rows[0]):
        raise ValueError('Empty table or duplicate/empty header')
    if any(len(row) != len(rows[0]) for row in rows[1:]):
        raise ValueError('Ragged table')
    return rows[0], rows[1:]


def unique_pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError('Duplicate JSON key: '+key)
        result[key] = value
    return result


def read_json(path):
    def reject(value):
        raise ValueError('Nonfinite JSON number: '+value)
    result = json.loads(Path(path).read_text(), object_pairs_hook=unique_pairs, parse_constant=reject)
    if not isinstance(result, dict):
        raise ValueError('Expected a JSON object')
    def finite(value):
        if isinstance(value, float) and not math.isfinite(value):
            raise ValueError('Nonfinite JSON value')
        if isinstance(value, dict):
            for item in value.values(): finite(item)
        elif isinstance(value, list):
            for item in value: finite(item)
    finite(result)
    return result


def compare(reference, candidate, rtol=RTOL, atol=ATOL):
    if not math.isfinite(rtol) or not math.isfinite(atol) or rtol < 0 or atol < 0:
        raise ValueError('Tolerances must be finite and nonnegative')
    reference, candidate = Path(reference), Path(candidate)
    report = {'schema_version': '1.0', 'status': 'match', 'rtol': rtol, 'atol': atol,
              'criterion': 'abs(a-b) <= max(atol, rtol * max(abs(a), abs(b)))',
              'files_compared': 0, 'numeric_values_compared': 0, 'within_tolerance_nonidentical_values': 0, 'difference_count': 0,
              'differences': [], 'difference_details_limit': 100, 'environment_differences': [],
              'file_sha256': {'reference': {}, 'candidate': {}},
              'scope': '14 CSVs and 2 JSON files; identities/counts exact, finite floating values tolerant; no PNG comparison or biological validation'}

    def issue(location, reason, expected=None, actual=None):
        report['difference_count'] += 1
        if len(report['differences']) < report['difference_details_limit']:
            report['differences'].append({'location': location, 'reason': reason, 'expected': expected, 'actual': actual})

    def numeric(a, b, location):
        try:
            x, y = float(a), float(b)
        except (TypeError, ValueError, OverflowError):
            issue(location, 'Expected finite numeric value', a, b); return
        if not math.isfinite(x) or not math.isfinite(y):
            issue(location, 'Nonfinite numeric value', str(a), str(b)); return
        report['numeric_values_compared'] += 1
        if not math.isclose(x, y, rel_tol=rtol, abs_tol=atol):
            issue(location, 'Numeric difference exceeds tolerance', a, b)
        elif x != y:
            report['within_tolerance_nonidentical_values'] += 1

    def json_values(a, b, location):
        if type(a) is not type(b):
            issue(location, 'JSON type changed', type(a).__name__, type(b).__name__); return
        if isinstance(a, dict):
            if set(a) != set(b):
                issue(location, 'JSON keys changed', sorted(a), sorted(b))
            for key in sorted(set(a) & set(b)):
                if location == 'summary.json' and key == 'package_versions':
                    if (not isinstance(a[key], dict) or not isinstance(b[key], dict)
                            or set(a[key]) != set(b[key])
                            or any(not isinstance(v, str) for obj in (a[key], b[key]) for v in obj.values())):
                        issue(location+'.'+key, 'Environment schema changed', a[key], b[key])
                    elif a[key] != b[key]:
                        report['environment_differences'].append({'reference': a[key], 'candidate': b[key]})
                else:
                    json_values(a[key], b[key], location+'.'+key)
        elif isinstance(a, list):
            if len(a) != len(b):
                issue(location, 'JSON list length changed', len(a), len(b))
            for i, (left, right) in enumerate(zip(a, b)):
                json_values(left, right, f'{location}[{i}]')
        elif isinstance(a, float):
            numeric(a, b, location)
        elif a != b:
            issue(location, 'Exact JSON value changed', a, b)

    try:
        header, rows = read_csv(reference/'sample_qc.csv')
        if header != SCHEMAS['sample_qc.csv']:
            raise ValueError('Reference sample QC schema changed')
        samples = [row[0] for row in rows]
        if any(not s for s in samples) or len(set(samples)) != len(samples):
            raise ValueError('Reference sample identities must be nonempty and unique')
    except (ValueError, OSError, csv.Error) as error:
        issue('sample_qc.csv', 'Cannot establish reference identities: '+str(error))
        report['status'] = 'different'
        return report

    for filename in FILES:
        parsed = []
        for side, directory in [('reference', reference), ('candidate', candidate)]:
            try:
                path = directory/filename
                report['file_sha256'][side][filename] = hashlib.sha256(path.read_bytes()).hexdigest()
                parsed.append(read_json(path) if filename in JSON_FILES else read_csv(path))
            except (ValueError, OSError, csv.Error) as error:
                issue(filename, side+' parse/read failure: '+str(error)); parsed.append(None)
        if any(value is None for value in parsed):
            continue
        report['files_compared'] += 1
        left, right = parsed
        if filename in JSON_FILES:
            if set(left) != JSON_KEYS[filename] or set(right) != JSON_KEYS[filename]:
                issue(filename, 'JSON root schema changed', sorted(JSON_KEYS[filename]), {'reference': sorted(left), 'candidate': sorted(right)})
            json_values(left, right, filename)
            continue
        expected_header = ['sample', *samples] if filename in MATRICES else SCHEMAS[filename]
        a_header, a_rows = left; b_header, b_rows = right
        if a_header != expected_header or b_header != expected_header:
            issue(filename, 'CSV schema or column order changed', expected_header, {'reference': a_header, 'candidate': b_header}); continue
        if len(a_rows) != len(b_rows):
            issue(filename, 'CSV row count changed', len(a_rows), len(b_rows))
        if filename in MATRICES:
            for side, table in [('reference', a_rows), ('candidate', b_rows)]:
                if [r[0] for r in table] != samples:
                    issue(filename, side+' matrix row identities/order changed', samples, [r[0] for r in table])
        for row_number, (a, b) in enumerate(zip(a_rows, b_rows), 2):
            for field, x, y in zip(expected_header, a, b):
                location = f'{filename}:{row_number}:{field}'
                if field in IDENTITY:
                    if not x or not y or x != y:
                        issue(location, 'Exact identity/category changed', x, y)
                elif field in COUNTS:
                    if not re.fullmatch(r'\d+', x) or not re.fullmatch(r'\d+', y) or int(x) != int(y):
                        issue(location, 'Nonnegative integer count changed or malformed', x, y)
                else:
                    numeric(x, y, location)
    report['status'] = 'different' if report['difference_count'] else 'match'
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', type=Path)
    parser.add_argument('candidate', type=Path)
    parser.add_argument('--rtol', type=float, default=RTOL)
    parser.add_argument('--atol', type=float, default=ATOL)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args(argv)
    try:
        protected = {(directory/name).resolve() for directory in (args.reference, args.candidate) for name in FILES}
        if args.report and args.report.resolve() in protected:
            raise ValueError('Comparison report cannot overwrite an input result')
        result = compare(args.reference, args.candidate, args.rtol, args.atol)
        text = json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n'
        if args.report:
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(text)
        print(text, end='')
    except (ValueError, OSError) as error:
        print(f'Comparison failed: {error}', file=sys.stderr); return 2
    return 0 if result['status'] == 'match' else 1


if __name__ == '__main__':
    raise SystemExit(main())
