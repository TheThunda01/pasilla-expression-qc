"""Link saved sample coordinates, scaling factors and omission sensitivity."""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import itertools
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METHODS = {'total_count': 'Total-count', 'median_ratio': 'Median-ratio'}


def load_results(directory: Path) -> tuple[dict, dict]:
    data = {}
    fingerprints = {}
    for name in ('pca_total_count.csv', 'pca_median_ratio.csv', 'normalization_factors.csv', 'leave_one_out_stability.csv'):
        path = directory / name
        fingerprints[name] = hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
        with path.open(encoding='utf-8-sig', newline='') as handle:
            data[name] = list(csv.DictReader(handle))
    samples = {}
    for raw in data['normalization_factors.csv']:
        sample = raw['sample']
        if sample in samples or raw['condition'] not in {'treated', 'untreated'} or raw['type'] not in {'single-read', 'paired-end'} or raw['file'] != sample + 'fb':
            raise ValueError('Invalid or duplicate sample metadata')
        row = dict(raw)
        for field in ('total_count_factor', 'median_ratio_factor', 'median_to_total_factor_ratio'):
            row[field] = float(raw[field])
            if not math.isfinite(row[field]) or row[field] <= 0:
                raise ValueError('Normalization factors must be finite and positive')
        if not math.isclose(row['median_to_total_factor_ratio'], row['median_ratio_factor'] / row['total_count_factor'], rel_tol=1e-8, abs_tol=1e-8):
            raise ValueError('Recorded ratio disagrees with its factors')
        samples[sample] = row
    if len(samples) != 7:
        raise ValueError('The bundled pasilla example requires seven samples')
    coordinates = {}
    for method in METHODS:
        coordinates[method] = {}
        for raw in data[f'pca_{method}.csv']:
            sample = raw['sample']
            if sample in coordinates[method] or sample not in samples:
                raise ValueError('Duplicate or unsupported coordinate sample')
            if any(raw[k] != samples[sample][k] for k in ('file', 'condition', 'type')):
                raise ValueError('Sample metadata disagree across tables')
            row = dict(raw, PC1=float(raw['PC1']), PC2=float(raw['PC2']))
            if not all(math.isfinite(row[k]) for k in ('PC1', 'PC2')):
                raise ValueError('PCA coordinates must be finite')
            coordinates[method][sample] = row
        if set(coordinates[method]) != set(samples):
            raise ValueError('Both coordinate tables must cover every sample')
    omissions = {}
    expected = set(itertools.product(samples, METHODS))
    for raw in data['leave_one_out_stability.csv']:
        key = (raw['omitted_sample'], raw['method'])
        if key in omissions or key not in expected or int(raw['remaining_samples']) != 6 or int(raw['fixed_pca_genes']) != 8148:
            raise ValueError('Invalid omission identity or analysis dimensions')
        residual = float(raw['procrustes_residual'])
        if not math.isfinite(residual) or not 0 <= residual <= 2:
            raise ValueError('Invalid normalized Procrustes residual')
        omissions[key] = dict(raw, procrustes_residual=residual)
    if set(omissions) != expected:
        raise ValueError('Every sample needs both omission methods')
    return dict(samples=samples, coordinates=coordinates, omissions=omissions), fingerprints


def render(data: dict, fingerprints: dict, output: Path) -> None:
    import plotly
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    samples, coordinates, omissions = data['samples'], data['coordinates'], data['omissions']
    figure = make_subplots(rows=2, cols=2, subplot_titles=(
        'Total-count normalization', 'Median-ratio normalization',
        'Sensitivity when each sample is omitted', 'Median-ratio / total-count factor'),
        vertical_spacing=0.22, horizontal_spacing=0.15)
    colors = {'untreated': '#245b82', 'treated': '#b24a2d'}
    for col, method in enumerate(METHODS, 1):
        for condition, read_type in itertools.product(colors, ('single-read', 'paired-end')):
            selected = [r for r in coordinates[method].values() if r['condition'] == condition and r['type'] == read_type]
            custom = [[r['sample'], r['condition'], r['type'], samples[r['sample']]['median_to_total_factor_ratio'], omissions[r['sample'], method]['procrustes_residual']] for r in selected]
            figure.add_trace(go.Scatter(x=[r['PC1'] for r in selected], y=[r['PC2'] for r in selected], mode='markers+text',
                text=[r['sample'] for r in selected], textposition=['bottom center' if r['sample'] in {'untreated1', 'untreated4', 'treated1', 'treated2'} else 'top center' for r in selected], name=f'{condition} · {read_type}',
                legendgroup=f'{condition}-{read_type}', showlegend=col == 1, customdata=custom,
                marker=dict(color=colors[condition], symbol='circle' if read_type == 'single-read' else 'diamond', size=13),
                hovertemplate='%{customdata[0]}<br>%{customdata[1]} · %{customdata[2]}<br>PC1 %{x:.4f}; PC2 %{y:.4f}<br>Factor ratio %{customdata[3]:.5f}<br>Residual when omitted %{customdata[4]:.5f}<extra>' + METHODS[method] + '</extra>'), row=1, col=col)
        figure.update_xaxes(title_text='PC1 score', row=1, col=col)
        figure.update_yaxes(title_text='PC2 score', row=1, col=col)
    xs = [r['PC1'] for table in coordinates.values() for r in table.values()]
    ys = [r['PC2'] for table in coordinates.values() for r in table.values()]
    for col in (1, 2):
        figure.update_xaxes(range=[min(xs) - 7, max(xs) + 7], row=1, col=col)
        figure.update_yaxes(range=[min(ys) - 8, max(ys) + 8], scaleanchor='x' if col == 1 else 'x2', scaleratio=1, constrain='domain', row=1, col=col)
    for method, color in (('total_count', '#5b798a'), ('median_ratio', '#76558b')):
        figure.add_trace(go.Bar(x=list(samples), y=[omissions[sample, method]['procrustes_residual'] for sample in samples],
            name=METHODS[method], marker_color=color,
            customdata=[[samples[sample]['condition'], samples[sample]['type']] for sample in samples],
            hovertemplate='Omit %{x}<br>%{customdata[0]} · %{customdata[1]}<br>Residual %{y:.5f}<extra>' + METHODS[method] + '</extra>'), row=2, col=1)
    figure.add_trace(go.Scatter(x=list(samples), y=[r['median_to_total_factor_ratio'] for r in samples.values()], mode='markers', showlegend=False,
        marker=dict(size=12, color=[colors[r['condition']] for r in samples.values()], symbol=['circle' if r['type'] == 'single-read' else 'diamond' for r in samples.values()]),
        customdata=[[r['condition'], r['type'], r['total_count_factor'], r['median_ratio_factor']] for r in samples.values()],
        hovertemplate='%{x}<br>%{customdata[0]} · %{customdata[1]}<br>Ratio %{y:.5f}<br>Total-count factor %{customdata[2]:.5f}<br>Median-ratio factor %{customdata[3]:.5f}<extra></extra>'), row=2, col=2)
    figure.add_hline(y=1, line_dash='dot', line_color='#64748b', row=2, col=2)
    figure.update_yaxes(title_text='Normalized Procrustes residual', rangemode='tozero', row=2, col=1)
    figure.update_yaxes(title_text='Factor ratio (1 = equal)', row=2, col=2)
    figure.update_xaxes(tickangle=-40, row=2, col=1)
    figure.update_xaxes(tickangle=-40, row=2, col=2)
    figure.update_layout(template='plotly_white', width=1250, height=1115, barmode='group',
        title='Seven samples: scaling, geometry and influence', font=dict(family='Arial', size=13),
        margin=dict(l=80, r=50, t=100, b=260), legend=dict(orientation='h', y=-0.15, yanchor='top'))
    figure.add_annotation(x=0, y=-0.23, xref='paper', yref='paper', xanchor='left', yanchor='top', showarrow=False, align='left', font=dict(size=12),
        text='Labels follow the package and remain unresolved against the source history; all seven samples are retained.<br>Geometry sensitivity is not a sample-rejection rule or a treatment-effect test.')
    caption = ('These views read the saved numerical tables and retain all seven samples. Package condition/read-type labels remain unresolved against the original source history. '
        'Each PCA panel uses its saved coordinate system; equal displayed ranges do not make axes biological measurements or independently align the two fits. '
        'Omission residuals measure geometry changes after alignment, with 8,148 genes and the reference depth held fixed. '
        'treated1 is the only treated single-read sample; influence is not a rejection criterion or treatment-effect test. '
        'Hover on samples to connect their labels, scaling factor ratio and omission residual. No model was refit.')
    output.mkdir(parents=True, exist_ok=True)
    fragment = figure.to_html(full_html=False, include_plotlyjs=True, div_id='pasilla-sample-explorer', config={'displaylogo': False, 'responsive': True})
    content = '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Pasilla sample explorer</title><style>body{font:16px/1.55 Arial,sans-serif;margin:24px auto;max-width:1300px;padding:0 20px;color:#20303c}pre{white-space:pre-wrap;font-size:12px}</style><h1>Inspect sample sensitivity</h1><p>' + caption + '</p>' + fragment + '<h2>Source fingerprints (LF-normalized SHA-256)</h2><pre>' + html.escape(json.dumps(fingerprints, indent=2)) + '</pre></html>'
    (output / 'sample-explorer.html').write_text(content, encoding='utf-8')
    figure.write_image(str(output / 'sample-explorer.png'), scale=1.5)
    (output / 'source-manifest.json').write_text(json.dumps({'schema_version': 1, 'plotly_version': plotly.__version__, 'source_sha256_lf_normalized': fingerprints, 'samples': len(samples), 'omission_rows': len(omissions), 'scope': caption}, indent=2) + '\n', encoding='utf-8')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-dir', type=Path, default=ROOT / 'results')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'docs/figures/sample-explorer')
    parser.add_argument('--check-only', action='store_true')
    args = parser.parse_args()
    data, fingerprints = load_results(args.results_dir)
    if not args.check_only:
        render(data, fingerprints, args.output_dir)
    print(json.dumps({'samples': len(data['samples']), 'omission_rows': len(data['omissions']), 'validated': True}))


if __name__ == '__main__':
    main()
