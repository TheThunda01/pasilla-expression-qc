# Inspect sample sensitivity together

The [interactive sample explorer](figures/sample-explorer/sample-explorer.html) puts the saved PCA coordinates, normalization factor ratios, and omission residuals in one view. Download the HTML using GitHub's **Download raw file** control and open it locally. Its plotting library is included, so it needs no server or network connection.

![Sample geometry, factor ratios and omission sensitivity](figures/sample-explorer/sample-explorer.png)

Hover on either PCA panel to see a sample's condition/read-type label, factor ratio and residual when that sample is omitted. Color identifies the package condition label; circle/diamond identifies single-read/paired-end. All seven samples remain visible. The lower panels expose the scale-factor ratio and both omission methods for each sample.

The PCA panels use equal coordinate units and shared displayed ranges. They retain their separately fitted coordinate systems; the builder does not independently align the two PCA fits or reinterpret an axis as a biological measurement. The saved omission residuals already use the alignment described in [Methods](methods.md).

`treated1` is the only treated single-read sample, so its omission changes sample composition as well as geometry. A large residual is not a sample-rejection criterion or a treatment-effect test. The original package labels are preserved, including their [unresolved source provenance](review-2026-10-21.md#source-reconciliation). This view does not resolve that discrepancy.

## Input and output contract

The [builder](../scripts/build_sample_explorer.py) reads four saved tables:

- `results/pca_total_count.csv` and `results/pca_median_ratio.csv`: both complete seven-sample coordinate sets.
- `results/normalization_factors.csv`: metadata, positive factors and their ratio.
- `results/leave_one_out_stability.csv`: both methods for all seven omissions, with six remaining samples and the fixed 8,148-gene PCA set.

Sample identities, file names, conditions and read types must agree across tables. Coordinates and omission residuals must be finite; factors must be positive and reproduce their recorded ratio within the tolerance needed for the nine-decimal CSV values. Missing/duplicate samples, incomplete omissions, invalid counts or contradictory metadata stop generation.

Outputs are `sample-explorer.html`, `sample-explorer.png`, and [source-manifest.json](figures/sample-explorer/source-manifest.json) under `docs/figures/sample-explorer/`. Source SHA-256 fingerprints normalize CRLF to LF only, making Git text checkout conventions explicit. Nothing under the audited `results/` package is changed; no normalization or PCA is rerun by this builder.

## Regenerate and verify

Table validation uses the Python standard library:

```sh
python scripts/build_sample_explorer.py --check-only
python -m unittest discover -s tests -v
python audit_results.py verify --manifest audits/example-artifacts.json
```

Use a separate environment for optional rendering:

```sh
python -m venv .venv-visualization
# macOS/Linux: source .venv-visualization/bin/activate
# Windows PowerShell: .venv-visualization\Scripts\Activate.ps1
python -m pip install -r requirements-visualization.txt
python scripts/build_sample_explorer.py
```

The optional renderer pins Plotly 7.1.0 and Kaleido 1.5.0; the generated HTML retains the embedded Plotly.js MIT license header. PNG export also needs an existing Chrome installation. See [Plotly's official export guide](https://plotly.com/python/static-image-export/). The builder does not download a browser. Use `--results-dir` for compatible rerun tables and `--output-dir` to keep a separate visualization; its three named output files are replaced on regeneration. Check the command's exit status before treating a bundle as complete.

Tests check missing coordinate rows, inconsistent read-type labels, missing omissions, duplicate samples, nonfinite coordinates and contradictory factor ratios. These validate the visualization's input contract. They do not independently establish the source labels or biological interpretation.
