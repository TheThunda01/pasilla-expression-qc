"""Descriptive QC for a pinned public Drosophila count matrix.

This module deliberately performs no differential-expression inference.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from sensitivity import principal_components, run_sensitivity
from design import inspect_design

ROOT = Path(__file__).resolve().parent
MIN_COUNT = 10
MIN_SAMPLES = 3


def verify_manifest(data_dir: Path) -> dict:
    """Fail before analysis if any bundled upstream file has changed."""
    manifest = json.loads((data_dir / "manifest.json").read_text())
    for record in manifest["files"]:
        path = data_dir / record["path"]
        raw = path.read_bytes()
        if len(raw) != record["bytes"] or hashlib.sha256(raw).hexdigest() != record["sha256"]:
            raise ValueError(f"Provenance check failed for {record['path']}")
    return manifest


def load_counts(path: Path) -> pd.DataFrame:
    """Read without allowing pandas to silently rename duplicated headers."""
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if len(header) < 3 or header[0] != "gene_id":
            raise ValueError("Counts need gene_id and at least two sample columns")
        samples = header[1:]
        if len(samples) != len(set(samples)) or any(not s.strip() for s in samples):
            raise ValueError("Sample names must be nonempty and unique")
        genes, values = [], []
        seen = set()
        for line, row in enumerate(reader, 2):
            if len(row) != len(header) or not row[0].strip() or row[0] in seen:
                raise ValueError(f"Invalid or duplicate gene row at line {line}")
            try:
                # Decimal strings, NaN and scientific notation are not integer counts.
                if any(not x.isascii() or not x.isdecimal() for x in row[1:]):
                    raise ValueError
                numbers = [int(x) for x in row[1:]]
                if any(x > np.iinfo(np.int64).max for x in numbers):
                    raise ValueError
            except ValueError as exc:
                raise ValueError(f"Nonnegative integer counts required at line {line}") from exc
            seen.add(row[0])
            genes.append(row[0])
            values.append(numbers)
    if not genes:
        raise ValueError("Count matrix has no genes")
    frame = pd.DataFrame(values, index=genes, columns=samples, dtype=np.int64)
    # Use float64 sum to avoid integer overflow for user-provided matrices.
    if np.any(frame.to_numpy(dtype=float).sum(axis=0) <= 0):
        raise ValueError("Each sample must have a positive library count total")
    return frame


def load_metadata(path: Path, samples: list[str]) -> pd.DataFrame:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if (not {"file", "condition", "type"}.issubset(reader.fieldnames or [])
                or len(reader.fieldnames) != len(set(reader.fieldnames))):
            raise ValueError("Metadata needs file, condition and type")
        records = list(reader)
    if not records:
        raise ValueError("Metadata is empty")
    if any(None in row or any(value is None for value in row.values()) for row in records):
        raise ValueError("Metadata rows must match the header")
    names = [r["file"] for r in records]
    if any(not n or not n.endswith("fb") for n in names):
        raise ValueError("Expected the documented pasilla file suffix 'fb'")
    mapped = [name[:-2] for name in names]
    if len(set(mapped)) != len(mapped):
        raise ValueError("Metadata sample identifiers must be unique")
    if set(mapped) != set(samples):
        raise ValueError("Metadata and counts must contain exactly the same samples")
    for record in records:
        if record["condition"] not in {"treated", "untreated"}:
            raise ValueError("Unknown condition in pasilla metadata")
        if record["type"] not in {"single-read", "paired-end"}:
            raise ValueError("Unknown read type in pasilla metadata")
    frame = pd.DataFrame(records, index=mapped)
    return frame.loc[samples, ["file", "condition", "type"]]


def make_pca_plot(scores: np.ndarray, ratios: np.ndarray, metadata: pd.DataFrame,
                  title: str, path: Path) -> None:
    colors = {"treated": "#b34d20", "untreated": "#1b6a8a"}
    shapes = {"single-read": "o", "paired-end": "s"}
    label_offsets = {"untreated1": (-8, -18), "untreated2": (8, 8),
                     "untreated3": (8, 8), "untreated4": (8, -18),
                     "treated1": (8, 8), "treated2": (8, -16), "treated3": (8, 8)}
    fig, ax = plt.subplots(figsize=(8, 5.5), layout="constrained")
    for i, (sample, row) in enumerate(metadata.iterrows()):
        ax.scatter(*scores[i], c=colors[row["condition"]], marker=shapes[row["type"]], s=70)
        offset = label_offsets.get(sample, (8, 8))
        ax.annotate(sample, scores[i], xytext=offset, textcoords="offset points", fontsize=9,
                    ha="right" if offset[0] < 0 else "left")
    ax.margins(.22)
    ax.set(xlabel=f"PC1 ({ratios[0]:.1%} of variance)",
           ylabel=f"PC2 ({ratios[1]:.1%} of variance)", title=title)
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=c, label=l)
               for l, c in colors.items()]
    handles += [Line2D([0], [0], marker=m, color="#666666", linestyle="none", label=l)
                for l, m in shapes.items()]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1))
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(path, dpi=160, metadata={"Software": "pasilla-expression-qc"})
    plt.close(fig)


def run(data_dir: Path, out: Path) -> dict:
    manifest = verify_manifest(data_dir)
    counts = load_counts(data_dir / "pasilla_gene_counts.tsv")
    metadata = load_metadata(data_dir / "pasilla_sample_annotation.csv", list(counts.columns))
    totals = counts.to_numpy(dtype=float).sum(axis=0)
    keep = (counts >= MIN_COUNT).sum(axis=1) >= MIN_SAMPLES
    retained = counts.loc[keep]
    if len(retained) < 2:
        raise ValueError("Fewer than two genes pass the exploratory filter")
    # Full count-table totals are the denominators, before filtering.
    log_cpm = np.log2(retained.to_numpy(dtype=float) / totals[None, :] * 1e6 + 1)
    log_raw = np.log2(retained.to_numpy(dtype=float) + 1)
    scores, ratios = principal_components(log_cpm.T)
    raw_scores, raw_ratios = principal_components(log_raw.T)
    out.mkdir(parents=True, exist_ok=True)
    cross, design_matrix, design = inspect_design(metadata)
    cross.to_csv(out / "condition_read_type_counts.csv")
    design_matrix.to_csv(out / "additive_design_matrix.csv")
    (out / "design_diagnostics.json").write_text(json.dumps(design, indent=2, sort_keys=True) + "\n")
    sensitivity = run_sensitivity(counts, keep.to_numpy(), metadata, out)
    qc = metadata.copy()
    qc.index.name = "sample"
    qc["count_table_total"] = totals.astype(np.int64)
    qc["zero_fraction_all_genes"] = (counts == 0).mean(axis=0)
    qc["genes_with_positive_count"] = (counts > 0).sum(axis=0)
    qc.to_csv(out / "sample_qc.csv", float_format="%.9f")
    for name, values in [("pca_logcpm", scores), ("pca_logcounts", raw_scores)]:
        result = metadata.copy()
        result.index.name = "sample"
        result["PC1"], result["PC2"] = values[:, 0], values[:, 1]
        result.to_csv(out / f"{name}.csv", float_format="%.9f")
    correlation = pd.DataFrame(np.corrcoef(log_cpm.T), index=counts.columns, columns=counts.columns)
    correlation.index.name = "sample"
    correlation.to_csv(out / "sample_correlations.csv", float_format="%.9f")
    make_pca_plot(scores, ratios, metadata, "Pasilla: PCA of log2(CPM + 1)", out / "pca_logcpm.png")
    make_pca_plot(raw_scores, raw_ratios, metadata, "Sensitivity check: log2(count + 1)", out / "pca_logcounts.png")
    fig, ax = plt.subplots(figsize=(8, 4.5), layout="constrained")
    ax.bar(counts.columns, totals / 1e6, color=["#b34d20" if c == "treated" else "#1b6a8a"
                                             for c in metadata.condition])
    ax.set(ylabel="Summed gene counts (millions)", title="Count-table totals, before filtering")
    ax.tick_params(axis="x", rotation=25)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(out / "count_totals.png", dpi=160, metadata={"Software": "pasilla-expression-qc"})
    plt.close(fig)
    summary = {
        "analysis": "Descriptive QC only; no differential expression or causal inference",
        "analysis_version": "0.3.0",
        "design_diagnostics": design,
        "normalization_sensitivity": sensitivity,
        "genes_before_filter": int(len(counts)), "genes_after_filter": int(keep.sum()),
        "samples": int(counts.shape[1]),
        "filter": {"minimum_count": MIN_COUNT, "minimum_samples": MIN_SAMPLES},
        "normalization": "CPM using unfiltered count-table totals, then log2(CPM + 1)",
        "pca": {"center_genes": True, "scale_genes": False,
                "logcpm_explained_variance_ratio": ratios.tolist(),
                "logcounts_explained_variance_ratio": raw_ratios.tolist()},
        "condition_counts": metadata.condition.value_counts().sort_index().to_dict(),
        "read_type_counts": metadata.type.value_counts().sort_index().to_dict(),
        "upstream_commit": manifest["upstream_commit"],
        "input_sha256": {r["path"]: r["sha256"] for r in manifest["files"]},
        "package_versions": {p: importlib.metadata.version(p) for p in ["numpy", "pandas", "matplotlib"]},
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    report = f"""# Pasilla expression QC: example run

This run uses {len(counts):,} genes and {counts.shape[1]} samples from the pinned public count table. {keep.sum():,} genes pass the exploratory filter (at least {MIN_COUNT} counts in at least {MIN_SAMPLES} samples).

## Sample design before interpreting the plots

| Condition | Single-read | Paired-end |
| --- | ---: | ---: |
| Untreated | {int(cross.loc['untreated', 'single-read'])} | {int(cross.loc['untreated', 'paired-end'])} |
| Treated | {int(cross.loc['treated', 'single-read'])} | {int(cross.loc['treated', 'paired-end'])} |

The additive design uses an intercept, a treated indicator, and a paired-end indicator. Its rank is **{design['rank']} of {design['column_count']} columns**, with residual dimension **{design['residual_dimension']}**. Assessment: **{design['assessment']}**.

The smallest intended condition/read-type cell contains {design['minimum_cell_size']} sample(s). Full rank makes the two coded contrasts algebraically distinguishable; it does not establish adequate power, eliminate other confounding, or validate an additive model. No expression model, treatment effect, or p-value is fitted here. See `additive_design_matrix.csv` and `design_diagnostics.json` for the exact coding and checks.

## Descriptive PCA views

The $`\\log_2(\\mathrm{{CPM}}+1)`$ PCA assigns {ratios[0]:.1%} of variance to PC1 and {ratios[1]:.1%} to PC2. The $`\\log_2(\\text{{count}}+1)`$ sensitivity analysis assigns {raw_ratios[0]:.1%} and {raw_ratios[1]:.1%}, respectively. These percentages describe two different transformed data matrices; they are not an accuracy comparison.

![PCA after total-count normalization](pca_logcpm.png)

![PCA without total-count normalization](pca_logcounts.png)

The plots display both condition and read type so potential technical structure is visible. They do not isolate the cause of sample separation or establish differential expression. No sample is dropped based on these plots.

![Count-table totals](count_totals.png)

Count-table totals sum the supplied gene counts. They are not the original total sequenced-read counts. The full numeric results are in `sample_qc.csv`, `sample_correlations.csv`, and `summary.json`.
"""
    report += f"""

## Normalization sensitivity with a shared count scale

The new comparison holds the retained genes, reference count depth ({sensitivity['reference_depth']:,.1f}), and +1 offset fixed. It contrasts total-count scaling with an independently implemented median-ratio calculation using {sensitivity['positive_reference_genes']:,} genes positive in all seven samples. The two factor vectors have geometric mean one. Their median-ratio/total-count factor ratios range from {sensitivity['factor_ratio_range'][0]:.3f} to {sensitivity['factor_ratio_range'][1]:.3f}.

![Comparable normalized PCA views](normalization_comparison.png)

The normalized two-dimensional Procrustes residual between the views is {sensitivity['between_method_procrustes_residual']:.4f}. Zero would mean identical plotted geometry after translation, overall scaling, rotation and reflection. There is no universal acceptance threshold, and neither close agreement nor disagreement establishes which normalization is biologically correct.

![Comparable sample distances](normalization_distances.png)

These distances summarize root-mean-square differences over the retained log2 gene profiles. They use all retained genes, rather than only the first two PCs. All samples remain included. The numeric distance and correlation matrices are provided separately for each method.

## How much does one sample influence the PCA picture?

The leave-one-out study refits normalization and PCA after each omission, while holding the original gene set and reference depth fixed. It compares the remaining samples with those same samples in the full-data two-dimensional picture.

![Leave-one-sample-out sensitivity](leave_one_out_stability.png)

For total-count scaling, residuals range from {sensitivity['leave_one_out_residuals']['total_count']['minimum']:.4f} to {sensitivity['leave_one_out_residuals']['total_count']['maximum']:.4f}; the largest change follows omission of `{sensitivity['leave_one_out_residuals']['total_count']['largest_change_omitted_sample']}`. For median-ratio scaling, they range from {sensitivity['leave_one_out_residuals']['median_ratio']['minimum']:.4f} to {sensitivity['leave_one_out_residuals']['median_ratio']['maximum']:.4f}; the largest change follows omission of `{sensitivity['leave_one_out_residuals']['median_ratio']['largest_change_omitted_sample']}`. An influential sample is not thereby an erroneous sample. These are descriptive refits, not cross-validation of a predictive model or evidence of causation.

See `normalization_factors.csv`, `leave_one_out_stability.csv`, and the [methods](../docs/methods.md) for exact definitions and limitations. The older CPM and raw-log views above use different overall count scales; the shared-scale comparison isolates normalization choice more carefully.
"""
    (out / "report.md").write_text(report)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--out", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    try:
        summary = run(args.data_dir, args.out)
    except (ValueError, OSError, KeyError, json.JSONDecodeError) as exc:
        print(f"Analysis stopped: {exc}", file=sys.stderr)
        return 2
    print(f"Validated {summary['genes_before_filter']:,} genes and {summary['samples']} samples; "
          f"retained {summary['genes_after_filter']:,} genes. Results: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
