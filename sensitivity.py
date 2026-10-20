"""Independent mathematical implementations for descriptive QC sensitivity checks.

No DESeq2 code is included or called. See docs/source-notes.md for method sources
and the distinction between this calculation and a fitted count-data model.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

METHODS = ("total_count", "median_ratio")
METHOD_LABELS = {"total_count": "Total-count scaling", "median_ratio": "Median-ratio scaling"}


def validate_matrix(counts: np.ndarray) -> np.ndarray:
    matrix = np.asarray(counts, dtype=float)
    if matrix.ndim != 2 or min(matrix.shape) < 2:
        raise ValueError("Need at least two genes and two samples")
    if not np.isfinite(matrix).all() or (matrix < 0).any():
        raise ValueError("Normalization needs finite nonnegative counts")
    if (matrix.sum(axis=0) <= 0).any():
        raise ValueError("Normalization needs positive library totals")
    return matrix


def geometric_mean(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if values.size == 0 or not np.isfinite(values).all() or (values <= 0).any():
        raise ValueError("Geometric means need finite positive values")
    return float(np.exp(np.mean(np.log(values))))


def size_factors(counts: np.ndarray) -> tuple[dict[str, np.ndarray], int]:
    """Return both factor vectors with geometric mean one.

    The median-ratio estimate uses only genes positive in every sample and is
    exp(median(log count - mean log count)). No pseudocount is used here.
    """
    matrix = validate_matrix(counts)
    totals = matrix.sum(axis=0)
    positive = (matrix > 0).all(axis=1)
    if not positive.any():
        raise ValueError("No gene is positive in every sample; median-ratio estimation is undefined")
    logs = np.log(matrix[positive])
    log_ratios = logs - logs.mean(axis=1, keepdims=True)
    median_logs = np.median(log_ratios, axis=0)
    median = np.exp(median_logs - median_logs.mean())
    return {"total_count": totals / geometric_mean(totals), "median_ratio": median}, int(positive.sum())


def transformed_counts(counts: np.ndarray, factors: np.ndarray,
                       reference_depth: float) -> np.ndarray:
    """Anchor either factor vector to a declared common count scale.

    Factors have geometric mean one in the supplied sample subset. Multiplying
    by subset geometric-mean library depth / fixed reference depth avoids an
    arbitrary global rescaling when a sample is omitted.
    """
    matrix = validate_matrix(counts)
    factors = np.asarray(factors, dtype=float)
    if factors.shape != (matrix.shape[1],) or not np.isfinite(factors).all() or (factors <= 0).any():
        raise ValueError("One finite positive size factor is required per sample")
    if not np.isfinite(reference_depth) or reference_depth <= 0:
        raise ValueError("Reference depth must be finite and positive")
    scale = geometric_mean(matrix.sum(axis=0)) / reference_depth
    return np.log2(matrix / (factors[None, :] * scale) + 1)


def principal_components(sample_by_gene: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Center genes without scaling their variance; retain two sample PCs."""
    matrix = np.asarray(sample_by_gene, dtype=float)
    if matrix.ndim != 2 or min(matrix.shape) < 2 or not np.isfinite(matrix).all():
        raise ValueError("PCA needs a finite matrix with at least two samples and genes")
    centered = matrix - matrix.mean(axis=0)
    u, singular, vt = np.linalg.svd(centered, full_matrices=False)
    total = np.square(singular).sum()
    if total <= 0:
        raise ValueError("PCA needs nonzero variation")
    for i in range(2):
        if vt[i, np.argmax(np.abs(vt[i]))] < 0:
            u[:, i] *= -1
    return u[:, :2] * singular[:2], np.square(singular[:2]) / total


def pairwise_summaries(sample_by_gene: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    matrix = np.asarray(sample_by_gene, dtype=float)
    if matrix.ndim != 2 or min(matrix.shape) < 2 or not np.isfinite(matrix).all():
        raise ValueError("Distances need finite sample-by-gene values")
    if np.any(np.std(matrix, axis=1) == 0):
        raise ValueError("Pearson correlation is undefined for a constant sample profile")
    distances = np.sqrt(np.mean(np.square(matrix[:, None, :] - matrix[None, :, :]), axis=2))
    return distances, np.corrcoef(matrix)


def procrustes_residual(reference: np.ndarray, candidate: np.ndarray) -> float:
    """Compare 2-D shape after centering, unit Frobenius norm and rotation/reflection."""
    reference, candidate = np.asarray(reference, dtype=float), np.asarray(candidate, dtype=float)
    if (reference.shape != candidate.shape or reference.ndim != 2
            or reference.shape[1] != 2 or reference.shape[0] < 3
            or not np.isfinite(reference).all() or not np.isfinite(candidate).all()):
        raise ValueError("Procrustes comparison needs matching finite n-by-2 matrices, n >= 3")
    a = reference - reference.mean(axis=0)
    b = candidate - candidate.mean(axis=0)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na <= 0 or nb <= 0:
        raise ValueError("Procrustes comparison needs nonzero configuration variation")
    a, b = a / na, b / nb
    u, _, vt = np.linalg.svd(b.T @ a)
    return float(np.linalg.norm(a - b @ (u @ vt)))


def leave_one_out(counts: np.ndarray, keep: np.ndarray, samples: list[str],
                  full_scores: dict[str, np.ndarray], reference_depth: float) -> pd.DataFrame:
    matrix = validate_matrix(counts)
    keep = np.asarray(keep, dtype=bool)
    if (matrix.shape[1] < 4 or len(samples) != matrix.shape[1]
            or keep.shape != (matrix.shape[0],) or keep.sum() < 2):
        raise ValueError("Leave-one-out study needs >=4 samples and >=2 fixed retained genes")
    rows = []
    for omitted, name in enumerate(samples):
        included = np.arange(matrix.shape[1]) != omitted
        subset = matrix[:, included]
        factors, eligible = size_factors(subset)
        for method in METHODS:
            transformed = transformed_counts(subset, factors[method], reference_depth)[keep]
            scores, ratios = principal_components(transformed.T)
            error = procrustes_residual(full_scores[method][included], scores)
            rows.append({"omitted_sample": name, "method": method, "remaining_samples": int(included.sum()),
                         "fixed_pca_genes": int(keep.sum()), "positive_reference_genes": eligible,
                         "pc1_variance_fraction": float(ratios[0]), "pc2_variance_fraction": float(ratios[1]),
                         "procrustes_residual": error})
    return pd.DataFrame(rows)


def _pca_comparison(scores: dict, ratios: dict, metadata: pd.DataFrame, path: Path) -> None:
    colors = {"treated": "#b34d20", "untreated": "#1b6a8a"}
    shapes = {"single-read": "o", "paired-end": "s"}
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.7), layout="constrained")
    offsets = {"untreated1": (-7, -18), "untreated2": (7, 8), "untreated3": (7, 8),
               "untreated4": (7, -18), "treated1": (7, 8), "treated2": (7, -18), "treated3": (7, 8)}
    for ax, method in zip(axes, METHODS):
        for i, (sample, row) in enumerate(metadata.iterrows()):
            ax.scatter(*scores[method][i], color=colors[row.condition], marker=shapes[row.type], s=65)
            offset = offsets.get(sample, (7, 8))
            ax.annotate(sample, scores[method][i], xytext=offset, textcoords="offset points", fontsize=8,
                        ha="right" if offset[0] < 0 else "left")
        ax.set(title=METHOD_LABELS[method], xlabel=f"PC1 ({ratios[method][0]:.1%})",
               ylabel=f"PC2 ({ratios[method][1]:.1%})")
        ax.margins(.35)
        ax.spines[["top", "right"]].set_visible(False)
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=c, label=l)
               for l, c in colors.items()]
    handles += [Line2D([0], [0], marker=m, linestyle="none", color="#666666", label=l)
                for l, m in shapes.items()]
    axes[1].legend(handles=handles, loc="upper left", fontsize=9)
    fig.suptitle("Normalization sensitivity: same genes, reference depth and +1 offset")
    fig.savefig(path, dpi=160, metadata={"Software": "pasilla-expression-qc"})
    plt.close(fig)


def run_sensitivity(counts: pd.DataFrame, keep: np.ndarray, metadata: pd.DataFrame, out: Path) -> dict:
    matrix = validate_matrix(counts.to_numpy())
    factors, eligible = size_factors(matrix)
    reference_depth = geometric_mean(matrix.sum(axis=0))
    samples = list(counts.columns)
    normalization = metadata.copy()
    normalization.index.name = "sample"
    normalization["total_count_factor"] = factors["total_count"]
    normalization["median_ratio_factor"] = factors["median_ratio"]
    normalization["median_to_total_factor_ratio"] = factors["median_ratio"] / factors["total_count"]
    normalization["largest_gene_fraction"] = matrix.max(axis=0) / matrix.sum(axis=0)
    normalization.to_csv(out / "normalization_factors.csv", float_format="%.9f")
    scores, ratios, distances = {}, {}, {}
    for method in METHODS:
        transformed = transformed_counts(matrix, factors[method], reference_depth)[keep]
        scores[method], ratios[method] = principal_components(transformed.T)
        distances[method], correlation = pairwise_summaries(transformed.T)
        pc = metadata.copy()
        pc.index.name = "sample"
        pc["PC1"], pc["PC2"] = scores[method][:, 0], scores[method][:, 1]
        pc.to_csv(out / f"pca_{method}.csv", float_format="%.9f")
        for prefix, values in [("distances", distances[method]), ("correlations", correlation)]:
            pd.DataFrame(values, index=samples, columns=samples).rename_axis("sample").to_csv(
                out / f"{prefix}_{method}.csv", float_format="%.9f")
    stability = leave_one_out(matrix, keep, samples, scores, reference_depth)
    stability.to_csv(out / "leave_one_out_stability.csv", index=False, float_format="%.9f")
    _pca_comparison(scores, ratios, metadata, out / "normalization_comparison.png")
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), layout="constrained")
    vmax = max(float(d.max()) for d in distances.values())
    for ax, method in zip(axes, METHODS):
        img = ax.imshow(distances[method], cmap="viridis", vmin=0, vmax=vmax)
        ax.set(xticks=range(len(samples)), yticks=range(len(samples)), xticklabels=samples,
               yticklabels=samples, title=METHOD_LABELS[method])
        ax.tick_params(axis="x", labelrotation=45, labelsize=9)
        ax.tick_params(axis="y", labelsize=9)
        for label in ax.get_xticklabels():
            label.set_ha("right")
    fig.colorbar(img, ax=axes, shrink=.85, label="RMS difference across retained log2 profiles")
    fig.suptitle("Sample distances at a common reference count scale")
    fig.savefig(out / "normalization_distances.png", dpi=160, metadata={"Software": "pasilla-expression-qc"})
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 4.8), layout="constrained")
    for delta, method, color in [(-.12, "total_count", "#1b6a8a"), (.12, "median_ratio", "#b34d20")]:
        values = stability[stability.method == method].set_index("omitted_sample").loc[samples]
        ax.scatter(np.arange(len(samples)) + delta, values.procrustes_residual, label=METHOD_LABELS[method], color=color)
    ax.set(xticks=np.arange(len(samples)), xticklabels=samples, xlabel="Omitted sample",
           ylabel="Normalized Procrustes residual", title="Leave-one-sample-out PCA sensitivity")
    ax.tick_params(axis="x", labelrotation=25)
    ax.set_ylim(bottom=0)
    ax.legend(loc="best")
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(out / "leave_one_out_stability.png", dpi=160, metadata={"Software": "pasilla-expression-qc"})
    plt.close(fig)
    errors = {}
    for method in METHODS:
        values = stability[stability.method == method]
        worst = values.loc[values.procrustes_residual.idxmax()]
        errors[method] = {"minimum": float(values.procrustes_residual.min()),
                          "maximum": float(worst.procrustes_residual), "largest_change_omitted_sample": worst.omitted_sample}
    return {"reference_depth": reference_depth, "positive_reference_genes": eligible,
            "factor_ratio_range": [float(normalization.median_to_total_factor_ratio.min()),
                                   float(normalization.median_to_total_factor_ratio.max())],
            "pca_variance_fractions": {k: v.tolist() for k, v in ratios.items()},
            "between_method_procrustes_residual": procrustes_residual(scores["total_count"], scores["median_ratio"]),
            "leave_one_out_residuals": errors,
            "leave_one_out_policy": "Fixed full-data PCA gene set and reference depth; re-estimate factors on remaining samples"}
