"""Inspect the declared sample design without fitting an expression model."""
from __future__ import annotations

import numpy as np
import pandas as pd

CONDITIONS = ("untreated", "treated")
READ_TYPES = ("single-read", "paired-end")


def inspect_design(metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Return cell counts, a reference-coded additive matrix, and diagnostics.

    The intended levels stay fixed even when a subset omits a level. This makes
    unidentifiable columns visible rather than silently dropping the problem.
    """
    if metadata.empty or not {"condition", "type"}.issubset(metadata.columns):
        raise ValueError("Design inspection needs sample condition and type")
    if metadata.index.has_duplicates:
        raise ValueError("Design inspection needs unique sample identifiers")
    if metadata[["condition", "type"]].isna().any().any():
        raise ValueError("Design annotations must not be missing")
    if (not set(metadata.condition).issubset(CONDITIONS)
            or not set(metadata.type).issubset(READ_TYPES)):
        raise ValueError("Design inspection expects the documented pasilla condition and read-type levels")
    cross = pd.crosstab(metadata.condition, metadata.type).reindex(
        index=CONDITIONS, columns=READ_TYPES, fill_value=0)
    cross.index.name, cross.columns.name = "condition", "read_type"
    matrix = pd.DataFrame({
        "intercept": np.ones(len(metadata), dtype=int),
        "condition_treated": (metadata.condition == "treated").astype(int),
        "type_paired_end": (metadata.type == "paired-end").astype(int),
    }, index=metadata.index)
    matrix.index.name = "sample"
    rank = int(np.linalg.matrix_rank(matrix.to_numpy(dtype=float)))
    missing = {
        "condition": [v for v in CONDITIONS if v not in set(metadata.condition)],
        "read_type": [v for v in READ_TYPES if v not in set(metadata.type)],
    }
    all_levels = not any(missing.values())
    full_rank = rank == len(matrix.columns)
    # For two observed levels in each of two binary factors, a deficient
    # intercept-plus-indicators matrix means their effects cannot be separated.
    perfect = bool(all_levels and not full_rank)
    if not all_levels:
        assessment = "Not assessable for both intended contrasts: at least one factor level is absent"
    elif perfect:
        assessment = "Perfectly confounded: condition and read type cannot be separated in this additive design"
    else:
        assessment = "Not perfectly confounded: the additive design columns are linearly independent"
    warnings = ["Matrix rank is an algebraic check, not evidence for a treatment effect or adequate statistical power."]
    if not full_rank:
        warnings.append("Do not interpret separate condition and read-type effects from this rank-deficient design.")
    if int(cross.to_numpy().min()) < 2:
        warnings.append("At least one condition/read-type cell contains fewer than two samples, including any empty cell.")
    if len(metadata) < 10:
        warnings.append("This is a small sample set; model assumptions and uncertainty require separate evaluation.")
    summary = {
        "formula": "intercept + indicator(condition=treated) + indicator(type=paired-end)",
        "reference_levels": {"condition": "untreated", "read_type": "single-read"},
        "samples": int(len(metadata)), "columns": list(matrix.columns),
        "rank": rank, "column_count": int(len(matrix.columns)),
        "full_column_rank": bool(full_rank),
        "residual_dimension": int(len(metadata) - rank),
        "all_intended_levels_observed": bool(all_levels), "missing_levels": missing,
        "perfect_confounding": perfect if all_levels else None,
        "assessment": assessment,
        "all_design_cells_observed": bool((cross.to_numpy() > 0).all()),
        "equal_cell_counts": bool(np.unique(cross.to_numpy()).size == 1),
        "minimum_cell_size": int(cross.to_numpy().min()),
        "cell_counts": {condition: {read_type: int(cross.loc[condition, read_type])
                                    for read_type in READ_TYPES} for condition in CONDITIONS},
        "warnings": warnings,
        "model_fitted": False,
    }
    return cross, matrix, summary
