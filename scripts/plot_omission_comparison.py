"""Compare the two saved leave-one-out residuals for each Pasilla sample."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from build_sample_explorer import METHODS, ROOT, load_results


def render(results: Path, output: Path) -> None:
    data, hashes = load_results(results)
    omissions = data["omissions"]
    order = sorted(data["samples"], key=lambda name: max(
        omissions[name, method]["procrustes_residual"] for method in METHODS), reverse=True)
    colors = ["#126681", "#aa4d20"]
    markers = ["o", "D"]
    with plt.rc_context({"svg.fonttype": "none", "svg.hashsalt": "pasilla-omission-comparison"}):
        fig, axis = plt.subplots(figsize=(8.6, 4.7), layout="constrained")
        for index, name in enumerate(order):
            values = [omissions[name, method]["procrustes_residual"] for method in METHODS]
            axis.hlines(index, min(values), max(values), color="#82929b", linewidth=2)
        for method, color, marker in zip(METHODS, colors, markers):
            values = [omissions[name, method]["procrustes_residual"] for name in order]
            circle = method == "total_count"
            axis.scatter(values, range(len(order)), label=METHODS[method], marker=marker,
                         facecolor="white" if circle else color, edgecolor=color,
                         linewidth=1.5 if circle else .6, s=88 if circle else 32, zorder=3)
            top_index = values.index(max(values))
            above = method == "total_count"
            axis.annotate(f"{values[top_index]:.3f}", (values[top_index], top_index),
                          xytext=(0, 9 if above else -10), textcoords="offset points",
                          ha="center", va="bottom" if above else "top", fontsize=10, color=color)
        axis.set(yticks=range(len(order)), yticklabels=order, xlabel="Normalized Procrustes residual",
                 ylabel="Sample omitted", ylim=(len(order) - .5, -.8))
        axis.set_xlim(left=0, right=max(r["procrustes_residual"] for r in omissions.values()) * 1.1 or 1)
        axis.grid(axis="x", alpha=.18)
        axis.spines[["top", "right", "left"]].set_visible(False)
        axis.legend(loc="lower right", frameon=False)
        axis.set_title("How much does the sample plot move after one library is omitted?",
                       loc="left", weight="bold", fontsize=14, pad=12)
        fig.supxlabel("Each point compares the six remaining samples with their positions in the full fit.\n"
                      "Larger values indicate a larger change after alignment.", fontsize=9)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, format="svg", metadata={"Date": None, "Description":
            "Source: results/leave_one_out_stability.csv. LF-normalized SHA-256: "
            + hashes["leave_one_out_stability.csv"]})
        plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=ROOT / "results")
    parser.add_argument("--output", type=Path, default=ROOT / "docs/figures/omission-comparison.svg")
    args = parser.parse_args()
    render(args.results, args.output)
