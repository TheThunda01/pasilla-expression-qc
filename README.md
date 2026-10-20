# Pasilla expression QC

**[Start here: project background, purpose, and key themes](docs/project-context.md).**

For a deeper explanation, read **[the biology and analysis concepts](docs/biology-and-concepts.md)**.

Before interpreting an expression experiment, I want to know whether the sample labels line up, how much the count totals vary, and whether a simple transformation changes the picture. This small project makes those checks visible using public Drosophila RNA-seq counts.

**Status: working educational analysis, version 0.3.0.** The inputs are real public data; this is a descriptive reanalysis, not a new biological finding or a reproduction of the original splicing study. It runs locally without downloading data after the Python dependencies are installed.

## Why this project

A useful expression analysis starts with knowing which samples and transformations are being compared. This project makes those choices inspectable before any inferential model is considered.

The `pasilla` Bioconductor package provides counts derived from the experiment described by Brooks and colleagues. Its annotation has three treated and four untreated samples, including single-read and paired-end libraries. The count table and original annotation are bundled unchanged, together with the upstream package description and a SHA-256 manifest. [Provenance and license details](data/README.md) explain the pinned source snapshot.

The program verifies those files before running. It is intentionally specific to this example: arbitrary datasets need a separate manifest, an appropriate metadata mapping, and a scientific review of the analysis choices.

## What it does

- Validate the count matrix and align the sample annotation by identifier rather than row order.
- Compare count-table totals, zero fractions, and numbers of detected genes across samples.
- Inspect condition-by-read-type cell counts and the rank of the explicitly coded additive sample design before interpreting the plots.
- View sample PCA on log2(CPM + 1), with condition shown by color and read type by shape.
- Compare that view with PCA on log2(count + 1), using the same genes in both analyses.
- Compare total-count and median-ratio normalization with the same reference depth and logarithm offset, then inspect the resulting sample distances and correlations.
- Refit both normalized PCA views after leaving out each sample, and quantify how much the remaining sample geometry changes.
- Inspect numeric outputs, input fingerprints, and the exact analysis settings.

Python data validation, reproducible scientific scripting, independent implementation of normalization formulas, matrix transformations, SVD-based PCA, geometric sensitivity analysis, clear plots, provenance checks, and tests for input failures. It does not run DESeq2, estimate gene-level differential expression, or test alternative splicing. [Methods and limitations](docs/methods.md) give the formulas and explain the boundaries. [Source notes](docs/source-notes.md) distinguish the original literature, current software documentation, and this implementation.

The separate `audit_results.py` command now checks completeness, source fingerprints, sample labels, and byte-level consistency of a result package. It does not rerun the statistics or resolve source-label uncertainty.

## Quick start

Use Python 3.12. From this repository's root:

```sh
python -m venv .venv
# macOS/Linux:
source .venv/bin/activate
# Windows PowerShell instead: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python analysis.py
```

The last command writes `results/`. To keep the included example outputs unchanged, use `python analysis.py --out my-results`. Output folders should be separate from the input folder. Run checks in CI with the included GitHub Actions workflow; a workflow file is not evidence that hosted CI has run.

Check the committed result package with:

```sh
python audit_results.py verify --manifest audits/example-artifacts.json
```

## Example results

[Read the example report](results/report.md), or open the [sample QC table](results/sample_qc.csv) and [machine-readable summary](results/summary.json).

![Normalization comparison with common reference scale](results/normalization_comparison.png)

The analysis begins with 14,599 genes and seven samples. It retains 8,148 genes with at least 10 counts in at least three samples for every PCA view. Count-table totals and zero fractions use all genes. These are declared exploratory settings, not universal thresholds for accepting an experiment.

The new normalization comparison estimates relative scale in two ways: from each sample's total counts, and from median log ratios across the 9,063 genes positive in every sample. Both use a reference depth of about 12.46 million counts before applying `log2(normalized count + 1)`. That common scale matters: changing the overall count scale while keeping a `+1` offset can itself change low-count relationships.

In the included run, median-ratio/total-count size-factor ratios range from **0.947 to 1.086**. The two-dimensional PCA configurations have a normalized Procrustes residual of **0.0548** after aligning their position, scale, rotation and reflection. These numbers describe sensitivity; they do not identify a winning normalization method.

![Leave-one-sample-out PCA sensitivity](results/leave_one_out_stability.png)

Omitting `treated1` produces the largest change in the two-dimensional configuration under both normalizations (residuals **0.2330** and **0.2187**, respectively). That is a reason to understand its influence, not a reason to discard it. With only seven samples, this exercise cannot distinguish experimental heterogeneity from a sample-quality problem. See the [full report](results/report.md) for every omission and the [methods](docs/methods.md) for what is held fixed.

The PCA pictures are summaries of the chosen transformed matrices. A cluster is not a significance test, and separation does not tell us whether condition, library preparation, or another factor caused it. Read type differs among these samples, so the plots display that annotation rather than concealing it. No sample is removed because of where it falls on a plot.

## Practical use and limitations

Use the [practical-use guide](docs/practical-use.md) to create and verify a snapshot for a new result directory. A byte match is distinct from numerical reproducibility and biological validity. The upstream article/package sample-allocation discrepancy remains unresolved, and this project retains the package labels. No differential-expression or splicing inference is added.

## Documentation

| Guide | Purpose |
| --- | --- |
| [Project context](docs/project-context.md) | Background, rationale, and key themes |
| [Biology and concepts](docs/biology-and-concepts.md) | Counts, normalization, PCA, and interpretation |
| [Methods](docs/methods.md) | Exact calculations and design checks |
| [References and source notes](docs/source-notes.md) | Scientific and software attribution |
| [Data provenance](data/README.md) | Pinned files, licenses, and unresolved labels |
| [Practical use](docs/practical-use.md) | Contracts, worked audit, and next improvements |
| [Literature review](docs/literature-review.md) | Bounded evidence review and implementation rationale |
| [Search log](docs/search-log.json) | Exact queries, counts, identifiers, and decisions |

## Development and reuse

Run `python -m unittest discover -s tests -v`. The scientific example remains version 0.3.0; the independent audit has schema version 1.0. Existing numeric and plot outputs are retained. See [contribution notes](CONTRIBUTING.md).

### References

1. Huber W, Reyes A. [pasilla data package](https://bioconductor.org/packages/release/data/experiment/html/pasilla.html). DOI: [10.18129/B9.bioc.pasilla](https://doi.org/10.18129/B9.bioc.pasilla).
2. Brooks AN, et al. (2011). Conservation of an RNA regulatory map between Drosophila and mammals. *Genome Research*, 21, 193–202. DOI: [10.1101/gr.108662.110](https://doi.org/10.1101/gr.108662.110).
3. Love MI, Huber W, Anders S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550. DOI: [10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8). Background on count-data analysis; its inferential method is not used here.
4. Anders S, Huber W. (2010). Differential expression analysis for sequence count data. *Genome Biology*, 11, R106. DOI: [10.1186/gb-2010-11-10-r106](https://doi.org/10.1186/gb-2010-11-10-r106). Source of the median-ratio normalization approach.
5. Robinson MD, Oshlack A. (2010). A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biology*, 11, R25. DOI: [10.1186/gb-2010-11-3-r25](https://doi.org/10.1186/gb-2010-11-3-r25). Explains composition effects; its TMM method is not implemented here.
6. Schönemann PH. (1966). A generalized solution of the orthogonal Procrustes problem. *Psychometrika*, 31, 1–10. DOI: [10.1007/BF02289451](https://doi.org/10.1007/BF02289451).

Original analysis code and documentation are MIT-licensed. The upstream data retain their original terms; see `data/README.md`.
