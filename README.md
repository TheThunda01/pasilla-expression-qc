# Pasilla expression QC

Before interpreting an expression experiment, I want to know whether the sample labels line up, how much the count totals vary, and whether a simple transformation changes the picture. This small project makes those checks visible using public Drosophila RNA-seq counts.

**Status: completed educational analysis.** The inputs are real public data; this is a descriptive reanalysis, not a new biological finding or a reproduction of the original splicing study. It runs locally without downloading data after the Python dependencies are installed.

## What you can explore

- Validate the count matrix and align the sample annotation by identifier rather than row order.
- Compare count-table totals, zero fractions, and numbers of detected genes across samples.
- View sample PCA on log2(CPM + 1), with condition shown by color and read type by shape.
- Compare that view with PCA on log2(count + 1), using the same genes in both analyses.
- Inspect numeric outputs, input fingerprints, and the exact analysis settings.

## Run it

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

## Start with the results

[Read the example report](results/report.md), or open the [sample QC table](results/sample_qc.csv) and [machine-readable summary](results/summary.json).

![Sample PCA after CPM normalization](results/pca_logcpm.png)

The analysis begins with 14,599 genes and seven samples. It retains 8,148 genes with at least 10 counts in at least three samples for both PCA views. Count-table totals and zero fractions use all genes. These are declared exploratory settings, not universal thresholds for accepting an experiment.

The PCA pictures are summaries of the chosen transformed matrices. A cluster is not a significance test, and separation does not tell us whether condition, library preparation, or another factor caused it. Read type differs among these samples, so the plots display that annotation rather than concealing it. No sample is removed because of where it falls on a plot.

## Where the data came from

The `pasilla` Bioconductor package provides counts derived from the experiment described by Brooks and colleagues. Its annotation has three treated and four untreated samples, including single-read and paired-end libraries. The count table and original annotation are bundled unchanged, together with the upstream package description and a SHA-256 manifest. [Provenance and license details](data/README.md) explain the pinned source snapshot.

The program verifies those files before running. It is intentionally specific to this example: arbitrary datasets need a separate manifest, an appropriate metadata mapping, and a scientific review of the analysis choices.

## What the project demonstrates

Python data validation, reproducible scientific scripting, matrix transformations, SVD-based PCA, clear plots, provenance checks, and tests for input failures. It does not implement DESeq2, estimate gene-level differential expression, or test alternative splicing. [Methods and limitations](docs/methods.md) give the formulas and explain the boundaries.

## References

1. Huber W, Reyes A. [pasilla data package](https://bioconductor.org/packages/release/data/experiment/html/pasilla.html). DOI: [10.18129/B9.bioc.pasilla](https://doi.org/10.18129/B9.bioc.pasilla).
2. Brooks AN, et al. (2011). Conservation of an RNA regulatory map between Drosophila and mammals. *Genome Research*, 21, 193–202. DOI: [10.1101/gr.108662.110](https://doi.org/10.1101/gr.108662.110).
3. Love MI, Huber W, Anders S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550. DOI: [10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8). Background on count-data analysis; its inferential method is not used here.

Original analysis code and documentation are MIT-licensed. The upstream data retain their original terms; see `data/README.md`.
