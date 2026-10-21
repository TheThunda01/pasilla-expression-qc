# Pasilla expression QC: the project in context

## Background and rationale

An RNA-sequencing (RNA-seq) count table can look tidy while still leaving important questions unanswered. Do the sample labels match? How different are the library totals? Does the picture change after normalization, or when one sample is omitted? This project explores those questions before attempting formal biological inference.

The example uses published Drosophila data from the [pasilla Bioconductor package](https://doi.org/10.18129/B9.bioc.pasilla), derived from [Brooks et al. (2011)](https://doi.org/10.1101/gr.108662.110). The supplied table contains 14,599 genes and seven samples, with treated and untreated conditions and two sequencing read types. Reusing a documented public example makes the analysis reproducible and lets the project focus on its methodological choices.

## What it does

The workflow checks input fingerprints, validates integer counts, aligns metadata by identifier, and summarizes count totals and zero fractions. It then compares sample principal component analysis (PCA), distances, and correlations under total-count and median-ratio normalization. The comparison uses a common reference depth and logarithm offset so an arbitrary change of scale does not obscure the normalization question.

A second analysis omits each sample in turn and refits normalization and PCA. It measures changes in the remaining samples' two-dimensional geometry after aligning translation, scale, rotation, and reflection. The original PCA gene set stays fixed, so this experiment examines sample inclusion rather than repeating feature selection.

## Why it matters

Normalization choices rely on assumptions. [Anders and Huber (2010)](https://doi.org/10.1186/gb-2010-11-10-r106) and [Robinson and Oshlack (2010)](https://doi.org/10.1186/gb-2010-11-3-r25) explain why library totals and expression composition can complicate comparisons. Here, disagreement between views is something to investigate; agreement does not prove that either method is correct.

The project demonstrates how to make analysis decisions inspectable through code, tests, provenance, numerical outputs, and plots. It does not identify differentially expressed genes, reproduce the original splicing study, or establish a treatment mechanism. An influential sample is not automatically a bad sample.

## Key themes

- **Reproducibility:** pinned public inputs and executable checks.
- **Normalization assumptions:** compare reasonable choices without declaring a winner.
- **Sensitivity:** examine how one sample changes a small dataset's visual summary.
- **Careful interpretation:** separate descriptive patterns from biological evidence.

The [glossary](glossary.md) defines the biological and numerical terms, and [assumptions](assumptions.md) distinguishes rationale from evidence. The [example report](../results/report.md) shows the results; the [methods](methods.md) explain how to read them.
