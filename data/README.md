# Data provenance and terms

These are public, derived Drosophila gene-count data, not patient records or synthetic fixtures.

The bundled `pasilla_gene_counts.tsv`, `pasilla_sample_annotation.csv`, and `UPSTREAM_DESCRIPTION` were retrieved unchanged on 2026-10-20 from the Bioconductor source mirror at commit `3c91e8df8add3549129e46b7a7a078fda4e2755a`. The preserved description identifies that development snapshot as package version 1.41.0. The release-page version may differ; the snapshot and hashes, rather than a moving release page, define these inputs.

See [the pinned source](https://github.com/bioconductor-source/pasilla/tree/3c91e8df8add3549129e46b7a7a078fda4e2755a), the [official package page](https://bioconductor.org/packages/release/data/experiment/html/pasilla.html), and `manifest.json` for the exact per-file URLs, byte counts, and SHA-256 values. The counts have 14,599 gene rows and seven sample columns. The annotation is left in its original row order and spelling; only the in-memory analysis maps its `fb` suffix to count-column names.

The data originate from Brooks AN, Yang L, Duff MO, Hansen KD, Park JW, Dudoit S, Brenner SE, Graveley BR. (2011). *Conservation of an RNA regulatory map between Drosophila and mammals.* Genome Research, 21, 193–202. DOI: [10.1101/gr.108662.110](https://doi.org/10.1101/gr.108662.110). Wolfgang Huber and Alejandro Reyes provide the derived `pasilla` package, DOI: [10.18129/B9.bioc.pasilla](https://doi.org/10.18129/B9.bioc.pasilla).

The preserved upstream metadata declares `License: LGPL` without a version qualifier. These files retain those upstream terms and attribution; this repository's MIT code license does not relicense them. Copies of LGPL version 3 and its incorporated GPL version 3 are provided as license reference texts. Their presence does not assert an additional version restriction absent from the upstream declaration. Consult the upstream package for authoritative licensing and any later changes.

The analysis output reports transformations of these inputs, not an independent reproduction of read alignment, gene annotation, or the original experiment.
