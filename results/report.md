# Pasilla expression QC: example run

This run uses 14,599 genes and 7 samples from the pinned public count table. 8,148 genes pass the exploratory filter (at least 10 counts in at least 3 samples).

The log2(CPM + 1) PCA assigns 40.4% of variance to PC1 and 32.4% to PC2. The log2(count + 1) sensitivity analysis assigns 79.6% and 10.6%, respectively. These percentages describe two different transformed data matrices; they are not an accuracy comparison.

![PCA after total-count normalization](pca_logcpm.png)

![PCA without total-count normalization](pca_logcounts.png)

The plots display both condition and read type so potential technical structure is visible. They do not isolate the cause of sample separation or establish differential expression. No sample is dropped based on these plots.

![Count-table totals](count_totals.png)

Count-table totals sum the supplied gene counts. They are not the original total sequenced-read counts. The full numeric results are in `sample_qc.csv`, `sample_correlations.csv`, and `summary.json`.
