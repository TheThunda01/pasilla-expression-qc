# Methods and interpretation

Key vocabulary is defined as it arises below and in the [glossary](glossary.md); [assumptions](assumptions.md) records what is chosen, checked, or unresolved.

## The question

How do seven samples in a published Drosophila expression count table look before formal modeling? How sensitive are the sample summaries to normalization choice and to the inclusion of one sample?

This is an exploratory workflow with fixed choices. It does not estimate a treatment effect. The original experiment and its derived package are described in the [pasilla package](https://bioconductor.org/packages/release/data/experiment/html/pasilla.html) and its [data-creation vignette](https://bioconductor.org/packages/release/data/experiment/vignettes/pasilla/inst/doc/create_objects.html).

## Input validation and alignment

`analysis.py` checks every upstream file against the bundled SHA-256 and byte-count manifest. This detects changed input bytes; it does not independently establish the correctness of the original data.

The tab-separated-value (TSV) parser requires unique sample names and gene identifiers, a rectangular table, and nonnegative integer counts. Missing values, duplicate rows, fractions, negative values, and zero-total samples stop the run. It does not silently fill or deduplicate records.

The annotation's `file` field ends in `fb`. For this dataset only, removing that trailing suffix maps values such as `treated1fb` to the count column `treated1`. The program checks that the resulting names are unique and match the count columns exactly, then reorders metadata to the count-column order. It retains the original `file` name in the output. Condition and read-type labels must match the known annotation categories.

## Inspect the sample design

`design.py` produces a condition-by-read-type table and an explicitly reference-coded additive matrix: numeric columns encode separate condition and read-type indicators without their interaction. The reference categories are `untreated` and `single-read`; the three columns are a constant intercept, `condition_treated`, and `type_paired_end` (zero/one indicators). Intended levels remain fixed if a subset omits a category, making missing contrasts visible rather than silently discarding them.

Matrix rank uses NumPy's singular-value-based rank calculation. With both levels of both binary factors observed, rank deficiency indicates that the two annotations are perfectly confounded for this additive design. If a level is absent, the result instead reports that the intended contrasts cannot both be assessed. Rank three means the columns are independent, not that any effect is significant or that a model is scientifically adequate. The reported residual dimension is `number of samples - rank`; no residual variance or p-value is estimated.

For the bundled annotation, all four cells are observed, but the treated/single-read cell has only one sample. The matrix has rank three and residual dimension four. Condition and read type are not perfectly confounded; their allocation is nevertheless uneven. These checks mirror the design-identifiability concern in the [official DESeq2 guidance](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#model-matrix-not-full-rank). They do not adjust the PCA for read type or fit an expression model.

## Count summaries

For each sample, `count_table_total` is the sum of its supplied gene counts, `zero_fraction_all_genes` is the fraction of rows with zero counts, and `genes_with_positive_count` is the number of rows with a positive count. The total is not the original number of sequenced reads or an estimate of absolute RNA abundance. Gene counts also should not be interpreted as directly comparable absolute expression levels across genes of different lengths.

## Transformation and PCA

Every principal component analysis (PCA) view uses the same genes: at least 10 counts in at least three samples. This is a simple prespecified exploratory filter, not a test of differential expression or a comprehensive quality threshold.

For the first view, each count is divided by its sample's **unfiltered** count-table total, multiplied by one million, then transformed as $`\log_2(\mathrm{CPM}+1)`$, where CPM means counts per million supplied assignments. The second view uses $`\log_2(\text{count}+1)`$ without total-count normalization. Filtering therefore does not change the CPM denominator.

PCA treats samples as observations and retained genes as variables. It centers each gene across samples and uses singular value decomposition (SVD), a factorization of the matrix into orthogonal directions and their strengths. It does not divide genes by their standard deviations. The reported percentage for a component is its squared singular value divided by the sum of squared singular values. Component signs are made reproducible by a fixed loading-sign convention; their signs have no biological meaning.

The original `sample_correlations.csv` uses the $`\log_2(\mathrm{CPM}+1)`$ values for the retained genes. The normalization comparison also exports correlations and distances separately for both normalization methods. High sample correlation is not proof of correct sample identity or experimental quality.

Total-count normalization only rescales sample totals. It does not generally resolve composition bias, read-type differences, or unmeasured technical variation. The $`+1`$ pseudocount (the offset added to permit a logarithm at zero) and gene filter can influence low-count patterns. [Love et al. (2014)](https://doi.org/10.1186/s13059-014-0550-8) discuss the mean–variance behavior of RNA-seq counts and specialized transformations. This lightweight illustration is not a replacement for those methods.

## Compare normalization without changing the count scale

Let $`K_{ij}`$ be the original count for gene $`i`$ in sample $`j`$, $`L_j`$ the sum over **all** supplied genes, and $`D`$ the geometric mean of the seven library totals. A geometric mean is the exponential of the mean log value, so all contributing totals must be positive.

The total-count factor is $`L_j/D`$. For the median-ratio factor, retain only genes positive in all samples. For each such gene, calculate its mean log count, subtract that mean from each sample's log count, and take the median log ratio within each sample. Exponentiate these medians and rescale the resulting factors to have geometric mean one. The numerator and denominator of a gene ratio therefore refer to the same gene; no gene-length adjustment is attempted. No $`+1`$ is used to estimate the factors.

This follows the log-median convention documented in the current [DESeq2 implementation](source-notes.md). For an even number of reference genes, exponentiating the median of logs is not exactly the same operation as taking an arithmetic median of ratios. The chosen convention is explicit and tested on small examples. The program raises an error if every gene has a zero in at least one sample rather than silently substituting another estimator.

For each factor vector $`s`$, the comparable exploratory matrix is $`\log_2(K_{ij}/s_j+1)`$. Both vectors have geometric mean one, so both normalizations refer to the same overall count depth $`D`$ and the offset has the same count units. The older CPM plot remains available for continuity, but it has a different overall scale; comparing it directly with normalized-count PCA would mix factor differences with pseudocount-scale differences.

Total-count normalization can be sensitive to a minority of highly expressed genes. Median-ratio estimation instead relies on a sufficiently stable central reference across genes; widespread directional shifts can violate that assumption. The papers by [Anders and Huber (2010)](https://doi.org/10.1186/gb-2010-11-10-r106) and [Robinson and Oshlack (2010)](https://doi.org/10.1186/gb-2010-11-3-r25) motivate this concern. Agreement between the two methods is not independent evidence that either assumption holds, and this dataset has no external absolute-RNA reference that resolves the issue.

`normalization_factors.csv` contains both factors, their ratio, and the largest single-gene fraction of each original library total. The latter is a descriptive composition summary; it is not a contamination score or a quality threshold.

## Distances and correlations

For each method, the sample distance is the root mean square of the difference between two retained log-normalized gene profiles. It equals their Euclidean distance divided by the square root of the number of retained genes. The heatmaps share a color scale and retain the input sample order; no clustering algorithm or reordered dendrogram is implied.

Pearson correlations, which measure linear association between profiles, use the same transformed profiles. A constant profile makes Pearson correlation undefined, so the program stops rather than inventing a value. Both summaries depend on the retained genes and transformation, and neither establishes sample identity, mechanism, or causal response.

## Leave-one-sample-out PCA stability

The omission analysis measures how much the **displayed two-dimensional sample geometry** changes after omitting a single sample. For each of seven omissions and each normalization method:

1. Keep the full analysis's 8,148-gene PCA set fixed. This isolates sensitivity to sample inclusion; it does not measure sensitivity to repeating feature selection.
2. Re-estimate size factors from all genes in the six remaining samples. The set positive in every remaining sample can grow; its size is reported for each refit.
3. Keep the original full-data reference depth $`D`$ fixed. Technically, the centered factor vector is multiplied by the geometric mean of the remaining library totals divided by $`D`$. This anchors the log-transform offset to the same count scale. For total-count normalization it makes each retained sample's pre-PCA transformed profile invariant to which other sample is omitted.
4. Center genes again and refit the two-component PCA on the remaining six samples.
5. Take those same six samples' coordinates from the full-data PCA. Center both coordinate matrices, scale each so its Frobenius norm (square root of the sum of squared entries) is one, and use SVD to find the orthogonal rotation or reflection minimizing their difference.
6. Report the Frobenius norm of that aligned difference as `procrustes_residual`.

The alignment follows the orthogonal solution of [Schönemann (1966)](https://doi.org/10.1007/BF02289451). A residual near numerical zero means the two configurations have the same shape after the allowed transformations. Larger values mean a larger change by this definition. There is **no prespecified biological pass/fail cutoff**. Overall scale and orientation are deliberately discarded, so this statistic does not capture every possible change in expression variability. Variance explained by the two refitted PCs is exported separately.

This is neither predictive cross-validation nor a statistical test. An influential sample might reflect ordinary biological or technical heterogeneity, not an error. In particular, `treated1` is the only treated single-read sample. Its omission changes both the sample composition and the available contrast between condition and read type. The workflow reports that sensitivity without declaring the sample unsuitable or proving a condition effect. No sample is automatically removed.

## What can be concluded

The outputs describe the supplied tables under these analysis settings. Differences between the PCA views demonstrate sensitivity to a transformation choice; they do not establish which view is biologically correct. With seven samples, descriptive patterns need careful interpretation. Both condition and sequencing read type are displayed, but neither is statistically adjusted for here. Library-level scaling does not remove a design covariate, and the grouping of plotted points is not a fitted treatment effect.

This workflow produces no p-values, gene-level effect estimates, pathway claims, mechanism claims, or splicing conclusions. A future inferential analysis would require an explicit design, checks for estimability and confounding, suitable count modeling, and multiple-testing control. That work is a planned extension, not a completed result.

## Output dictionary

| File | Contents |
| --- | --- |
| `sample_qc.csv` | Original file name, condition, type, count total, zero fraction and positive-count genes, in count-column order. |
| `pca_logcpm.csv` | Per-sample PC1 and PC2 coordinates after $`\log_2(\mathrm{CPM}+1)`$. |
| `pca_logcounts.csv` | Per-sample PC1 and PC2 coordinates after $`\log_2(\text{count}+1)`$. |
| `sample_correlations.csv` | Sample Pearson correlations across retained $`\log_2(\mathrm{CPM}+1)`$ gene values. |
| `summary.json` | Dimensions, settings, variance ratios, input hashes and package versions. |
| `condition_read_type_counts.csv` | Sample counts in each intended condition/read-type combination. |
| `additive_design_matrix.csv` | The exact three-column reference-coded sample matrix. |
| `design_diagnostics.json` | Matrix rank, missing levels, residual dimension, and interpretation warnings; no fitted model. |
| `normalization_factors.csv` | Centered total-count and median-ratio factors, their ratio, and largest-gene fraction. |
| `pca_total_count.csv`, `pca_median_ratio.csv` | Coordinates from the comparison at a common count scale. |
| `distances_total_count.csv`, `distances_median_ratio.csv` | Root-mean-square sample distances over retained transformed genes. |
| `correlations_total_count.csv`, `correlations_median_ratio.csv` | Pearson correlations over the same profiles. |
| `leave_one_out_stability.csv` | Fourteen refits, reference-gene counts, retained-gene counts, variance fractions and Procrustes residuals. |
| `report.md` and PNGs | Human-readable run summary, original views, normalization comparison, distance heatmaps, and omission-sensitivity plot. |

## Failure modes and verification

The tests cover table shape, malformed counts, duplicate identifiers, sample alignment, file-integrity failures, centering, lack of per-gene standardization, and degenerate PCA inputs. Small synthetic matrices check known sample multipliers, an influential high-count gene, zero-count exclusion from the reference, geometric centering, fixed-depth invariance, hand-calculated distances, and the invariance of Procrustes alignment to translation, scale, rotation and reflection. Other tests confirm that a true shape change is detected and every sample omission is evaluated.

The GitHub Actions workflow runs the tests and the full public-data analysis. These checks assess the implemented numerical behavior and executable workflow rather than the scientific truth of the source dataset. Source verification, dataset checksums, and tested calculations serve different purposes; none substitutes for the others.
