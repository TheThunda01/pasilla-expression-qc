# Glossary

These definitions describe how the terms are used in this project. The [methods](methods.md) give the exact calculations, and the [assumptions](assumptions.md) explain their limits.

## Biology and experimental evidence

| Term | Meaning here |
| --- | --- |
| DNA, RNA, and protein | Deoxyribonucleic acid (DNA) stores sequence information; ribonucleic acid (RNA) is produced from DNA by transcription; proteins are molecules made by translating some RNAs. RNA abundance does not directly measure protein abundance. |
| Gene and transcript | A gene is a genomic feature; a transcript is an RNA product. A gene can produce multiple transcript forms. |
| Exon, intron, and splicing | Exons are retained segments of a processed transcript; introns are removed during splicing. Alternative splicing changes which segments or splice sites are used. Gene-total counts cannot resolve these choices. |
| Pasilla (`ps`) | The fly RNA-binding protein studied in the source experiment. An RNA-binding protein interacts with RNA and can influence its processing. |
| S2 cells | Cultured cells derived from *Drosophila melanogaster*, the fruit fly; they are not human or whole-animal samples. |
| NOVA | Mammalian RNA-binding proteins compared with Pasilla in the source study; the acronym is part of their established name. |
| RNA interference (RNAi) and knockdown | RNA-mediated reduction of a target's expression. Knockdown is a reduction, not necessarily complete removal or deletion of a gene. |
| RNA sequencing (RNA-seq) | Sequencing of material prepared from RNA to investigate the transcripts represented in a sample. |
| Library, read, and fragment | A library is prepared material for sequencing; a read is a reported sequence; a fragment is a piece of the prepared material. Count assignments are not necessarily original molecule counts. |
| Single-read and paired-end | Reading one end or both ends of a library fragment. Two mates from one fragment are not two biological replicates. |
| Biological and technical replicate | A biological replicate comes from a separately sampled biological unit; a technical replicate repeats measurement of the same material. Sequencing lanes do not by themselves establish biological independence. |
| Differential expression | A statistical comparison of gene abundance between conditions, accounting for a specified design and variability. This project does not perform it. |
| Causal mechanism | An explanation of how an intervention produces an effect. A cluster, correlation, or normalization difference alone does not establish one. |

## Mathematics and statistics

| Term | Meaning here |
| --- | --- |
| Count-table total and sequencing depth | The former is the sum of supplied gene assignments; the latter usually describes sequencing output. The quantities can differ because of processing and assignment rules. |
| CPM | Counts per million supplied assignments: a gene count divided by its unfiltered column total, multiplied by one million. |
| Normalization and size factor | Rescaling to compare samples. A size factor is the sample-specific divisor; normalization does not remove every technical effect. |
| Composition bias | A change in the relative abundance of some features changes the share of sequencing output available to others. |
| Median-ratio factor | A relative scale estimate based on each gene's ratio to a geometric reference, combined across genes. This implementation takes the exponential of the median log ratio. |
| TMM | Trimmed mean of M-values, a relative normalization method based on weighted averages of selected log expression ratios (M-values). It motivates the composition-bias discussion but is not implemented here. |
| Dispersion and negative-binomial model | Dispersion describes variation beyond a simple mean-dependent count model. A negative-binomial model permits count variance greater than its mean. These are part of the cited inferential methods, not models fitted here. |
| VST and rlog | Variance-stabilizing transformation and regularized logarithm are DESeq2 transformations designed for count-data visualization. This project's `log2(count + 1)` calculation does not reproduce them. |
| Geometric mean and logarithm | The geometric mean of positive values is the exponential of their mean natural logarithm. `log2` expresses values on a base-two logarithmic scale. |
| Pseudocount | The added value in `log2(normalized count + 1)`, which permits zeros but influences low-count differences. |
| PCA and principal component (PC) | Principal component analysis summarizes directions of variation. PC1 and PC2 are the first two directions; their coordinates depend on the input transformation. |
| Centering and standardization | Centering subtracts a mean. Standardization also divides by a standard deviation; this project's PCA centers genes but does not standardize them. |
| Loading and variance explained | A loading describes a feature's contribution to a PC. Variance explained is the component's share of variation in the analyzed matrix, not biological explanatory power. |
| Singular value decomposition (SVD) | A matrix factorization used to calculate PCA and the alignment between two coordinate configurations. |
| Root mean square (RMS) and Euclidean distance | RMS distance is the square root of the mean squared gene-wise difference. It equals Euclidean distance divided by the square root of the number of genes here. |
| Pearson correlation | A measure of linear association between two profiles. It is undefined when a profile has no variation and does not establish sample identity. |
| Design matrix, intercept, and indicator | Numeric encoding of annotations: a constant intercept plus zero/one columns for the selected condition and read-type categories. |
| Rank, estimability, and perfect confounding | Rank counts independent matrix columns. A separate coefficient is estimable only if the design provides independent information for it. Perfect confounding makes coded effects inseparable. |
| Additive design and interaction | An additive design encodes separate contributions without a condition-by-read-type interaction. An interaction would allow one contribution to depend on the other. No coefficients are fitted here. |
| Residual dimension | Number of samples minus design rank. It is a design diagnostic here, not an estimated residual variance or a power calculation. |
| Leave-one-sample-out sensitivity | Recalculating summaries after omitting each sample. This assesses influence on the display; it is not predictive cross-validation. |
| Orthogonal Procrustes alignment | Aligning two coordinate configurations by rotation or reflection after centering and scaling. The reported residual measures remaining shape difference. |
| Frobenius norm | Square root of the sum of squared matrix entries. Unit-Frobenius scaling removes overall configuration size before alignment. |
| p-value and multiple testing | A p-value measures incompatibility with a specified null model; multiple-testing methods address the many simultaneous gene comparisons in an inferential analysis. Neither is calculated here. |

## Files, software, and reproducibility

| Term | Meaning here |
| --- | --- |
| Metadata and annotation | Descriptions attached to samples or features, such as condition and read type. Matching these descriptions to identifiers does not establish that the descriptions are biologically correct. |
| TSV, CSV, and JSON | Tab-separated tables, comma-separated tables, and structured text objects, respectively. A schema defines expected fields and types. |
| Manifest, provenance, and commit | A manifest inventories files; provenance records where data and settings came from; a Git commit identifies a recorded repository revision. |
| SHA-256 and fingerprint | A hash algorithm and its byte-dependent identifier. A matching hash supports equality to the chosen bytes, not scientific correctness or authenticated ownership. |
| Floating point and tolerance | Finite-precision numerical representation and an allowed difference when comparing results. `rtol` is a relative allowance; `atol` is an absolute floor. |
| Nonfinite value | NaN (not a number) or an infinity. These values are rejected in numerical result comparisons. |
| CLI and exit status | A command-line interface and its integer completion code. Check the documented code before relying on generated files. |
| CI and dependency | Continuous integration runs checks automatically; a dependency is software required by the project. Local tests do not prove a hosted CI run occurred. |
| LGPL and MIT | Software-license names. The project code and upstream data have separate notices; consult the preserved license material rather than assuming one license covers both. |

