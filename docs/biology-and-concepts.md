# Biology and concepts behind the example

## What biological question led to these data?

Pasilla is an RNA-binding protein encoded by the Drosophila melanogaster gene `ps`. RNA-binding proteins can influence how an RNA molecule is processed after it is transcribed. [Brooks et al. (2011)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3032923/) investigated Pasilla-related splicing regulation in fly-derived S2 cells and compared aspects of that regulation with mammalian NOVA proteins. These are cultured model-organism cells, not human samples or whole animals observed in a clinical study. The [NCBI gene record](https://www.ncbi.nlm.nih.gov/gene/44258) is retained as a background annotation link; it could not be reopened during this review. The general RNA-regulation context is supported by the source study, and an annotation is not a new result generated here.

This project starts much later in the evidence chain: with a public teaching count table and its annotation. It asks whether the inputs are coherent and how analysis choices affect sample summaries. It does not repeat the experiment or reconstruct the original biological conclusions.

## From transcription to a table of counts

Transcription produces RNA from a DNA template. The amount of an RNA present at a particular time reflects both production and subsequent processing or removal. A measured RNA abundance therefore is not a direct measurement of transcription rate, and RNA abundance need not track protein abundance one-for-one.

RNA sequencing (RNA-seq) samples material from a prepared library. After sequence reads are mapped and assigned to features, a count table records how many assignments were associated with each gene in each sample. Those values depend on both the biological material and the measurement workflow. A larger gene count can reflect more captured RNA, greater sequencing depth, or other technical differences. It should not automatically be read as an absolute number of molecules per cell. The [DESeq2 methods paper](https://doi.org/10.1186/s13059-014-0550-8) explains why count properties and sample variability matter for statistical analysis.

Here, each row is a gene identifier and each column is a library. Summing a column gives the supplied **count-table total**, which need not equal the original number of sequenced reads.

## What does knockdown mean?

The source experiment used RNA interference (RNAi), RNA-mediated reduction of a target's expression, to reduce Pasilla expression. A knockdown reduces the abundance or activity of a target; it is not equivalent to deleting the gene, and a response after perturbation does not by itself identify every direct molecular target. This repository contains no experimental knockdown procedure and does not estimate knockdown efficiency. It uses the public `treated` and `untreated` labels as supplied. Those labels describe the derived dataset used here, not an independently verified reconstruction of every experimental detail in the paper.

## Gene abundance and splicing are different questions

Many initial RNA transcripts contain exons separated by introns. Splicing removes introns and joins exons; alternative splicing can produce different transcript forms from the same gene. Total abundance and the relative use of an exon can change independently.

That distinction is important because the original study concerned splicing regulation, while this project's main input aggregates assignments at **gene level**. A gene-total table cannot establish which exon or splice junction changed, nor demonstrate an isoform switch. The [pasilla data-creation vignette](https://bioconductor.org/packages/release/data/experiment/vignettes/pasilla/inst/doc/create_objects.html) describes the broader derived package, which includes exon-oriented material. This workflow uses only the stated gene-count table and annotation; it performs no exon-usage or splice-junction analysis.

## Why normalization has assumptions

Total-count scaling compares samples after accounting for their column sums. Median-ratio scaling estimates a relative sample factor from many genes, reducing the influence of a small number of large counts. Neither produces absolute RNA measurements. Both depend on what changed across samples; a widespread directional shift can make a supposedly stable reference inappropriate. [Anders and Huber (2010)](https://doi.org/10.1186/gb-2010-11-10-r106) and [Robinson and Oshlack (2010)](https://doi.org/10.1186/gb-2010-11-3-r25) motivate these issues.

The project compares both methods at a shared count scale before applying $`\log_2(\text{normalized count}+1)`$. The offset prevents taking the logarithm of zero but affects low-count relationships. Agreement between methods is useful to inspect, not proof that their assumptions hold.

## Read type and condition in the sample design

The supplied annotation contains two untreated single-read libraries, two untreated paired-end libraries, one treated single-read library, and two treated paired-end libraries. A paired-end library reads both ends of a fragment; it does not provide two biological replicates. A biological replicate is a separate biological sample, not an extra read or sequencing lane. A technical replicate repeats a measurement of the same biological material.

There is a provenance caveat: the original article's Results text describes four treated and three untreated cultures, whereas this package annotation labels three treated and four untreated libraries. This QC exercise preserves the supplied labels and does not resolve that discrepancy. Verify the primary-data mapping before using these labels for treatment-effect conclusions.

The program turns the two annotations into an intercept, a treated indicator, and a paired-end indicator. Rank is the number of independent columns. These three columns have rank three: condition and read type are **not perfectly confounded** in this additive design. If every treated sample instead used one read type and every control the other, separate coded effects would not be identifiable. The current [DESeq2 design guidance](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) explains this distinction.

Full rank does not imply adequate power or absence of other confounding. There are only seven libraries, one cell contains a single sample, and interactions or unrecorded technical variation may matter. No gene-level model is fitted here.

## Reading the visual summaries and related tools

Principal component analysis (PCA) summarizes directions of variation in a chosen transformed matrix. Its axes are not automatically biological mechanisms. Removing one sample can change that geometry without showing that the sample was erroneous. The omission analysis is descriptive sensitivity analysis, not predictive validation.

[DESeq2's official repository](https://github.com/thelovelab/DESeq2) and [Bioconductor package](https://bioconductor.org/packages/DESeq2) provide a fuller, established modeling workflow. This project adapts mathematical ideas and reporting practices, not their source code or validation claims. The reviewed source revision is recorded in [source notes](source-notes.md). Consult [the results](../results/report.md) alongside [the exact methods](methods.md); no PCA plot here establishes differential expression, alternative splicing, or a causal mechanism.
