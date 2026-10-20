# Rapid literature and practice review

## Question and scope

What small operational check would make a descriptive public RNA-seq reanalysis easier to transfer and inspect without adding unsupported inference? This bounded review focused on sample metadata, provenance, and reproducible result packaging. It did not reassess normalization performance across experiments or conduct a complete systematic review.

### Round 2 update — 2026-10-20

This same-day follow-up asks how to compare numerical reruns while keeping sample identities and method settings exact. It is an additional implementation scope, not a new surveillance date or a repeat of the broad literature review.

## Search record

Search date: **2026-10-20**. The [machine-readable log](search-log.json) records exact PubMed queries, date filters, sort order, returned identifiers, and screening decisions. The initial recent query returned 37 records; only its first 10, ordered by publication date, were screened. A narrower metadata/provenance query returned six records, all screened. Those pages contain 14 distinct records. Two directly relevant software papers were retained, with an older source-study citation and official documentation added for context. No inference about the 27 unreturned broad-search records is made.

A quoted-title lookup for the foundational paper returned zero; the known DOI was then resolved through the PubMed identifier service and its record retrieved. The seeker paper was read through its complete rights-eligible PMC body. SampleExplorer and the source experiment were appraised from citation/abstract metadata in this round; the original source discrepancy also remains documented in the existing data notes. Supplementary official documentation and repository files were inspected directly. No protected abstract is republished here.

### Round 2 update — 2026-10-20

A narrow PubMed query for computational biology/bioinformatics, floating point, and reproducibility returned **zero records** in the 2022-01-01 to 2026-10-20 publication window. This is not evidence that the subject lacks literature. The engineering decision was instead checked against official Python 3.12 documentation. Exact syntax and the zero-result record are preserved in [search-update-02.json](search-update-02.json); the original search log is unchanged.

## What the sources support

| Source | Type and access | Relevant finding | Limit | Consequence here |
| --- | --- | --- | --- | --- |
| [Schoenbachler and Hughey, 2022, PMID 36389425](https://pubmed.ncbi.nlm.nih.gov/36389425/) | Primary software paper; complete PMC body read | seeker records parameters, session information, dependencies, and produces sample metadata with expression outputs. | Its end-to-end R workflow is not this Python QC script; recording provenance does not establish biological correctness. | Keep the existing analysis/environment record and add a transferable byte-level artifact snapshot. |
| [Chin and Lassmann, PMID 39786428](https://pubmed.ncbi.nlm.nih.gov/39786428/) | Primary software paper; abstract-only appraisal; online 2024, volume 41 issue 1 | SampleExplorer uses transcriptome metadata to retrieve relevant samples. | Retrieval relevance is not verification of sample labels; the underlying model and benchmark were not re-evaluated. | Keep source labels visible and compare output annotations with the chosen source, without treating agreement as truth. |
| [Brooks et al., 2011, PMID 20921232](https://pubmed.ncbi.nlm.nih.gov/20921232/) | Primary experimental source; abstract reviewed in this round | The original work concerns Pasilla-dependent splicing. | The available gene-count QC does not reproduce exon-level inference or resolve sample allocation by itself. | Retain descriptive scope and the unresolved article/package discrepancy. |
| [pasilla data-creation vignette, version 1.40.0, dated 5 May 2026](https://bioconductor.org/packages/release/data/experiment/vignettes/pasilla/inst/doc/create_objects.html) | Official source documentation; sample table inspected | Documents the derived package sample annotations and count-generation context. | A derived package is not independent evidence that all original labels are correct. | Audit sample names/order against the exact pinned annotation, not a silently edited table. |
| [seeker source, commit ec636b2d55ffa93265cfb3150df5d9295c35f07c](https://github.com/hugheylab/seeker/tree/ec636b2d55ffa93265cfb3150df5d9295c35f07c) | Official repository; README and DESCRIPTION inspected; version 1.1.6; MIT | Declares a standard interface and detailed logging for transcriptome processing. | No seeker execution or code-equivalence test was performed. | Adopt the general practice of explicit records; no code was copied and no new dependency was added. |

The implementation choice is an engineering judgment informed by those sources: preserving associations between metadata and outputs is useful, while a larger pipeline would not resolve the existing scientific uncertainty. None of the papers validates this new audit tool.

### Round 2 update — 2026-10-20

| Source | Type/access | Actual support | Limit and implementation consequence |
| --- | --- | --- | --- |
| [Python 3.12 `math.isclose`](https://docs.python.org/3.12/library/math.html#math.isclose), inspected documentation version 3.12.15 | Official API documentation | Defines comparison using relative and absolute tolerances; nonfinite values require explicit treatment. | Does not prescribe a scientific tolerance. The project declares `rtol=1e-7`, `atol=1e-8`, rejects nonfinite inputs, and tests allowed roundoff versus larger changes. |
| [Python 3.12 JSON documentation](https://docs.python.org/3.12/library/json.html) | Official API documentation | Describes parsing hooks and default handling of constants and repeated names. | Default parsing is not a scientific schema check. The comparator rejects duplicate keys and nonfinite numbers and compares key/type structure explicitly. |

These are original standard-library checks; no external implementation was copied. The tolerance is an engineering boundary for recorded numerical outputs, not a measure of experimental uncertainty.

## Changes made

Added `audit_results.py`, six tests, and an example manifest for the existing 23-file package. The audit rejects missing/empty artifacts, changed source fingerprints, incompatible sample order or labels, changed snapshot bytes, and accidental manifest overwrite. The README now uses consistent navigation, with an explicit [practical-use contract](practical-use.md). The original scientific code, source hashes, and checked-in numeric and plot results remain unchanged.

### Round 2 update — 2026-10-20

Added `compare_results.py`, eleven targeted tests, a numerical comparison example, and a CI rerun comparison. It compares 14 CSVs and two JSON files, preserves exact identities and integer counts, records environment-version differences separately, and reports the location of mismatches. The earlier artifact audit remains a separate byte-integrity tool. Original scientific outputs are unchanged.

## Remaining work

The next useful additions are a tolerance-based numerical comparator and more matrix-level consistency tests. Definitive source-label reconciliation remains separate evidence-gathering work. The audit should not be expanded into a claim of differential-expression readiness while that issue remains unresolved.

### Round 2 update — 2026-10-20

The numerical comparator proposed in round 1 is now implemented. Matrix symmetry/diagonal checks, a second operating-system comparison, and definitive source-label reconciliation remain open; current acceptance criteria are in [practical use](practical-use.md).

## Limits of this review

One reviewer screened bounded PubMed pages and a small set of official resources. The search was not registered, duplicated, risk-of-bias scored, or exhaustive across databases. The initial query was intentionally narrowed after mostly off-topic returns. Abstract-only access limited appraisal of two publications; repository inspection was limited to the named files. This review supports a modest workflow improvement, not a comparative claim that this tool outperforms existing RNA-seq systems.

### Round 2 update — 2026-10-20

No new biomedical publication was selected from the zero-result query. Official documentation supports API semantics, not empirical validation of the chosen tolerance. This update did not test a second operating system, redo the original experiment, or resolve the source-label discrepancy.

