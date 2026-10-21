# Practical use

## Where this fits

This is a small, inspectable check of a published count table before formal expression modeling. It is useful for learning how sample annotations, count totals, normalization (relative rescaling), and principal component analysis (PCA) sensitivity fit together. It does not settle the original article/package label discrepancy or test differential expression.

The artifact audit addresses a separate operational problem: a folder can accidentally contain tables and figures from different runs. It checks the expected 23 files, their exact bytes, the input fingerprints recorded in the summary, and the sample ordering and annotations in five sample-level tables. It does not recompute the analysis.

The [assumptions](assumptions.md) and [glossary](glossary.md) provide the scientific and technical context for these checks.

## Input and output contract

`analysis.py` continues to accept the pinned `data/` directory and a separate output directory. Its validation and scientific outputs are unchanged. `audit_results.py` accepts a result directory, the corresponding pinned data directory, and an audit JSON path.

`create` requires every expected artifact, rejects empty files and symlinked artifacts, verifies the three upstream file fingerprints, and compares the five sample-level tables against the source annotation. It creates a new manifest exclusively; an existing manifest is never overwritten. `verify` repeats those checks and compares the entire snapshot. Both commands return exit code 0 on success and 2 on an expected validation or file error. The audit JSON contains filenames, byte counts, SHA-256 values, ordered sample names, analysis version, source commit, and a statement of scope.

The supplied [example snapshot](../audits/example-artifacts.json) pins the committed example's bytes. It is not a cross-platform numerical tolerance test. A newly generated analysis should receive its own snapshot; PNG bytes can differ between environments even when the scientific values agree.

`compare_results.py REFERENCE CANDIDATE` supplies the separate numerical check. It expects the 14 named CSV tables and `summary.json` plus `design_diagnostics.json`. Column schemas, ordered categories/sample IDs, integer counts, JSON keys/types, and recorded source hashes must match exactly. Floating values must be finite and satisfy `abs(a-b) <= max(atol, rtol * max(abs(a), abs(b)))`, using the [Python `math.isclose` definition](https://docs.python.org/3.12/library/math.html#math.isclose). Defaults are `rtol=1e-7` and `atol=1e-8`. The absolute floor is larger than one unit of the CSVs' nine-decimal printing precision; the relative allowance accommodates small floating-point variation in larger values. These are engineering tolerances, not biological or statistical uncertainty limits.

The comparator rejects missing files, changed row counts, malformed tables, duplicate JSON keys, and nonfinite values. Matrix row/column sample identities are checked; matrix symmetry is not independently tested here. It reports every difference count and details the first 100 with file/row/field locations. Version-string differences for the same package keys in `summary.json` are recorded separately without failing numeric agreement; malformed environment metadata still fails. No other metadata difference is ignored. Exit codes are 0 for agreement, 1 for differences or malformed/missing result inputs, and 2 for invalid options or report-writing errors. `--report PATH` saves the machine-readable comparison and cannot overwrite a compared input.

## Worked use case

Check the supplied package before copying it to another location:

```sh
python audit_results.py verify --results results --data data --manifest audits/example-artifacts.json
```

Expected output: `verify: 23 files and 7 sample labels verified`.

For a fresh local run, keep the original results and create a separate audit:

```sh
python analysis.py --out my-results
python compare_results.py results my-results --report my-comparison.json
python audit_results.py create --results my-results --data data --manifest my-audit.json
python audit_results.py verify --results my-results --data data --manifest my-audit.json
```

Run the last command again after transferring the folder and manifest. A changed figure, missing CSV, reordered sample table, or changed sample label causes failure. Files outside the declared inventory are ignored; the audit does not claim a clean or exclusive output directory.

The [example numerical comparison](../examples/numerical-comparison.json) records a fresh local reproduction of all 16 numeric files. A comparison does not inspect PNG pixels or prove the reference analysis is correct; the separate artifact snapshot still covers the files chosen for transfer.

## Checks before relying on a result

The software checks file integrity and selected table relationships. The operator must separately inspect the design, provenance, and analysis assumptions. In particular, keep the package labels as supplied and read the [source discrepancy](../data/README.md) before interpreting condition groups. A matching checksum proves equality to the chosen reference bytes, not that those labels are biologically correct.

Use a new output directory for a new run, retain its parameters and environment record, and avoid editing files concurrently with an audit. Independently compare numeric results when reproducing across environments. The existing tests cover normalization and PCA calculations; this audit adds neither a significance test nor a universal sample-acceptance threshold.

## Known limits

The audit is a snapshot, not a signed chain of custody: someone could deliberately replace both data and the manifest. Creating a snapshot does not prove that all files were generated together; already mixed artifacts with matching labels can pass. The operator must first choose a reviewed run. It validates sample metadata in five tables but does not verify every numerical relationship in every output. The numerical comparator uses the selected reference as its authority; two identically wrong results can agree, and a scientifically meaningful change smaller than its declared tolerance can pass. Both tools read files sequentially, so concurrent changes are outside their contracts. The analysis still has only seven experimental libraries, unresolved source-label allocation, and no absolute-RNA reference. A full-rank design does not remove unmeasured confounding.

## Next improvements

1. **More matrix relationships:** independently check correlation/distance symmetry and diagonals; acceptance means asymmetric matrices fail even when they match a deliberately corrupted reference.
2. **Source-label reconciliation:** locate and archive definitive sample-accession evidence; acceptance requires a documented mapping and review before any relabeling, not an inferred correction.
3. **Environment recreation:** test a second supported Python platform; acceptance means the full numeric comparator passes and environment differences are recorded without demanding identical PNG bytes.
4. **Tolerance review:** test documented near-boundary changes across the supported environments; acceptance means the adopted tolerance has explicit examples of allowed rounding and rejected scientific changes.
