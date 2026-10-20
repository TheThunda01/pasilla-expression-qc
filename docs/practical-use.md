# Practical use

## Where this fits

This is a small, inspectable check of a published count table before formal expression modeling. It is useful for learning how sample annotations, count totals, normalization, and PCA sensitivity fit together. It does not settle the original article/package label discrepancy or test differential expression.

The new artifact audit addresses a separate operational problem: a folder can accidentally contain tables and figures from different runs. It checks the expected 23 files, their exact bytes, the input fingerprints recorded in the summary, and the sample ordering and annotations in five sample-level tables. It does not recompute the analysis.

## Input and output contract

`analysis.py` continues to accept the pinned `data/` directory and a separate output directory. Its validation and scientific outputs are unchanged. `audit_results.py` accepts a result directory, the corresponding pinned data directory, and an audit JSON path.

`create` requires every expected artifact, rejects empty files and symlinked artifacts, verifies the three upstream file fingerprints, and compares the five sample-level tables against the source annotation. It creates a new manifest exclusively; an existing manifest is never overwritten. `verify` repeats those checks and compares the entire snapshot. Both commands return exit code 0 on success and 2 on an expected validation or file error. The audit JSON contains filenames, byte counts, SHA-256 values, ordered sample names, analysis version, source commit, and a statement of scope.

The supplied [example snapshot](../audits/example-artifacts.json) pins the committed example's bytes. It is not a cross-platform numerical tolerance test. A newly generated analysis should receive its own snapshot; PNG bytes can differ between environments even when the scientific values agree.

## Worked use case

Check the supplied package before copying it to another location:

```sh
python audit_results.py verify --results results --data data --manifest audits/example-artifacts.json
```

Expected output: `verify: 23 files and 7 sample labels verified`.

For a fresh local run, keep the original results and create a separate audit:

```sh
python analysis.py --out my-results
python audit_results.py create --results my-results --data data --manifest my-audit.json
python audit_results.py verify --results my-results --data data --manifest my-audit.json
```

Run the last command again after transferring the folder and manifest. A changed figure, missing CSV, reordered sample table, or changed sample label causes failure. Files outside the declared inventory are ignored; the audit does not claim a clean or exclusive output directory.

## Checks before relying on a result

The software checks file integrity and selected table relationships. The operator must separately inspect the design, provenance, and analysis assumptions. In particular, keep the package labels as supplied and read the [source discrepancy](../data/README.md) before interpreting condition groups. A matching checksum proves equality to the chosen reference bytes, not that those labels are biologically correct.

Use a new output directory for a new run, retain its parameters and environment record, and avoid editing files concurrently with an audit. Independently compare numeric results when reproducing across environments. The existing tests cover normalization and PCA calculations; this audit adds neither a significance test nor a universal sample-acceptance threshold.

## Known limits

The audit is a snapshot, not a signed chain of custody: someone could deliberately replace both data and the manifest. Creating a snapshot does not prove that all files were generated together; already mixed artifacts with matching labels can pass. The operator must first choose a reviewed run. It validates sample metadata in five tables but does not verify every numerical relationship in every output. It reads files sequentially, so concurrent changes are outside its contract. The analysis still has seven related experimental samples, unresolved source-label allocation, and no absolute-RNA reference. A full-rank design does not remove unmeasured confounding.

## Next improvements

1. **Numerical result comparison:** add a documented tolerance-based comparator; acceptance means expected rounding passes while a changed normalization factor fails.
2. **More cross-file relationships:** check correlation/distance row and column identities and matrix symmetry; acceptance means swapped columns and asymmetric matrices fail independently of checksums.
3. **Source-label reconciliation:** locate and archive definitive sample-accession evidence; acceptance requires a documented mapping and review before any relabeling, not an inferred correction.
4. **Environment recreation:** test a second supported Python platform; acceptance means the full numeric comparator passes and environment differences are recorded without demanding identical PNG bytes.
