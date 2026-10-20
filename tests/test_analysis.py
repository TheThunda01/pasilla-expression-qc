import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

import analysis


class InputChecks(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def counts(self, content):
        path = self.root / "counts.tsv"
        path.write_text(content)
        return analysis.load_counts(path)

    def test_valid_integer_counts_and_order(self):
        frame = self.counts("gene_id\tb\ta\ng1\t10\t0\ng2\t0\t5\n")
        self.assertEqual(list(frame.columns), ["b", "a"])
        self.assertEqual(frame.loc["g1", "b"], 10)

    def test_bad_counts_are_rejected(self):
        for value in ["-1", "1.2", "", "NaN", "inf", "1e2"]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.counts(f"gene_id\ta\tb\ng1\t{value}\t2\n")

    def test_duplicate_genes_and_samples_are_rejected(self):
        for content in ["gene_id\ta\ta\ng1\t1\t2\n", "gene_id\ta\tb\ng1\t1\t2\ng1\t3\t4\n"]:
            with self.assertRaises(ValueError):
                self.counts(content)

    def test_zero_library_is_rejected(self):
        with self.assertRaises(ValueError):
            self.counts("gene_id\ta\tb\ng1\t0\t2\ng2\t0\t1\n")

    def test_ragged_row_is_rejected(self):
        with self.assertRaises(ValueError):
            self.counts("gene_id\ta\tb\ng1\t1\n")

    def test_metadata_is_aligned_and_not_positionally_joined(self):
        path = self.root / "meta.csv"
        path.write_text("file,condition,type\nbfb,treated,paired-end\nafb,untreated,single-read\n")
        meta = analysis.load_metadata(path, ["a", "b"])
        self.assertEqual(list(meta.index), ["a", "b"])
        self.assertEqual(meta.loc["a", "condition"], "untreated")

    def test_metadata_extra_missing_and_duplicates_fail(self):
        path = self.root / "meta.csv"
        for text in ["afb,treated,single-read\n", "afb,treated,single-read\nafb,treated,single-read\n",
                     "a,treated,single-read\nbfb,treated,single-read\n"]:
            path.write_text("file,condition,type\n" + text)
            with self.assertRaises(ValueError):
                analysis.load_metadata(path, ["a", "b"])

    def test_pinned_data_integrity_and_shape(self):
        analysis.verify_manifest(analysis.ROOT / "data")
        frame = analysis.load_counts(analysis.ROOT / "data/pasilla_gene_counts.tsv")
        self.assertEqual(frame.shape, (14599, 7))

    def test_modified_data_is_rejected(self):
        (self.root / "x").write_text("changed")
        (self.root / "manifest.json").write_text(json.dumps({"files": [{"path": "x", "bytes": 0, "sha256": "abc"}]}))
        with self.assertRaises(ValueError):
            analysis.verify_manifest(self.root)


class PCAChecks(unittest.TestCase):
    def test_centering_and_variance(self):
        values = np.array([[0, 0], [1, 1], [2, 2]])
        scores, ratios = analysis.principal_components(values)
        np.testing.assert_allclose(scores.mean(axis=0), 0, atol=1e-12)
        self.assertAlmostEqual(ratios[0], 1)
        self.assertAlmostEqual(ratios.sum(), 1)

    def test_no_gene_standardization(self):
        values = np.array([[0, 0], [10, 0], [0, 1], [10, 1]])
        _, ratios = analysis.principal_components(values)
        self.assertGreater(ratios[0], 0.98)

    def test_invalid_pca_matrices_fail(self):
        for values in [np.ones((3, 2)), np.array([[1, np.nan], [2, 1]]), np.ones((1, 3))]:
            with self.assertRaises(ValueError):
                analysis.principal_components(values)


if __name__ == "__main__":
    unittest.main()
