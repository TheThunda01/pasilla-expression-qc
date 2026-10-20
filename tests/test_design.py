import unittest

import pandas as pd

from design import inspect_design


def metadata(conditions, read_types):
    return pd.DataFrame({"condition": conditions, "type": read_types},
                        index=[f"sample{i}" for i in range(len(conditions))])


class DesignChecks(unittest.TestCase):
    def test_crossed_design_is_full_rank_with_explicit_baselines(self):
        data = metadata(["untreated", "untreated", "treated", "treated"],
                        ["single-read", "paired-end", "single-read", "paired-end"])
        cross, matrix, result = inspect_design(data)
        self.assertEqual(cross.to_numpy().tolist(), [[1, 1], [1, 1]])
        self.assertEqual(matrix.to_numpy().tolist(), [[1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]])
        self.assertEqual(result["rank"], 3)
        self.assertFalse(result["perfect_confounding"])
        self.assertEqual(result["residual_dimension"], 1)
        self.assertFalse(result["model_fitted"])

    def test_perfect_confounding_is_detected(self):
        data = metadata(["untreated", "untreated", "treated", "treated"],
                        ["single-read", "single-read", "paired-end", "paired-end"])
        cross, _, result = inspect_design(data)
        self.assertEqual(result["rank"], 2)
        self.assertTrue(result["perfect_confounding"])
        self.assertFalse(result["full_column_rank"])
        self.assertEqual(int(cross.loc["treated", "single-read"]), 0)

    def test_reverse_perfect_confounding_is_detected(self):
        data = metadata(["untreated", "untreated", "treated", "treated"],
                        ["paired-end", "paired-end", "single-read", "single-read"])
        self.assertTrue(inspect_design(data)[2]["perfect_confounding"])

    def test_missing_level_is_not_misreported_as_no_confounding(self):
        data = metadata(["untreated", "untreated"], ["single-read", "paired-end"])
        _, matrix, result = inspect_design(data)
        self.assertEqual(matrix.shape[1], 3)
        self.assertEqual(result["missing_levels"]["condition"], ["treated"])
        self.assertIsNone(result["perfect_confounding"])
        self.assertFalse(result["full_column_rank"])

    def test_sparse_nonperfect_design_can_still_have_full_rank(self):
        data = metadata(["untreated", "untreated", "treated"],
                        ["single-read", "paired-end", "paired-end"])
        _, _, result = inspect_design(data)
        self.assertEqual(result["rank"], 3)
        self.assertFalse(result["perfect_confounding"])
        self.assertFalse(result["all_design_cells_observed"])
        self.assertEqual(result["residual_dimension"], 0)

    def test_unknown_missing_and_duplicate_annotations_fail(self):
        fixtures = [metadata(["control"], ["single-read"]),
                    metadata([None], ["single-read"]), pd.DataFrame(),
                    pd.DataFrame({"condition": ["treated", "untreated"],
                                  "type": ["paired-end", "single-read"]}, index=["a", "a"])]
        for fixture in fixtures:
            with self.subTest(fixture=fixture), self.assertRaises(ValueError):
                inspect_design(fixture)


if __name__ == "__main__":
    unittest.main()
