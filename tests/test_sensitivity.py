"""Check mathematical behavior with small independent, hand-checkable examples."""
import unittest

import numpy as np

from sensitivity import (geometric_mean, leave_one_out, pairwise_summaries,
                         principal_components, procrustes_residual, size_factors,
                         transformed_counts)


class NormalizationTests(unittest.TestCase):
    def test_known_library_multiplier_is_removed(self):
        matrix = np.array([[10, 20], [30, 60], [50, 100]], dtype=float)
        factors, eligible = size_factors(matrix)
        self.assertEqual(eligible, 3)
        for factor in factors.values():
            np.testing.assert_allclose(factor, [1 / np.sqrt(2), np.sqrt(2)])
            transformed = transformed_counts(matrix, factor, geometric_mean(matrix.sum(axis=0)))
            np.testing.assert_allclose(transformed[:, 0], transformed[:, 1])

    def test_one_high_count_gene_affects_total_more_than_median(self):
        matrix = np.array([[100, 100], [200, 200], [300, 300], [10, 10000]], dtype=float)
        factors, _ = size_factors(matrix)
        np.testing.assert_allclose(factors['median_ratio'], [1, 1])
        self.assertGreater(factors['total_count'][1] / factors['total_count'][0], 17)

    def test_zero_count_genes_do_not_enter_reference(self):
        matrix = np.array([[10, 20], [30, 60], [0, 500]], dtype=float)
        factors, eligible = size_factors(matrix)
        self.assertEqual(eligible, 2)
        np.testing.assert_allclose(factors['median_ratio'], [1 / np.sqrt(2), np.sqrt(2)])

    def test_no_positive_reference_is_explicit_failure(self):
        with self.assertRaisesRegex(ValueError, 'positive in every sample'):
            size_factors(np.array([[0, 3], [4, 0]]))

    def test_invalid_matrices_fail(self):
        for values in [np.array([[1, -1], [2, 3]]), np.array([[1, np.nan], [2, 3]]),
                       np.array([[0, 3], [0, 2]]), np.ones((1, 4))]:
            with self.subTest(values=values), self.assertRaises(ValueError):
                size_factors(values)

    def test_reference_depth_controls_only_declared_scale(self):
        matrix = np.array([[10, 20, 40], [30, 60, 120]], dtype=float)
        factors, _ = size_factors(matrix)
        # Total-count normalization at a fixed reference is invariant to which
        # other samples are present, before any PCA refit.
        full = transformed_counts(matrix, factors['total_count'], 100)
        subset_factors, _ = size_factors(matrix[:, :2])
        subset = transformed_counts(matrix[:, :2], subset_factors['total_count'], 100)
        np.testing.assert_allclose(full[:, :2], subset)

    def test_factor_scale_is_centered(self):
        matrix = np.array([[20, 45, 80], [50, 100, 180], [70, 30, 15]])
        factors, _ = size_factors(matrix)
        for factor in factors.values():
            self.assertAlmostEqual(geometric_mean(factor), 1)


class GeometryTests(unittest.TestCase):
    def test_pairwise_rms_matches_hand_calculation(self):
        distance, correlation = pairwise_summaries(np.array([[1, 3], [2, 5], [4, 2]]))
        self.assertAlmostEqual(distance[0, 1], np.sqrt(2.5))
        np.testing.assert_allclose(distance, distance.T)
        np.testing.assert_allclose(np.diag(distance), 0)
        np.testing.assert_allclose(np.diag(correlation), 1)
        self.assertAlmostEqual(correlation[0, 1], 1)

    def test_constant_profile_correlation_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'constant'):
            pairwise_summaries(np.array([[2, 2], [1, 3]]))

    def test_rotation_reflection_translation_and_scale_do_not_create_instability(self):
        reference = np.array([[0., 0.], [1., 0.], [0., 2.], [3., 2.]])
        rotation_reflection = np.array([[0., 1.], [1., 0.]])
        candidate = 7 * reference @ rotation_reflection + [50, -30]
        self.assertLess(procrustes_residual(reference, candidate), 1e-12)

    def test_shape_change_is_detected(self):
        reference = np.array([[0., 0.], [1., 0.], [0., 2.], [3., 2.]])
        candidate = reference.copy()
        candidate[0, 1] = 3
        self.assertGreater(procrustes_residual(reference, candidate), .1)

    def test_procrustes_rejects_invalid_inputs(self):
        for candidate in [np.ones((2, 2)), np.ones((4, 2)), np.full((4, 2), np.nan)]:
            with self.assertRaises(ValueError):
                procrustes_residual(np.arange(8).reshape(4, 2), candidate)

    def test_leave_one_out_refits_every_sample_and_keeps_gene_set(self):
        counts = np.array([[10, 22, 13, 40], [30, 12, 33, 28], [50, 90, 51, 60], [5, 8, 4, 6.]])
        keep = np.array([True, True, True, False])
        factors, _ = size_factors(counts)
        reference = geometric_mean(counts.sum(axis=0))
        full = {method: principal_components(transformed_counts(counts, sf, reference)[keep].T)[0]
                for method, sf in factors.items()}
        result = leave_one_out(counts, keep, ['a', 'b', 'c', 'd'], full, reference)
        self.assertEqual(len(result), 8)
        self.assertEqual(result.omitted_sample.value_counts().to_dict(), {'a': 2, 'b': 2, 'c': 2, 'd': 2})
        self.assertTrue((result.fixed_pca_genes == 3).all())
        self.assertTrue((result.remaining_samples == 3).all())
        self.assertTrue(np.isfinite(result.procrustes_residual).all())


if __name__ == '__main__':
    unittest.main()
