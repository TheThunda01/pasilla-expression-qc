import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

import analysis


class FullPublicDataPipeline(unittest.TestCase):
    def test_outputs_are_complete_and_numerically_consistent(self):
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary)
            summary = analysis.run(analysis.ROOT / 'data', out)
            self.assertEqual(summary['genes_after_filter'], 8148)
            self.assertEqual(summary['samples'], 7)
            factors = pd.read_csv(out / 'normalization_factors.csv')
            for field in ['total_count_factor', 'median_ratio_factor']:
                self.assertAlmostEqual(float(np.exp(np.log(factors[field]).mean())), 1, places=8)
            stability = pd.read_csv(out / 'leave_one_out_stability.csv')
            self.assertEqual(len(stability), 14)
            self.assertTrue((stability.fixed_pca_genes == 8148).all())
            for method in ['total_count', 'median_ratio']:
                matrix = pd.read_csv(out / f'distances_{method}.csv', index_col=0).to_numpy()
                np.testing.assert_allclose(matrix, matrix.T, atol=1e-8)
                np.testing.assert_allclose(np.diag(matrix), 0, atol=1e-8)
                self.assertTrue(np.isfinite(matrix).all())
            disk_summary = json.loads((out / 'summary.json').read_text())
            self.assertEqual(disk_summary, summary)
            self.assertGreater((out / 'normalization_comparison.png').stat().st_size, 1000)
            self.assertGreater((out / 'leave_one_out_stability.png').stat().st_size, 1000)


if __name__ == '__main__':
    unittest.main()
