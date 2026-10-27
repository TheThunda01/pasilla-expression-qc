"""Reject misaligned sample labels and incomplete normalization comparisons."""
import csv
import importlib.util
from pathlib import Path
import shutil
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('sample_explorer', ROOT / 'scripts/build_sample_explorer.py')
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class SampleExplorerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)
        for name in ('pca_total_count.csv', 'pca_median_ratio.csv', 'normalization_factors.csv', 'leave_one_out_stability.csv'):
            shutil.copyfile(ROOT / 'results' / name, self.path / name)

    def edit(self, name, change):
        path = self.path / name
        with path.open(newline='') as stream:
            rows = list(csv.DictReader(stream))
        fields = list(rows[0])
        change(rows)
        with path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    def test_reference_keeps_all_samples_and_omissions(self):
        data, hashes = MODULE.load_results(self.path)
        self.assertEqual((len(data['samples']), len(data['omissions']), len(hashes)), (7, 14, 4))

    def test_missing_coordinate_is_rejected(self):
        self.edit('pca_total_count.csv', lambda rows: rows.pop())
        with self.assertRaisesRegex(ValueError, 'every sample'):
            MODULE.load_results(self.path)

    def test_conflicting_read_type_is_rejected(self):
        self.edit('pca_median_ratio.csv', lambda rows: rows[0].update(type='paired-end'))
        with self.assertRaisesRegex(ValueError, 'metadata disagree'):
            MODULE.load_results(self.path)

    def test_missing_omission_is_rejected(self):
        self.edit('leave_one_out_stability.csv', lambda rows: rows.pop())
        with self.assertRaisesRegex(ValueError, 'both omission'):
            MODULE.load_results(self.path)

    def test_wrong_factor_ratio_is_rejected(self):
        self.edit('normalization_factors.csv', lambda rows: rows[0].update(median_to_total_factor_ratio='1.9'))
        with self.assertRaisesRegex(ValueError, 'ratio disagrees'):
            MODULE.load_results(self.path)

    def test_nonfinite_coordinate_is_rejected(self):
        self.edit('pca_total_count.csv', lambda rows: rows[0].update(PC1='inf'))
        with self.assertRaisesRegex(ValueError, 'finite'):
            MODULE.load_results(self.path)

    def test_duplicate_sample_is_rejected(self):
        self.edit('normalization_factors.csv', lambda rows: rows.append(rows[0]))
        with self.assertRaisesRegex(ValueError, 'duplicate sample'):
            MODULE.load_results(self.path)
