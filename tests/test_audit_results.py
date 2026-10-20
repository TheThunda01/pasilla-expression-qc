import json
from pathlib import Path
import shutil
import tempfile
import unittest
import audit_results as audit

ROOT = Path(__file__).resolve().parents[1]


class ArtifactAuditTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        shutil.copytree(ROOT/'results', self.root/'results')
        shutil.copytree(ROOT/'data', self.root/'data')
        self.manifest = self.root/'audit.json'

    def test_complete_snapshot_can_be_verified(self):
        result = audit.create(self.root/'results', self.root/'data', self.manifest)
        self.assertEqual(len(result['files']), 23)
        self.assertEqual(len(result['samples_in_order']), 7)
        self.assertEqual(result, audit.verify(self.root/'results', self.root/'data', self.manifest))
        with self.assertRaises(FileExistsError):
            audit.create(self.root/'results', self.root/'data', self.manifest)

    def test_changed_result_is_detected(self):
        audit.create(self.root/'results', self.root/'data', self.manifest)
        with (self.root/'results/report.md').open('a') as stream:
            stream.write('\nChanged\n')
        with self.assertRaisesRegex(ValueError, 'differs'):
            audit.verify(self.root/'results', self.root/'data', self.manifest)

    def test_mixed_sample_order_is_rejected_even_before_snapshot(self):
        path = self.root/'results/pca_logcpm.csv'
        lines = path.read_text().splitlines()
        lines[1], lines[2] = lines[2], lines[1]
        path.write_text('\n'.join(lines)+'\n')
        with self.assertRaisesRegex(ValueError, 'Sample order'):
            audit.snapshot(self.root/'results', self.root/'data')

    def test_changed_label_is_rejected(self):
        path = self.root/'results/sample_qc.csv'
        path.write_text(path.read_text().replace('treated1fb,treated,', 'treated1fb,untreated,', 1))
        with self.assertRaisesRegex(ValueError, 'Annotation'):
            audit.snapshot(self.root/'results', self.root/'data')

    def test_missing_output_or_changed_source_is_rejected(self):
        (self.root/'results/count_totals.png').unlink()
        with self.assertRaises(ValueError):
            audit.snapshot(self.root/'results', self.root/'data')
        shutil.copyfile(ROOT/'results/count_totals.png', self.root/'results/count_totals.png')
        path = self.root/'data/pasilla_sample_annotation.csv'
        path.write_bytes(path.read_bytes()+b'\n')
        with self.assertRaisesRegex(ValueError, 'Input fingerprint'):
            audit.snapshot(self.root/'results', self.root/'data')

    def test_snapshot_cannot_overwrite_a_result(self):
        path = self.root/'results/report.md'
        original = path.read_bytes()
        with self.assertRaisesRegex(ValueError, 'overwrite'):
            audit.create(self.root/'results', self.root/'data', path)
        self.assertEqual(original, path.read_bytes())


if __name__ == '__main__':
    unittest.main()
