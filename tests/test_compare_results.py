from contextlib import redirect_stderr
import csv
import io
import json
from pathlib import Path
import shutil
import tempfile
import unittest
import compare_results as compare

ROOT = Path(__file__).resolve().parents[1]


class NumericComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.target = Path(self.temp.name)/'candidate'
        shutil.copytree(ROOT/'results', self.target)

    def edit(self, filename, change):
        path = self.target/filename
        with path.open(newline='') as stream:
            reader = csv.DictReader(stream); fields, rows = reader.fieldnames, list(reader)
        change(rows)
        with path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fields); writer.writeheader(); writer.writerows(rows)

    def run_comparison(self):
        return compare.compare(ROOT/'results', self.target)

    def test_exact_reference_and_harmless_roundoff_pass(self):
        result = self.run_comparison()
        self.assertEqual(result['status'], 'match')
        self.assertEqual(result['files_compared'], 16)
        self.edit('normalization_factors.csv', lambda rows: rows[0].update(
            total_count_factor=str(float(rows[0]['total_count_factor'])+1e-9)))
        self.assertEqual(self.run_comparison()['status'], 'match')

    def test_meaningful_numeric_change_is_located(self):
        self.edit('normalization_factors.csv', lambda rows: rows[0].update(total_count_factor='99'))
        result = self.run_comparison()
        self.assertEqual(result['status'], 'different')
        self.assertTrue(any(d['location']=='normalization_factors.csv:2:total_count_factor' for d in result['differences']))

    def test_identities_and_integer_counts_are_exact(self):
        for field, value in [('sample', '01'), ('count_table_total', '13972512.0')]:
            with self.subTest(field=field):
                shutil.copyfile(ROOT/'results/sample_qc.csv', self.target/'sample_qc.csv')
                self.edit('sample_qc.csv', lambda rows: rows[0].update({field:value}))
                self.assertEqual(self.run_comparison()['status'], 'different')

    def test_matrix_row_order_cannot_be_hidden_by_numeric_tolerance(self):
        self.edit('distances_total_count.csv', lambda rows: rows.reverse())
        result = self.run_comparison()
        self.assertTrue(any('matrix row identities' in d['reason'] for d in result['differences']))

    def test_nonfinite_csv_values_fail(self):
        for value in ['NaN', 'inf', '-inf', '1e999']:
            with self.subTest(value=value):
                self.edit('normalization_factors.csv', lambda rows: rows[0].update(total_count_factor=value))
                self.assertEqual(self.run_comparison()['status'], 'different')

    def test_missing_file_and_schema_change_are_reported(self):
        (self.target/'pca_logcounts.csv').unlink()
        self.edit('pca_logcpm.csv', lambda rows: rows.pop())
        result = self.run_comparison()
        self.assertEqual(result['status'], 'different')
        self.assertTrue(any('row count' in d['reason'] for d in result['differences']))
        self.assertTrue(any('read failure' in d['reason'] for d in result['differences']))

    def test_json_duplicate_keys_nonfinite_and_boolean_type_fail(self):
        for text in ['{"rank":1,"rank":1}', '{"rank":NaN}', '{"rank":1e999}']:
            with self.subTest(text=text):
                (self.target/'design_diagnostics.json').write_text(text)
                self.assertEqual(self.run_comparison()['status'], 'different')
        original = json.loads((ROOT/'results/summary.json').read_text())
        original['samples'] = True
        (self.target/'summary.json').write_text(json.dumps(original))
        self.assertEqual(self.run_comparison()['status'], 'different')

    def test_environment_difference_is_recorded_without_a_scientific_mismatch(self):
        path = self.target/'summary.json'; data=json.loads(path.read_text())
        data['package_versions']['numpy']='different-environment'
        path.write_text(json.dumps(data))
        result=self.run_comparison()
        self.assertEqual(result['status'],'match')
        self.assertEqual(len(result['environment_differences']),1)

    def test_environment_schema_changes_are_not_ignored(self):
        path=self.target/'summary.json'; data=json.loads(path.read_text())
        data['package_versions']=[];path.write_text(json.dumps(data))
        self.assertEqual(self.run_comparison()['status'],'different')

    def test_matching_malformed_json_roots_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            reference=Path(temp)/'reference';shutil.copytree(ROOT/'results',reference)
            for directory in (reference,self.target):
                (directory/'design_diagnostics.json').write_text('{"rank":3}')
            self.assertEqual(compare.compare(reference,self.target)['status'],'different')

    def test_bad_tolerance_and_input_overwrite_return_two(self):
        with redirect_stderr(io.StringIO()):
            self.assertEqual(compare.main([str(ROOT/'results'),str(self.target),'--rtol','nan']),2)
            path=self.target/'summary.json'; before=path.read_bytes()
            self.assertEqual(compare.main([str(ROOT/'results'),str(self.target),'--report',str(path)]),2)
            self.assertEqual(before,path.read_bytes())


if __name__ == '__main__':
    unittest.main()
