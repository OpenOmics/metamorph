import unittest
import numpy as np
from expand_taxonomy import (canonicalize, parse_path, build_tree, reconcile,
                             exclusive_species, RESERVED)


def calculate(paths, values, atol=0.001):
    nodes, depths, parents, source = build_tree([canonicalize(p) for p in paths])
    residual, totals, issues = reconcile(np.array(values, dtype=float), nodes,
                                         depths, parents, source, atol, 1e-8)
    keys, species = exclusive_species(residual, nodes)
    return nodes, residual, totals[parents < 0].sum(axis=0), issues, keys, species


class HierarchyTests(unittest.TestCase):
    def test_no_parent_double_count(self):
        _, _, total, issues, _, species = calculate(
            ['k__Bacteria', 'k__Bacteria|p__A', 'k__Bacteria|p__B'],
            [[100, 100], [60, 20], [40, 80]])
        np.testing.assert_allclose(total, [100, 100])
        np.testing.assert_allclose(species.sum(axis=0), total)
        self.assertEqual(issues, [])

    def test_missing_ancestors(self):
        nodes, _, total, issues, keys, species = calculate(
            ['k__Bacteria|p__A|c__C', 'k__Bacteria|p__B'], [[25], [75]])
        self.assertIn('k__Bacteria', nodes)
        np.testing.assert_allclose(total, [100])
        np.testing.assert_allclose(species.sum(axis=0), [100])
        self.assertTrue(all(len(k) == 7 for k in keys))
        self.assertFalse(issues)

    def test_parent_residual_preserved(self):
        _, _, total, issues, keys, species = calculate(
            ['k__Bacteria', 'k__Bacteria|p__A'], [[100], [70]])
        self.assertEqual(issues[0][2], 'unresolved_parent_mass')
        np.testing.assert_allclose(species.sum(axis=0), total)
        self.assertTrue(any(RESERVED in '|'.join(k) for k in keys))

    def test_overrun_fails(self):
        _, _, _, issues, _, _ = calculate(
            ['k__Bacteria', 'k__Bacteria|p__A'], [[50], [60]])
        self.assertTrue(any(i[2].startswith('ERROR') for i in issues))

    def test_rounding_logged_and_reconciled(self):
        _, _, total, issues, _, species = calculate(
            ['k__Bacteria', 'k__Bacteria|p__A'], [[99.99999], [100]])
        self.assertEqual(issues[0][2], 'rounding_parent_raised')
        np.testing.assert_allclose(total, [100])
        np.testing.assert_allclose(species.sum(axis=0), [100])

    def test_multiple_kingdoms_and_global_unknown(self):
        _, _, total, _, _, species = calculate(
            ['k__Bacteria', 'k__Archaea', 'UNCLASSIFIED'], [[85], [5], [10]])
        np.testing.assert_allclose(total, [100])
        np.testing.assert_allclose(species.sum(axis=0), total)

    def test_strain_and_species_not_double_counted(self):
        prefix = 'k__Bacteria|p__P|c__C|o__O|f__F|g__G|s__S'
        _, _, total, _, keys, species = calculate(
            [prefix, prefix + '|t__SGB1', prefix + '|t__SGB2'], [[100], [30], [50]])
        self.assertEqual(len(keys), 1)
        np.testing.assert_allclose(total, [100])
        np.testing.assert_allclose(species, [[100]])

    def test_missing_ranks_and_taxonkit_fill(self):
        p = canonicalize('k__Bacteria|g__Example',
                         ['Bacteria', 'P', 'C', 'O', 'F', 'Example', 'Unassigned'])
        self.assertEqual(p, 'k__Bacteria|p__P|c__C|o__O|f__F|g__Example')
        self.assertIn('p__' + RESERVED, canonicalize('k__Bacteria|p__|g__Example'))
        self.assertEqual(canonicalize('k__Bacteria|p__Old', ['Bacteria', 'New']),
                         'k__Bacteria|p__Old')

    def test_malformed_paths_rejected(self):
        for path in ['Bacteria', 'g__G|k__Bacteria', 'k__Bacteria|k__Other', 'd__Bacteria']:
            with self.assertRaises(ValueError):
                parse_path(path)

    def test_normalization_is_shared_across_ranks(self):
        _, residual, total, _, _, species = calculate(
            ['k__Bacteria', 'k__Bacteria|p__A'], [[80, 40], [60, 30]])
        normalized = species * (100 / total)
        np.testing.assert_allclose(normalized.sum(axis=0), [100, 100])
        # Every taxonomy prefix is a partition of these exclusive contributions.
        np.testing.assert_allclose(residual.sum(axis=0), total)


# This integration test runs automatically when R + phyloseq + data.table exist.
# It is skipped in environments without those dependencies.
import csv
import json
import shutil
import subprocess
import tempfile
from pathlib import Path


def has_r_packages():
    if not shutil.which('Rscript'):
        return False
    return subprocess.run(['Rscript', str(Path(__file__).with_name('aggregate_phyloseq.R')), '--check'],
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0


@unittest.skipUnless(has_r_packages(), 'R, phyloseq, and data.table required')
class PhyloseqIntegrationTests(unittest.TestCase):
    def test_full_pipeline_batches_and_totals(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            source, out = d / 'input.csv', d / 'out'
            source.write_text('clade_name,A,B\nk__Bacteria,80,100\n'
                              'k__Bacteria|p__P,50,60\n'
                              'k__Bacteria|p__Q|c__C,20,40\n')
            subprocess.run(['python3', str(Path(__file__).with_name('expand_taxonomy.py')),
                            str(source), '--out', str(out), '--batch-size', '1'], check=True)
            self.assertEqual(json.loads((out / 'run.json').read_text())['status'], 'complete')
            for batch, sample in [('batch_0001', 'A'), ('batch_0002', 'B')]:
                for rank in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
                    with (out / batch / (rank + '.csv')).open() as f:
                        rows = list(csv.DictReader(f))
                    self.assertAlmostEqual(sum(float(r[sample]) for r in rows), 100, places=7)
                with (out / batch / 'all_ranks.csv').open() as f:
                    all_rows = list(csv.DictReader(f))
                self.assertTrue(any(r['clade_name'] == 'k__Bacteria' for r in all_rows))

if __name__ == '__main__':
    unittest.main()
