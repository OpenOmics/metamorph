#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""Unit tests for src.utils.valid_input(), metamorph's sample sheet parser/validator."""

import os
import sys
import csv
import tempfile
import unittest
from argparse import ArgumentTypeError

# Make `src` importable when this file is run directly or via
# `python3 -m unittest discover` from the repository root.
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from src.utils import valid_input


class ValidInputTestCase(unittest.TestCase):
    """Base class that sets up a scratch directory with placeholder fastq
    files and helpers for writing sample sheets that reference them."""

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        # Placeholder fastq files: contents are irrelevant, valid_input()
        # only checks that each path in the sheet exists.
        self.files = {}
        for name in (
            'dna_s1_R1', 'dna_s1_R2', 'dna_s2_R1', 'dna_s2_R2',
            'rna_s1_R1', 'rna_s1_R2', 'rna_s2_R1', 'rna_s2_R2',
        ):
            path = os.path.join(self.tmpdir.name, f'{name}.fastq.gz')
            open(path, 'w').close()
            self.files[name] = path

    def _write_sheet(self, rows, header=('DNA',), ext='.tsv'):
        """Writes a delimited sample sheet and returns its path.
        `rows` is a list of tuples, one value per column in `header`.
        """
        delim = ',' if ext == '.csv' else '\t'
        path = os.path.join(self.tmpdir.name, f'samplesheet{ext}')
        with open(path, 'w', newline='') as fh:
            writer = csv.writer(fh, delimiter=delim)
            writer.writerow(header)
            for row in rows:
                writer.writerow(row)
        return path


class TestValidSheets(ValidInputTestCase):
    """Sheets that should parse successfully."""

    def test_dna_only_sheet(self):
        sheet = self._write_sheet([
            (self.files['dna_s1_R1'],),
            (self.files['dna_s1_R2'],),
            (self.files['dna_s2_R1'],),
            (self.files['dna_s2_R2'],),
        ])
        data, rna_included = valid_input(sheet)
        self.assertEqual(len(data), 4)
        self.assertFalse(rna_included)

    def test_dna_and_rna_sheet(self):
        sheet = self._write_sheet(
            header=('DNA', 'RNA'),
            rows=[
                (self.files['dna_s1_R1'], self.files['rna_s1_R1']),
                (self.files['dna_s1_R2'], self.files['rna_s1_R2']),
                (self.files['dna_s2_R1'], self.files['rna_s2_R1']),
                (self.files['dna_s2_R2'], self.files['rna_s2_R2']),
            ],
        )
        data, rna_included = valid_input(sheet)
        self.assertEqual(len(data), 4)
        self.assertTrue(rna_included)

    def test_csv_delimiter(self):
        sheet = self._write_sheet(
            [(self.files['dna_s1_R1'],), (self.files['dna_s1_R2'],)],
            ext='.csv',
        )
        data, rna_included = valid_input(sheet)
        self.assertEqual(len(data), 2)
        self.assertFalse(rna_included)

    def test_paths_are_resolved_to_absolute(self):
        sheet = self._write_sheet([(self.files['dna_s1_R1'],)])
        data, _ = valid_input(sheet)
        self.assertTrue(os.path.isabs(data[0]['DNA']))

    def test_all_unique_files_across_multiple_samples_pass(self):
        # Sanity check: a sheet where every DNA/RNA file is genuinely
        # distinct must not trip the uniqueness check.
        sheet = self._write_sheet(
            header=('DNA', 'RNA'),
            rows=[
                (self.files['dna_s1_R1'], self.files['rna_s1_R1']),
                (self.files['dna_s1_R2'], self.files['rna_s1_R2']),
                (self.files['dna_s2_R1'], self.files['rna_s2_R1']),
                (self.files['dna_s2_R2'], self.files['rna_s2_R2']),
            ],
        )
        try:
            valid_input(sheet)
        except ArgumentTypeError as e:
            self.fail(f'valid_input() raised unexpectedly on a valid sheet: {e}')


class TestSheetLevelErrors(ValidInputTestCase):
    """Errors about the sheet file itself, not its row contents."""

    def test_missing_sheet_raises(self):
        missing = os.path.join(self.tmpdir.name, 'does_not_exist.tsv')
        with self.assertRaises(ArgumentTypeError):
            valid_input(missing)

    def test_missing_dna_column_raises(self):
        sheet = self._write_sheet(header=('RNA',), rows=[(self.files['rna_s1_R1'],)])
        with self.assertRaises(ArgumentTypeError):
            valid_input(sheet)


class TestRowLevelErrors(ValidInputTestCase):
    """Errors about individual row values."""

    def test_nonexistent_dna_file_raises(self):
        sheet = self._write_sheet([
            (os.path.join(self.tmpdir.name, 'nope_R1.fastq.gz'),),
        ])
        with self.assertRaises(ArgumentTypeError):
            valid_input(sheet)

    def test_nonexistent_rna_file_raises(self):
        sheet = self._write_sheet(
            header=('DNA', 'RNA'),
            rows=[(self.files['dna_s1_R1'], os.path.join(self.tmpdir.name, 'nope_R1.fastq.gz'))],
        )
        with self.assertRaises(ArgumentTypeError):
            valid_input(sheet)


class TestDuplicateFileDetection(ValidInputTestCase):
    """The sample-sheet uniqueness check: every file in the sheet must
    appear exactly once, whether reused as R1/R2 for different samples
    or reused across the DNA/RNA columns."""

    def test_same_file_reused_as_dna_in_two_rows_raises(self):
        sheet = self._write_sheet([
            (self.files['dna_s1_R1'],),
            (self.files['dna_s1_R2'],),
            (self.files['dna_s2_R1'],),
            (self.files['dna_s1_R1'],),  # reused from row 1
        ])
        with self.assertRaises(ArgumentTypeError) as ctx:
            valid_input(sheet)
        self.assertIn('used more than once', str(ctx.exception))

    def test_file_used_as_r1_for_one_sample_and_r2_for_another_raises(self):
        # The exact scenario this check exists for: the same fastq
        # accidentally used as R1 for sample 1 and R2 for sample 2.
        sheet = self._write_sheet([
            (self.files['dna_s1_R1'],),
            (self.files['dna_s1_R2'],),
            (self.files['dna_s2_R1'],),
            (self.files['dna_s1_R1'],),  # should have been dna_s2_R2
        ])
        with self.assertRaises(ArgumentTypeError):
            valid_input(sheet)

    def test_same_file_used_as_dna_and_rna_raises(self):
        sheet = self._write_sheet(
            header=('DNA', 'RNA'),
            rows=[
                (self.files['dna_s1_R1'], self.files['rna_s1_R1']),
                (self.files['dna_s1_R2'], self.files['dna_s1_R1']),  # reused as RNA
            ],
        )
        with self.assertRaises(ArgumentTypeError) as ctx:
            valid_input(sheet)
        self.assertIn('used more than once', str(ctx.exception))

    def test_error_message_identifies_both_locations(self):
        sheet = self._write_sheet([
            (self.files['dna_s1_R1'],),
            (self.files['dna_s1_R2'],),
            (self.files['dna_s2_R1'],),
            (self.files['dna_s1_R1'],),
        ])
        with self.assertRaises(ArgumentTypeError) as ctx:
            valid_input(sheet)
        message = str(ctx.exception)
        self.assertIn('sheet row 1', message)
        self.assertIn('sheet row 4', message)


class TestPairingChecks(ValidInputTestCase):
    """Every paired-end sample must resolve to exactly one R1 and one
    R2 file, identified via the same rename() transform metamorph uses
    to stage inputs (src.run.sym_safe() -> src.run.rename())."""

    def test_sample_missing_r2_raises(self):
        # dna_s2 only contributes an R1 -- dna_s1 (which has both mates)
        # is what makes this column "paired", so the missing R2 for s2
        # must be flagged rather than silently treated as single-end.
        sheet = self._write_sheet([
            (self.files['dna_s1_R1'],),
            (self.files['dna_s1_R2'],),
            (self.files['dna_s2_R1'],),
        ])
        with self.assertRaises(ArgumentTypeError) as ctx:
            valid_input(sheet)
        message = str(ctx.exception)
        self.assertIn('missing its R2 file', message)
        self.assertIn(self.files['dna_s2_R1'], message)

    def test_sample_missing_r1_raises(self):
        sheet = self._write_sheet([
            (self.files['dna_s1_R1'],),
            (self.files['dna_s1_R2'],),
            (self.files['dna_s2_R2'],),
        ])
        with self.assertRaises(ArgumentTypeError) as ctx:
            valid_input(sheet)
        message = str(ctx.exception)
        self.assertIn('missing its R1 file', message)
        self.assertIn(self.files['dna_s2_R2'], message)

    def test_all_single_end_column_is_not_flagged_as_missing_pairs(self):
        # No R2 anywhere in the column at all -- this is single-end
        # data, not a broken pair, and must not raise.
        # Single-end inputs still use the R1 naming convention -- there
        # is simply no corresponding R2 file anywhere in the sheet.
        single_end = os.path.join(self.tmpdir.name, 'sample_only_R1.fastq.gz')
        open(single_end, 'w').close()
        sheet = self._write_sheet([(single_end,)])
        try:
            valid_input(sheet)
        except ArgumentTypeError as e:
            self.fail(f'valid_input() raised on single-end-only data: {e}')

    def test_r1_and_r2_resolving_to_same_file_raises(self):
        # Different path strings that point at the same file on disk
        # (e.g. one path is a symlink to the other) must be rejected
        # even though the plain string-uniqueness check would miss it.
        real_file = self.files['dna_s1_R1']
        alias_dir = os.path.join(self.tmpdir.name, 'alias')
        os.makedirs(alias_dir, exist_ok=True)
        alias = os.path.join(alias_dir, 'dna_s1_R2.fastq.gz')
        os.symlink(real_file, alias)
        sheet = self._write_sheet([
            (real_file,),
            (alias,),
        ])
        with self.assertRaises(ArgumentTypeError) as ctx:
            valid_input(sheet)
        self.assertIn('resolve to the same file', str(ctx.exception))


class TestStagedSymlinkNameCollisions(ValidInputTestCase):
    """metamorph stages every input by renaming its basename to
    <sample>_R[12].fastq.gz and symlinking it into <output>/dna (or
    <output>/rna). Two different input files that rename to the same
    staged filename would silently collide/overwrite in that staging
    directory, so this must be rejected at sample-sheet validation time."""

    def _make_dir_with_file(self, dirname, filename):
        d = os.path.join(self.tmpdir.name, dirname)
        os.makedirs(d, exist_ok=True)
        path = os.path.join(d, filename)
        open(path, 'w').close()
        return path

    def test_same_basename_in_different_directories_raises(self):
        # Two genuinely different files (different directories) that
        # happen to share a basename both rename() to the identical
        # staged filename -- this is real data loss, not a copy-paste
        # duplicate, so the plain path-uniqueness check can't catch it.
        path_a = self._make_dir_with_file('runA', 'sample1_R1.fastq.gz')
        path_b = self._make_dir_with_file('runB', 'sample1_R1.fastq.gz')
        r2 = self._make_dir_with_file('runA', 'sample1_R2.fastq.gz')
        sheet = self._write_sheet([
            (path_a,),
            (r2,),
            (path_b,),
        ])
        with self.assertRaises(ArgumentTypeError) as ctx:
            valid_input(sheet)
        message = str(ctx.exception)
        self.assertIn('same renamed filename', message)
        self.assertIn(path_a, message)
        self.assertIn(path_b, message)

    def test_different_naming_conventions_collapsing_to_same_name_raises(self):
        # `_1.fastq.gz` and `.R1.fastq.gz` are both valid input naming
        # conventions that rename() normalizes to the same `_R1.fastq.gz`
        # staged name -- if two different source files with the same
        # sample prefix use different conventions, they still collide.
        path_a = self._make_dir_with_file('runA', 'sample2_1.fastq.gz')
        path_b = self._make_dir_with_file('runB', 'sample2.R1.fastq.gz')
        sheet = self._write_sheet([
            (path_a,),
            (path_b,),
        ])
        with self.assertRaises(ArgumentTypeError) as ctx:
            valid_input(sheet)
        self.assertIn('same renamed filename', str(ctx.exception))

    def test_dna_and_rna_columns_do_not_collide_with_each_other(self):
        # DNA and RNA are staged into separate directories (<output>/dna
        # vs <output>/rna), so a shared basename across the two columns
        # is not a staging collision -- it's covered (if at all) by the
        # separate cross-column path-uniqueness check, not this one.
        dna_path = self._make_dir_with_file('dna_run', 'shared_name_R1.fastq.gz')
        rna_path = self._make_dir_with_file('rna_run', 'shared_name_R1.fastq.gz')
        dna_r2 = self._make_dir_with_file('dna_run', 'shared_name_R2.fastq.gz')
        rna_r2 = self._make_dir_with_file('rna_run', 'shared_name_R2.fastq.gz')
        sheet = self._write_sheet(
            header=('DNA', 'RNA'),
            rows=[
                (dna_path, rna_path),
                (dna_r2, rna_r2),
            ],
        )
        try:
            valid_input(sheet)
        except ArgumentTypeError as e:
            self.fail(f'valid_input() raised on a cross-column basename match: {e}')


if __name__ == '__main__':
    unittest.main()
