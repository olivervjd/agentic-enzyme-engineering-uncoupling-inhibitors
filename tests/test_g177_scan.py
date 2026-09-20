"""Software integrity tests; synthetic fixtures are not biological evidence."""
import json
from pathlib import Path
import tempfile
import unittest

from herbicide_desensitization_agent.examples.g177_scan import AMINO_ACIDS, mutate, replace_query
from herbicide_desensitization_agent.examples.g177_apo_foldseek import alignment_counts, pack


class G177SequenceIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.sequence = 'A' * 176 + 'G' + 'A' * 343

    def test_all_twenty_variants_change_only_canonical_177(self):
        self.assertEqual(len(AMINO_ACIDS), 20)
        for aa in AMINO_ACIDS:
            result = mutate(self.sequence, aa)
            self.assertEqual(result[:176], self.sequence[:176])
            self.assertEqual(result[177:], self.sequence[177:])
            self.assertEqual(result[76:][100], aa)
            self.assertEqual(len(result[76:]), 444)

    def test_msa_replacement_preserves_homolog_bytes(self):
        tail = '>homolog1\nAAAAaaA--\n>homolog2\nGGGG\n'
        original = '>query\n' + self.sequence[76:] + '\n' + tail
        result = replace_query(original, self.sequence[76:], mutate(self.sequence, 'W')[76:])
        self.assertTrue(result.endswith(tail))
        self.assertEqual(result.splitlines()[1][100], 'W')

    def test_wrong_reference_or_invalid_substitution_is_rejected(self):
        for sequence, aa in [(self.sequence[:176]+'A'+self.sequence[177:], 'W'),
                             (self.sequence[:-1], 'W'), (self.sequence, 'X')]:
            with self.assertRaises(ValueError):
                mutate(sequence, aa)
        with self.assertRaises(ValueError):
            replace_query('>q\nWRONG\n>h\nAAAA\n', self.sequence[76:], self.sequence[76:])
        with self.assertRaises(ValueError):
            replace_query('>q\n'+self.sequence[76:]+'\n', self.sequence[76:], self.sequence[76:])


class FoldseekAlignmentTests(unittest.TestCase):
    def row(self, q='AC-DE', t='A-CDE', **changes):
        row = dict(qaln=q, taln=t, alnlen='5', qlen='4', tlen='4', qcov='1', tcov='1')
        row.update(changes)
        return row

    def test_gaps_do_not_count_as_paired_residues(self):
        counts = alignment_counts(self.row())
        self.assertEqual(counts, dict(alignment_columns=5, paired_residues=3, paired_coverage=.75))

    def test_self_alignment_pairs_all_residues(self):
        counts = alignment_counts(self.row(q='ACDE', t='ACDE', alnlen='4'))
        self.assertEqual(counts['paired_residues'], 4)
        self.assertEqual(counts['paired_coverage'], 1)

    def test_bad_alignment_or_coverage_is_rejected(self):
        for row in [self.row(alnlen='4'), self.row(q='ACDE'), self.row(qcov='.5'),
                    self.row(q='AC-DE', t='AC-DE'), self.row(qlen='2')]:
            with self.subTest(row=row), self.assertRaises(ValueError):
                alignment_counts(row)


class G177PackagingTests(unittest.TestCase):
    def test_incomplete_campaign_cannot_be_packaged(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'results.json').write_text(json.dumps(dict(complete=False, predictions=59)))
            (root/'validation.json').write_text(json.dumps(dict(status='PASSED')))
            with self.assertRaisesRegex(ValueError, 'incomplete'):
                pack(root)

    def test_changed_coordinates_cannot_be_packaged(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'model.cif').write_text('changed coordinates')
            (root/'results.json').write_text(json.dumps(dict(complete=True, predictions=60,
                artifacts=[dict(path='model.cif', sha256='0'*64)])))
            (root/'validation.json').write_text(json.dumps(dict(status='PASSED')))
            with self.assertRaisesRegex(ValueError, 'hash changed'):
                pack(root)


if __name__ == '__main__':
    unittest.main()
