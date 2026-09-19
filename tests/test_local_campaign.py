import unittest

from herbicide_desensitization_agent.examples.epsps_local_campaign import CONDITIONS, mutate


class LocalCampaignTests(unittest.TestCase):
    def test_wt_and_exact_single_substitution(self):
        sequence = "A" * 520
        self.assertEqual(mutate(sequence, "WT"), sequence)
        result = mutate(sequence, "A288V")
        self.assertEqual(result[287], "V")
        self.assertEqual(sum(a != b for a, b in zip(sequence, result)), 1)

    def test_mutations_outside_domain_or_wrong_identity_fail(self):
        for mutation in ("A76V", "M288L", "A521V", "A288A"):
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                mutate("A" * 520, mutation)

    def test_native_substrate_affinities_keep_same_complex_context(self):
        self.assertEqual(CONDITIONS["native"], ("pep", "B"))
        self.assertEqual(CONDITIONS["native_s3p"], ("pep", "C"))
        self.assertEqual(CONDITIONS["herbicide"], ("glyphosate", "B"))
