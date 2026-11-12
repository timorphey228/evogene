"""Regression tests: user-facing two-FASTA multi-gene scenario."""
import unittest

from evogene.user_sequences import analyze_user_fastas


FOCAL = """>geneA|dmel
ATGAAACCCGGGTTTAAACCCGGGTTTTAA
>geneB|dmel
ATGAAGCCCGGGTATCCCGGGAAACCCCAA
"""
COMPARE = """>geneA|dsim|locus1
ATGAAACCCGGGTTTAAACCCGGGTTTTAA
>geneA|dyak|locus2
ATGAAACCCGGGTTTAAACCCGGCTTTTAA
>geneB|dsim|locus3
ATGAAGCCCGGGTATCCCGGGAAACCCCAA
>unknown|dyak|locus4
ATGAAGCCCGGGTATCCCGGGAAACCCCAA
"""


class UserFastaTests(unittest.TestCase):
    def test_multiple_targets_with_comparative_species(self):
        result = analyze_user_fastas(FOCAL, COMPARE)
        self.assertEqual(result["gene_count"], 2)
        self.assertEqual(result["summary"]["failed"], 0)
        self.assertEqual(len(result["summary"]["unmatched_comparison_headers"]), 1)
        self.assertEqual(len(result["source"]["sha256"]), 64)
        self.assertEqual(len(result["results"]), 2)
        for entry in result["results"]:
            self.assertEqual(entry["conclusion"]["category"], "insufficient_data")
            self.assertEqual(len(entry["case_requirements"]), 6)
            self.assertEqual(entry["phylogeny"]["status"], "not_assessed")
            self.assertEqual(len(entry["evidence"]), 2)
            self.assertIn("cross_method_comparison", [x["stage"] for x in entry["analysis_steps"]])
        self.assertEqual(len(result["raw_reports"]["geneA"]["hits"]), 2)

    def test_missing_comparisons_remain_unknown(self):
        result = analyze_user_fastas(FOCAL, """>unmatched|dsim
ATGAAACCCGGGTTTAAACCCGGGTTTTAA
""")
        self.assertEqual(result["summary"]["completed"], 2)
        self.assertEqual(result["raw_reports"]["geneA"]["summary"]["record_count"], 0)
        self.assertEqual(result["results"][0]["conclusion"]["category"], "insufficient_data")

    def test_validation_does_not_silently_accept_invalid_data(self):
        with self.assertRaises(ValueError):
            analyze_user_fastas(""">geneA
AAACTG
>geneA
AAAGGG
""", COMPARE)
        with self.assertRaises(ValueError):
            analyze_user_fastas(FOCAL, """>no_species
ATGAAAACCC
""")
        with self.assertRaises(ValueError):
            analyze_user_fastas(FOCAL, """>geneA|dsim
ATGC--GTAC
""")


if __name__ == "__main__":
    unittest.main()
