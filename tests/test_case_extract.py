"""Regression tests for published three-section extract, species panel, and annotations."""
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from evogene.case_data import (
    sections, index_file, analyze, export_gene_fasta, header_info, compare
)
from evogene.ncbi_auto import PANEL, SPECIES_BY_CODE, analyze as auto_analyze
from evogene.analysis_notes import annotations_for_report


EXAMPLE = """Section I. Raw hits extracted from progressive cactus alignments
raw sequence for dmel_FBgn0000001
>dmel|FBgn0000001|score_1E10|shift_No|1
ATGAAACCCGGGTAA
>dsim|NW_000001.1+:100:115
ATGAAACCCGGGTAA
raw sequence for dmel_FBgn0000002
>dmel|FBgn0000002|score_1E10|shift_No|1
ATGAAATAG
Section II. Spliced align by Spaln
spaln for dmel_FBgn0000001
>dmel|FBgn0000001|score_1E10|shift_No|1
ATGAAACCCGGGTAA
>dsim|NW_000001.1+:100:115|score_55.0|shift_Yes|1
ATGAAATAA
spaln for dmel_FBgn0000002
>dmel|FBgn0000002|score_1E10|shift_No|1
ATGAAATAG
Section III. Spliced align by Genewise
genewise for dmel_FBgn0000001
>dmel|FBgn0000001|score_1E10|shift_No|1
ATGAAACCCGGGTAA
>dsim|NW_000001.1+:100:115|score_8.2|shift_No|1
ATGAAGTAA
genewise for dmel_FBgn0000002
>dmel|FBgn0000002|score_1E10|shift_No|1
ATGAAATAG
"""


class CaseExtractTests(unittest.TestCase):
    def setUp(self):
        self.folder = TemporaryDirectory()
        self.path = Path(self.folder.name) / "case.fasta"
        self.path.write_text(EXAMPLE, encoding="utf8")

    def tearDown(self):
        self.folder.cleanup()

    def test_index_uses_all_sections_and_genes(self):
        result = index_file(self.path)
        self.assertEqual(len(result["genes"]), 2)
        self.assertEqual(result["record_counts"], {"raw": 3, "spaln": 3, "genewise": 3})
        self.assertEqual(set(result["species"]), {"dmel", "dsim"})

    def test_method_boundaries_do_not_corrupt_dna(self):
        for method in ("raw", "spaln", "genewise"):
            gene, records = next(sections(self.path, method, "FBgn0000001"))
            self.assertEqual(gene, "dmel_FBgn0000001")
            self.assertEqual(len(records), 2)
            self.assertEqual(records[0][1], "ATGAAACCCGGGTAA")
            self.assertNotIn("SPALN", records[0][1])
            self.assertNotIn("GENEWISE", records[0][1])
        self.assertEqual(next(sections(self.path, "raw", "FBgn0000001"))[1][1][1],
                         "ATGAAACCCGGGTAA")
        self.assertEqual(next(sections(self.path, "spaln", "FBgn0000001"))[1][1][1],
                         "ATGAAATAA")

    def test_annotated_and_coordinate_headers(self):
        info = header_info("dsim|NW_000001.1+:100:115|score_55.0|shift_Yes|1")
        self.assertEqual(info["header_kind"], "genomic_coordinates")
        self.assertEqual(info["prediction_score_raw"], "55.0")
        self.assertEqual(info["shift_flag"], "Yes")

    def test_case_gene_analysis_does_not_claim_origin(self):
        result = analyze(self.path, "FBgn0000001", method="spaln",
                         selected_species=["dsim"])
        self.assertEqual(result["method"], "spaln")
        self.assertEqual(len(result["hits"]), 1)
        self.assertEqual(len(result["method_comparison"]), 3)
        self.assertEqual(result["phylogeny"]["status"], "not_assessed")
        self.assertNotIn("mechanism_report", result)

    def test_export_unaligned_fasta(self):
        text = export_gene_fasta(self.path, "FBgn0000001", "genewise",
                                 species=["dmel", "dsim"])
        self.assertEqual(text.count(">"), 2)
        self.assertIn("ATGAAGTAA", text)
        self.assertNotIn("genewise for", text)

    def test_quadratic_guard(self):
        self.assertEqual(compare("A" * 5000, "T" * 5000)["comparison_status"],
                         "skipped_quadratic_cost")


class PanelTests(unittest.TestCase):
    def test_species_panel_not_limited_to_five(self):
        self.assertEqual(len(SPECIES_BY_CODE), 20)
        self.assertEqual(len(PANEL["Drosophila melanogaster"]), 19)
        self.assertIn("Bactrocera dorsalis", PANEL["Drosophila melanogaster"])
        self.assertNotIn("Drosophila melanogaster", PANEL["Drosophila melanogaster"])

    def test_invalid_species_and_duplicate_rejected(self):
        with self.assertRaises(ValueError):
            auto_analyze("missing.fasta", "Drosophila melanogaster", "g",
                         relative_species=["Drosophila yakuba", "Drosophila yakuba"])
        with self.assertRaises(ValueError):
            auto_analyze("missing.fasta", "Drosophila melanogaster", "g",
                         relative_species=["Drosophila melanogaster"])


class AnnotationTests(unittest.TestCase):
    def test_missing_remains_not_assessed(self):
        result = annotations_for_report({})
        states = {x["key"]: x["state"] for x in result}
        self.assertEqual(states["phylogeny"], "not_assessed")
        self.assertEqual(states["protein_homology"], "not_assessed")

    def test_case_method_notes(self):
        result = annotations_for_report({"input_type":
                 "published_preselected_hits_not_whole_genome_or_msa",
                 "method_comparison": [{"method": "raw"}]})
        keys = {x["key"] for x in result}
        self.assertTrue({"case_raw", "case_spaln", "case_genewise"} <= keys)
        self.assertNotIn("download", keys)


if __name__ == "__main__":
    unittest.main()
