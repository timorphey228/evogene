"""Regression tests: all visualized results must come from genuine evidence."""
import unittest

from evogene.visualization import (
    alignment_blocks, synteny_blocks, orf_blocks, mechanism_rows,
    parse_aligned_fasta, alignment_cells, newick_to_dot, alignment_gap_rows,
)


class VisualizationTests(unittest.TestCase):
    def test_alignment_only_from_real_query_coordinates(self):
        report = {"relatives": {
            "D. simulans": [{
                "query_start": 3, "query_end": 45, "query_coverage": 0.7,
                "identity_approx": 0.8, "chrom": "2L", "start": 10,
                "end": 52, "strand": "+", "coding_overlap": "overlaps_cds",
            }],
            "D. yakuba": [{"chrom": "2L"}],
        }}
        result = alignment_blocks(report)
        self.assertEqual(len(result), 1)
        self.assertEqual((result[0]["start"], result[0]["end"]), (3, 45))
        self.assertEqual(result[0]["identity"], 80)
        self.assertEqual(result[0]["annotation"], "overlaps_cds")

    def test_synteny_tracks_require_both_anchors(self):
        anchors = {"status": "flanks_collinear", "locus_status": "homolog_no_cds_overlap",
                   "left_anchor": {"chrom": "2L", "start": 100, "end": 200, "strand": "+"},
                   "right_anchor": {"chrom": "2L", "start": 500, "end": 600, "strand": "+"},
                   "syntenic_hits": [{"chrom": "2L", "start": 280, "end": 320, "strand": "+"}]}
        rows = synteny_blocks({"synteny_evidence": {
            "a": anchors, "b": {"status": "not_established"},
        }})
        self.assertEqual(len(rows), 3)
        self.assertEqual([r["feature"] for r in rows],
                         ["left flank", "right flank", "candidate DNA homolog"])
        self.assertEqual(rows[-1]["start"], 180)

    def test_orf_is_only_recorded_if_measured(self):
        self.assertEqual(orf_blocks({}), [])
        self.assertEqual(orf_blocks({"query_orfs": [{
            "strand": "-", "frame": 2, "start": 10, "end": 100, "aa_length": 29
        }]})[0]["frame"], "-3")

    def test_mechanism_data_has_no_probabilities(self):
        result = mechanism_rows({"mechanism_report": {"tests": [{
            "mechanism": "de_novo", "label": "De novo",
            "status": "partial_signal", "signal_count": 1
        }]}})
        self.assertEqual(result[0]["status"], "partial_signal")
        self.assertNotIn("support_percent", result[0])

    def test_aligned_fasta_residue_grid(self):
        aligned = parse_aligned_fasta(">one\nAT-G\n>two\nACGG\n")
        self.assertEqual(len(alignment_cells(aligned, start=0, width=2)), 4)
        self.assertEqual(alignment_cells(aligned)[2]["kind"], "gap")
        protein = parse_aligned_fasta(">a\nM-EPQ\n>b\nMIKPQ\n")
        self.assertEqual(len(protein), 2)
        with self.assertRaises(ValueError):
            parse_aligned_fasta(">a\nACT\n>b\nACTG\n")

    def test_paf_cigar_insertion_deletion(self):
        from evogene.genomes import cigar_gap_summary
        self.assertEqual(cigar_gap_summary([""] * 12 + ["cg:Z:12M3I2D10M"]),
                         {"alignment_insertion_nt": 3, "alignment_deletion_nt": 2})
        self.assertIsNone(cigar_gap_summary([""] * 12)["alignment_deletion_nt"])
        rows = alignment_gap_rows({"relatives": {"a": [{
            "alignment_insertion_nt": 3, "alignment_deletion_nt": 2,
            "chrom": "X"
        }]}})
        self.assertEqual(len(rows), 2)

    def test_newick_tree_only_from_supplied_data(self):
        dot = newick_to_dot("((dmel:0.1,dsim:0.2)95:0.1,dyak:0.3);")
        self.assertIn("dmel", dot)
        self.assertIn("95", dot)
        self.assertIn("dyak", dot)
        self.assertIn('color="#FFFFFF"', dot)
        self.assertIn('fontcolor="#FFFFFF"', dot)
        self.assertIn('fillcolor="transparent"', dot)
        with self.assertRaises(ValueError):
            newick_to_dot(";")

if __name__ == "__main__":
    unittest.main()
