"""UI smoke validation: Altair specs and supplied Newick render without Streamlit server."""
from contextlib import nullcontext
from unittest.mock import patch
import unittest
from evogene.visualization import render_report_visualizations


class DummyUpload:
    def getvalue(self):
        return b"((dmel:0.1,dsim:0.2)95:0.1,dyak:0.3);"


class VisualSmokeTests(unittest.TestCase):
    def test_visual_panels_build_valid_altair_specs(self):
        import streamlit as st
        data = {
            "gene_id": "synthetic_gene",
            "relatives": {"dsim": [{
                "query_start": 0, "query_end": 100, "query_coverage": .9,
                "identity_approx": .8, "chrom": "2L", "start": 100,
                "end": 200, "strand": "+", "coding_overlap": "overlaps_cds",
                "alignment_insertion_nt": 3, "alignment_deletion_nt": 0,
            }]},
            "synteny_evidence": {"dsim": {
                "status": "flanks_collinear", "locus_status": "homolog_overlaps_cds",
                "left_anchor": {"chrom": "2L", "start": 0, "end": 90, "strand": "+"},
                "right_anchor": {"chrom": "2L", "start": 220, "end": 300, "strand": "+"},
                "syntenic_hits": [{"chrom": "2L", "start": 100, "end": 200, "strand": "+"}],
            }},
            "query_orfs": [{
                "strand": "+", "frame": 0, "start": 0, "end": 99, "aa_length": 32
            }],
            "mechanism_report": {
                "tests": [{"label": "De novo", "mechanism": "de_novo",
                           "status": "partial_signal", "signal_count": 1}],
                "completeness": {"assessed": 1, "total": 5,
                                 "not_assessed": ["fusion"]},
            },
        }
        def validate(chart, **kwargs):
            result = chart.to_dict(validate=True)
            self.assertIn("mark", result)
        names = ("header", "caption", "dataframe", "info", "warning",
                 "metric", "markdown", "error", "download_button")
        from contextlib import ExitStack
        with ExitStack() as stack:
            for name in names:
                stack.enter_context(patch.object(st, name))
            stack.enter_context(patch.object(st, "tabs",
                                             return_value=[nullcontext() for _ in range(6)]))
            stack.enter_context(patch.object(st, "file_uploader",
                                             side_effect=[None, DummyUpload()]))
            charts = stack.enter_context(patch.object(st, "altair_chart",
                                                       side_effect=validate))
            tree = stack.enter_context(patch.object(st, "graphviz_chart"))
            render_report_visualizations(data)
        self.assertGreaterEqual(charts.call_count, 4)
        self.assertEqual(tree.call_count, 1)
        self.assertIn("dmel", tree.call_args.args[0])


if __name__ == "__main__":
    unittest.main()
