"""UI smoke tests for redesigned EvoGene research dashboard."""
import unittest

from evogene.dashboard import _bundle, _hypothesis_summary_rows, _style


class DashboardTests(unittest.TestCase):
    def test_comparison_species_labels_have_contrasting_black_text(self):
        from unittest.mock import patch
        with patch("evogene.dashboard.st.markdown") as render:
            _style()
        css=render.call_args.args[0]
        self.assertIn('[data-testid="stMultiSelect"] [data-baseweb="tag"]',css)
        self.assertIn("background-color: var(--evo-lime) !important",css)
        self.assertIn("color: #0F1010 !important",css)

    def test_zip_exports_tree_and_timestamps(self):
        from io import BytesIO
        from zipfile import ZipFile
        report={
            "run_id":"synthetic01", "audit_log":[{
                "gene_id":"geneA","stage":"phylogeny",
                "timestamp_utc":"2026-11-13T01:00:00+00:00",
                "state":"completed","detail":"tree"}],
            "genes":{"geneA":{
                "phylogeny":{"newick":"(focal:0.1,relative1:0.2,relative2:0.3);",
                             "alignment_fasta":">focal\nACGT\n>relative1\nACGA\n>relative2\nATGA\n"}
            }}
        }
        with ZipFile(BytesIO(_bundle(report))) as archive:
            names=archive.namelist()
            self.assertIn("run_report.json",names)
            self.assertIn("audit.jsonl",names)
            self.assertIn("geneA/gene_tree.nwk",names)
            self.assertIn("geneA/msa.fasta",names)
            self.assertIn("timestamp_utc",archive.read("audit.jsonl").decode())

    def test_hypothesis_summary_counts_and_qualitative_confidence(self):
        hypotheses = [
            {"label": "Дупликация", "status": "candidate",
             "arguments": ["Дополнительная копия", "Белковый гомолог"]},
            {"label": "De novo", "status": "partial_signal",
             "arguments": ["Синтенный некодирующий участок"]},
            {"label": "Ретропозиция", "status": "no_signal", "arguments": []},
            {"label": "Транспозон", "status": "not_assessed", "arguments": []},
            {"label": "Слияние", "status": "failed", "arguments": None},
        ]
        rows = _hypothesis_summary_rows(hypotheses)
        self.assertEqual(len(rows), 5)
        self.assertEqual([r["Поддерживающих признаков"] for r in rows],
                         [2, 1, 0, 0, 0])
        self.assertEqual(rows[0]["Уверенность системы"], "Предварительная")
        self.assertIn("Низкая", rows[1]["Уверенность системы"])
        self.assertEqual(rows[3]["Уверенность системы"], "Не определена")
        self.assertTrue(all("%" not in str(r["Уверенность системы"]) for r in rows))

    def test_initial_app_renders_without_uploaded_data(self):
        try:
            from streamlit.testing.v1 import AppTest
        except ImportError:
            self.skipTest("Streamlit AppTest not installed")
        from pathlib import Path
        app=AppTest.from_file(str(Path(__file__).resolve().parents[1]/"auto_app.py")).run(timeout=25)
        self.assertFalse(app.exception)
        self.assertTrue(any("EvoGene" in msg.value for msg in app.markdown))


if __name__=="__main__":
    unittest.main()
