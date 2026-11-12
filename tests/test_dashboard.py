"""UI smoke tests for redesigned EvoGene research dashboard."""
import unittest

from evogene.dashboard import _bundle


class DashboardTests(unittest.TestCase):
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
