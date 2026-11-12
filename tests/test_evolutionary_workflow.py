"""Regression tests for multi-gene evolutionary workflow and real phylogeny."""
import json
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from Bio import Phylo

from evogene.phylogeny import infer_tree, collect_locus_sequences
from evogene.evolutionary_workflow import candidate_records, run_project


class TreeTests(unittest.TestCase):
    def test_real_dna_distance_tree_without_invented_support(self):
        seqs={
            "Drosophila_melanogaster":"ACGT"*35+"AAAA"*8,
            "Drosophila_simulans":"ACGT"*35+"AAAT"*8,
            "Drosophila_yakuba":"ACGT"*30+"TGCA"*13,
            "Drosophila_ananassae":"TGCA"*30+"ACGT"*13,
        }
        r=infer_tree(seqs,try_external=False)
        self.assertEqual(r["status"],"inferred_exploratory")
        self.assertIsNone(r["support_type"])
        self.assertIn("Jaccard",r["method"])
        tree=Phylo.read(StringIO(r["newick"]),"newick")
        self.assertEqual({x.name for x in tree.get_terminals()},set(seqs))
        self.assertTrue(any("NOT a validated" in x for x in r["warnings"]))
        self.assertFalse(r["rooted"])

    def test_fewer_than_three_sequences_is_not_fake_tree(self):
        r=infer_tree({"a":"ACGT"*12,"b":"ACGT"*12},try_external=False)
        self.assertEqual(r["status"],"not_inferred")
        self.assertIsNone(r["newick"])

    def test_homology_selection_records_ambiguous_loci(self):
        with TemporaryDirectory() as folder:
            genome=Path(folder)/"dsim.fa"
            genome.write_text(">chr1\n"+"A"*300+"\n")
            q="A"*110
            hit={"species":"Drosophila simulans","chrom":"chr1",
                 "start":10,"end":120,"strand":"+",
                 "query_coverage":1.0,"identity_approx":.99,
                 "mapq":0,"coding_overlap":"overlaps_cds"}
            report={"target_species":"Drosophila melanogaster",
                    "requested_relative_species":["Drosophila simulans"],
                    "ncbi_assemblies":[{"species":"Drosophila simulans",
                                       "accession":"GCF_TEST.1","fasta":str(genome)}],
                    "relatives":{"Drosophila simulans":[hit,dict(hit,start=130,end=240)]}}
            seqs, details=collect_locus_sequences(q,report)
            self.assertEqual(len(seqs),1)
            self.assertIn("ambiguous",details[1]["reason"])
            report["relatives"]["Drosophila simulans"][0]["mapq"]=60
            seqs, details=collect_locus_sequences(q,report)
            self.assertEqual(len(seqs),2)
            self.assertEqual(details[-1]["genome_accession"],"GCF_TEST.1")


class CandidateAndWorkflowTests(unittest.TestCase):
    def test_candidate_parser_rejects_duplicates_and_invalid_dna(self):
        self.assertEqual(len(candidate_records(">g1\n"+"ACGT"*20+"\n>g2\n"+"ACGT"*20)),2)
        with self.assertRaises(ValueError):
            candidate_records(">g1\n"+"ACGT"*20+"\n>g1\n"+"ACGT"*20)
        with self.assertRaises(ValueError):
            candidate_records(">g1\n"+"ACTX"*20)

    @patch("evogene.evolutionary_workflow.infer_report_tree")
    @patch("evogene.evolutionary_workflow.ncbi_auto.analyze")
    def test_audited_multi_gene_and_downloads(self,analyze,phylo):
        def fake_auto(path, focal, gene, **kwargs):
            kwargs["progress"]({
                "stage":"dna_alignment","state":"completed",
                "detail":"1 candidate hit","species":"Drosophila simulans",
                "timestamp_utc":"2026-11-13T00:00:00+00:00"
            })
            if gene=="bad":
                raise RuntimeError("Failed analysis in controlled test")
            return {
                "gene_id":gene,"requested_relative_species":["Drosophila simulans",
                                                              "Drosophila yakuba"],
                "relatives":{"Drosophila simulans":[]},
                "ncbi_assemblies":[],
                "mechanism_report":{"tests":[{
                    "mechanism":"dna_duplication","status":"candidate",
                    "signals":["Additional mapped locus"],"observations":{"copy_count":2},
                    "limitations":["Synteny not verified"]
                }]},
                "coding_species":[],
            }
        analyze.side_effect=fake_auto
        phylo.return_value={
            "status":"inferred_exploratory","newick":"(a:0.1,b:0.2,c:0.3);",
            "method":"EXPLORATORY_NJ_from_7mer_Jaccard","taxa":3,
            "support_type":None}
        dna="ATG"+"AAA"*30+"TAA"
        with TemporaryDirectory() as folder:
            report=run_project(
                ">good\n"+dna+"\n>bad\n"+dna,
                "Drosophila melanogaster",
                ["Drosophila simulans","Drosophila yakuba"],
                output_folder=folder)
            self.assertEqual(report["summary"]["requested"],2)
            self.assertEqual(report["summary"]["completed"],1)
            self.assertEqual(report["summary"]["failed"],1)
            self.assertEqual(report["summary"]["tree_computed"],1)
            self.assertTrue(report["input_sha256"])
            events=report["audit_log"]
            self.assertTrue(all("timestamp_utc" in e for e in events))
            self.assertTrue(any(e["state"]=="failed" for e in events))
            self.assertEqual(report["genes"]["good"]["hypotheses"]["most_likely"],
                             "not_established")
            self.assertEqual(report["genes"]["good"]["hypotheses"]["candidates"],
                             ["duplication"])
            self.assertTrue((Path(folder)/"good"/"gene_tree.nwk").exists())
            self.assertTrue((Path(folder)/"audit.jsonl").exists())
            loaded=json.loads((Path(folder)/"run_report.json").read_text())
            self.assertEqual(loaded["summary"]["requested"],2)
        self.assertEqual(analyze.call_count,2)

if __name__=="__main__":
    unittest.main()
