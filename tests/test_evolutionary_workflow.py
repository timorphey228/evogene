"""Regression tests for multi-gene evolutionary workflow and real phylogeny."""
import json
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from Bio import Phylo

from evogene.phylogeny import infer_tree, collect_locus_sequences, alignment_quality
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

    def test_real_aligned_ml_output_is_preserved(self):
        from types import SimpleNamespace
        from Bio.Align import MultipleSeqAlignment
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        rows={
            "a":"ATG"*34+"AAA"*8,
            "b":"ATG"*34+"AAG"*8,
            "c":"ATG"*30+"CCC"*12,
            "d":"ATG"*28+"TAA"*14
        }
        def executable(name):
            return "/usr/bin/"+name
        def fake_run(cmd, **kwargs):
            if "mafft" in cmd[0]:
                return SimpleNamespace(returncode=0,stderr="MAFFT successful",
                    stdout="".join(f">{name}\\n{seq}\\n" for name,seq in rows.items()))
            prefix=cmd[cmd.index("--prefix")+1]
            Path(prefix+".treefile").write_text(
                "((a:0.01,b:0.02)99/96:0.03,(c:0.04,d:0.05)98/94:0.06);")
            return SimpleNamespace(returncode=0,stderr="IQ-TREE successful",stdout="")
        with (patch("evogene.phylogeny.find_tool",side_effect=executable),
              patch("evogene.phylogeny.subprocess.run",side_effect=fake_run)):
            result=infer_tree(rows,try_external=True)
        self.assertEqual(result["status"],"inferred")
        self.assertIn("maximum_likelihood",result["method"])
        self.assertEqual(result["support_type"],
                         "UFBoot 1000 + SH-aLRT 1000 (IQ-TREE)")
        self.assertEqual(result["alignment_quality"]["aligned_columns"],126)
        self.assertEqual(len(result["diagnostics"]),2)
        self.assertEqual(len(list(Phylo.read(StringIO(result["newick"]),"newick").get_terminals())),4)

    def test_alignment_quality_flags_uninformative_sites(self):
        from Bio.Align import MultipleSeqAlignment
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        alignment=MultipleSeqAlignment([
            SeqRecord(Seq("A"*70),id=name) for name in ("a","b","c","d")])
        q=alignment_quality(alignment)
        self.assertTrue(q["low_information"])
        self.assertEqual(q["variable_sites"],0)

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
            runfolder=Path(folder)/report["run_id"]
            self.assertTrue((runfolder/"good"/"gene_tree.nwk").exists())
            self.assertTrue((runfolder/"audit.jsonl").exists())
            loaded=json.loads((runfolder/"run_report.json").read_text())
            self.assertEqual(loaded["summary"]["requested"],2)
        self.assertEqual(analyze.call_count,2)

if __name__=="__main__":
    unittest.main()
