import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from evogene.local_rescue import rescue_in_interval

SEQ='ATGCGTACCTGACATCGGATCAGTACGATCAGTTGCAAGCTAGCTACGACGTACGATGCTAGTACCGGTGATCGATGCTAGTCAGCTAGCTGACTAGCTAGTCGATACGATGCTGCA'

class LocalSyntenyRescueTests(unittest.TestCase):
    def test_short_genomic_homolog_is_rescued_in_syntenic_window(self):
        with TemporaryDirectory() as folder:
            p=Path(folder)
            genome=p/'ref.fna'; genome.write_text('>chr1\n'+'C'*200+SEQ+'G'*600+'\n')
            gff=p/'ref.gff3';gff.write_text('chr1\tX\tCDS\t701\t800\t.\t+\t0\tParent=rna1\n')
            result=rescue_in_interval(SEQ,genome,'chr1',170,450,gff)
        self.assertEqual(result['status'],'candidate_local_homology')
        self.assertGreater(result['hit']['query_coverage'],.9)
        self.assertEqual(result['hit']['coding_overlap'],'no_cds_overlap_not_proof_of_noncoding')

    def test_unavailable_gff_does_not_prove_noncoding(self):
        with TemporaryDirectory() as folder:
            p=Path(folder)/'ref.fna'
            p.write_text('>chr1\n'+'C'*200+SEQ+'G'*500+'\n')
            result=rescue_in_interval(SEQ,p,'chr1',170,450)
        self.assertEqual(result['status'],'candidate_local_homology')
        self.assertEqual(result['hit']['coding_overlap'],'unknown_no_gff')

    def test_unbounded_alignment_is_skipped(self):
        result=rescue_in_interval(SEQ,'nonexistent','chr1',0,100000)
        self.assertEqual(result['status'],'skipped_size')
