import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from evogene.protein_screen import query_orfs, protein_homology
from evogene.mechanism_v3 import classify

class ProteinHomologyTests(unittest.TestCase):
    def test_six_frame_orf_scan(self):
        dna='ATG'+('GCTGATGAACAA'*18)+'TAA'
        seqs=query_orfs(dna)
        self.assertTrue(seqs)
        self.assertGreaterEqual(seqs[0]['length_aa'],55)

    def test_protein_parent_can_support_retroposition(self):
        dna='ATG'+('GCTGATGAACAA'*18)+'TAA'
        peptide=query_orfs(dna)[0]['aa']
        with TemporaryDirectory() as folder:
            folder=Path(folder)
            fasta=folder/'proteins.faa'
            fasta.write_text('>XP_parent.1 putative_parent\n'+peptide+'\n')
            gff=folder/'genes.gff3'
            gff.write_text(
                '##gff-version 3\n'
                'chr2\tE\tCDS\t101\t300\t.\t+\t0\tParent=rna-child;protein_id=XP_child.1\n'
                'chrX\tE\tCDS\t1001\t1100\t.\t+\t0\tParent=rna-parent;protein_id=XP_parent.1\n'
                'chrX\tE\tCDS\t1301\t1400\t.\t+\t0\tParent=rna-parent;protein_id=XP_parent.1\n')
            focal={'chrom':'chr2','start':100,'end':300}
            evidence=protein_homology(dna,fasta,gff,focal)
        self.assertEqual(evidence['status'],'searched')
        self.assertTrue(evidence['candidates'])
        top=evidence['candidates'][0]
        self.assertEqual(top['gene_model']['cds_exons'],2)
        self.assertFalse(top['overlaps_target'])
        report={'target_locus':focal,'focal_reference':{'accession':'test'},
                'retrocopy_evidence':{'target_model':{'status':'annotated','exon_count':1}},
                'protein_evidence':evidence}
        m=classify(report)
        self.assertEqual(m['leading_candidate'],'retroposition')

    def test_no_protein_access(self):
        self.assertFalse(query_orfs('ATGAAATAG'))
