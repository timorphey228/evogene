import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from evogene.mechanism_v3 import classify
from evogene.gene_models import compare_retrocopy, models_for_loci

def row(chrom='2L',start=100,end=200,coverage=.9,identity=.9):
    return {'chrom':chrom,'start':start,'end':end,
            'query_coverage':coverage,'identity_approx':identity}

def result(report,key):
    return next(v for v in classify(report)['options'] if v['mechanism']==key)

class ObjectiveClassifierTests(unittest.TestCase):
    def test_no_forced_100_percent_unknown(self):
        r={'relatives':{'simulans':[dict(row(),coding_overlap='overlaps_cds')],
                        'sechellia':[dict(row(),coding_overlap='overlaps_cds')]}}
        m=classify(r)
        self.assertNotIn('support_percent',m['options'][0])
        self.assertEqual(m['candidate_mechanisms'],[])
        self.assertIsNone(m['leading_candidate'])
        self.assertEqual(len(m['conservation']['cds_overlap_species']),2)
        self.assertEqual(m['completeness']['assessed'],1)

    def test_denovo_one_syntenic_no_cds_is_partial_signal(self):
        r={'synteny_evidence':{'yakuba':{'status':'flanks_collinear',
                  'locus_status':'homolog_no_cds_overlap'}}}
        self.assertEqual(result(r,'de_novo')['status'],'partial_signal')

    def test_denovo_not_inferred_from_missing_hits(self):
        r={'relatives':{'yakuba':[]},'synteny_evidence':{'yakuba':{'status':'not_established'}}}
        self.assertEqual(result(r,'de_novo')['status'],'no_signal')

    def test_duplication_signal_independent_of_retro(self):
        r={'focal_reference':{'accession':'X'},'target_locus':row(),
           'extra_target_copies':[row(chrom='3R')],
           'retrocopy_evidence':{'target_model':{'status':'annotated','exon_count':1},
              'parent_candidates':[{'parent_multi_exon_and_target_single_exon':True}]}}
        m=classify(r)
        self.assertIn('dna_duplication',m['candidate_mechanisms'])
        self.assertIn('retroposition',m['candidate_mechanisms'])
        self.assertIsNone(m['leading_candidate'])

    def test_retrogen_structure(self):
        r={'focal_reference':{'accession':'X'},'target_locus':row(),
           'protein_evidence':{'status':'searched','candidates':[]},
           'retrocopy_evidence':{'target_model':{'status':'annotated','exon_count':1},
              'parent_candidates':[{'parent_multi_exon_and_target_single_exon':True}]}}
        self.assertEqual(result(r,'retroposition')['status'],'candidate')

    def test_protein_parent_requires_annotated_target(self):
        protein={'status':'searched','candidates':[{
          'gene_model':{'cds_exons':4,'chrom':'X'},'overlaps_target':False,
          'identity':.8,'query_coverage':.8}]}
        r={'protein_evidence':protein,'retrocopy_evidence':{
          'target_model':{'status':'unavailable'}}}
        self.assertEqual(result(r,'retroposition')['status'],'partial_signal')

    def test_relocation(self):
        r={'synteny_evidence':{'a':{'status':'flanks_collinear',
                                   'locus_status':'homolog_elsewhere'}}}
        self.assertEqual(result(r,'transposition')['status'],'candidate')

    def test_fusion(self):
        a=dict(row(chrom='2L',coverage=.4),query_length=1000,query_start=0,query_end=400)
        b=dict(row(chrom='X',coverage=.4),query_length=1000,query_start=600,query_end=1000)
        self.assertEqual(result({'relatives':{'a':[a,b]}},'fusion')['status'],'candidate')

    def test_unassessed_vs_no_signal(self):
        self.assertEqual(result({},'de_novo')['status'],'not_assessed')
        self.assertEqual(result({},'dna_duplication')['status'],'not_assessed')
        self.assertEqual(result({},'retroposition')['status'],'not_assessed')

    def test_gene_identifiers_have_no_effect(self):
        r={'relatives':{'a':[dict(row(),coding_overlap='overlaps_cds')]}}
        a=classify({**r,'gene_id':'Zeus'})
        b=classify({**r,'gene_id':'CG32582'})
        self.assertEqual(a['options'],b['options'])

class RetrocopyAnnotationTests(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.path=Path(self.tmp.name)/'genes.gff3'
        self.path.write_text(
            '##gff-version 3\n'
            'chr2\tX\tCDS\t101\t200\t.\t+\t0\tID=cds1;Parent=rna-child\n'
            'chrX\tX\tCDS\t1001\t1080\t.\t+\t0\tID=cds2;Parent=rna-parent\n'
            'chrX\tX\tCDS\t1201\t1280\t.\t+\t0\tID=cds3;Parent=rna-parent\n')
    def tearDown(self):self.tmp.cleanup()
    def test_multi_exon_parent(self):
        r=compare_retrocopy(row(chrom='chr2'),[row(chrom='chrX',start=1000,end=1280)],self.path)
        self.assertTrue(r['retrocopy_architecture_signal'])
    def test_no_gff_is_unknown(self):
        r=compare_retrocopy(row(),[row(chrom='X')],None)
        self.assertFalse(r['retrocopy_architecture_signal'])
        self.assertEqual(r['target_model']['status'],'unavailable')
    def test_no_cds_is_not_intronless(self):
        self.assertIsNone(models_for_loci(self.path,[row(chrom='chr3')])[0]['exon_count'])
