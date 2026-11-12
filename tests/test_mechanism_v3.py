import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from evogene.mechanism_v3 import classify
from evogene.gene_models import compare_retrocopy,models_for_loci

def row(chrom='2L', start=100, end=200, coverage=.9, identity=.9):
    return {'chrom':chrom,'start':start,'end':end,
            'query_coverage':coverage,'identity_approx':identity}

class MechanismSeparationTests(unittest.TestCase):
    def test_conserved_dna_does_not_count_as_mechanism(self):
        report={'relatives':{'simulans':[dict(row(),coding_overlap='overlaps_cds')],
                             'sechellia':[dict(row(),coding_overlap='overlaps_cds')]}}
        m=classify(report)
        self.assertEqual(m['leading_candidate'],None)
        self.assertEqual(m['options'][-1]['support_percent'],100.0)
        self.assertEqual(len(m['conservation']['cds_overlap_species']),2)

    def test_retrocopy_beats_non_specific_duplication(self):
        report={'extra_target_copies':[row(chrom='X')],
                'focal_reference':{'accession':'test'},
                'retrocopy_evidence':{'target_model':{'status':'annotated','exon_count':1},
                    'parent_candidates':[{'parent_multi_exon_and_target_single_exon':True}]}}
        m=classify(report)
        self.assertEqual(m['leading_candidate'],'retroposition')
        support={x['mechanism']:x['evidence_weight'] for x in m['tests']}
        self.assertGreater(support['retroposition'],support['dna_duplication'])

    def test_duplication_signal(self):
        m=classify({'extra_target_copies':[row(chrom='3R')],'focal_reference':{'accession':'a'}})
        self.assertEqual(m['leading_candidate'],'dna_duplication')

    def test_de_novo_needs_syntenic_noncoding_outgroups(self):
        report={'focal_reference':{'accession':'x'},'relatives':{'a':[],'b':[]}}
        self.assertIsNone(classify(report)['leading_candidate'])
        report['synteny_evidence']={
            'a':{'status':'flanks_collinear','locus_status':'homolog_no_cds_overlap'},
            'b':{'status':'flanks_collinear','locus_status':'homolog_no_cds_overlap'}}
        self.assertEqual(classify(report)['leading_candidate'],'de_novo')

    def test_transposition_signal(self):
        m=classify({'focal_reference':{'accession':'x'},
            'synteny_evidence':{'a':{'status':'flanks_collinear','locus_status':'homolog_elsewhere'}}})
        self.assertEqual(m['leading_candidate'],'transposition')

    def test_fusion_split_homology(self):
        a=dict(row(chrom='2L',coverage=.40),query_length=1000,query_start=0,query_end=400)
        b=dict(row(chrom='X',coverage=.40),query_length=1000,query_start=600,query_end=1000)
        m=classify({'focal_reference':{'accession':'x'},'relatives':{'a':[a,b]}})
        self.assertEqual(m['leading_candidate'],'fusion')

    def test_distribution_sums_exactly_100(self):
        for n in range(6):
            m=classify({'relatives':{str(i):[row()] for i in range(n)}})
            self.assertEqual(round(sum(x['support_percent'] for x in m['options']),1),100.0)
            self.assertFalse(m['is_calibrated_probability'])

class RetrocopyAnnotationTests(unittest.TestCase):
    def setUp(self):
        self.folder=TemporaryDirectory()
        self.gff=Path(self.folder.name)/'genes.gff3'
        self.gff.write_text(
            '##gff-version 3\n'
            'chr2\tX\tCDS\t101\t200\t.\t+\t0\tID=cds1;Parent=rna-child\n'
            'chrX\tX\tCDS\t1001\t1080\t.\t+\t0\tID=cds2;Parent=rna-parent\n'
            'chrX\tX\tCDS\t1201\t1280\t.\t+\t0\tID=cds3;Parent=rna-parent\n')
    def tearDown(self):
        self.folder.cleanup()
    def test_gene_model_detects_spliced_parent(self):
        target=row(chrom='chr2',start=100,end=200)
        parent=row(chrom='chrX',start=1000,end=1280,coverage=.7)
        r=compare_retrocopy(target,[parent],self.gff)
        self.assertEqual(r['target_model']['exon_count'],1)
        self.assertEqual(r['parent_candidates'][0]['candidate_parent_model']['exon_count'],2)
        self.assertTrue(r['retrocopy_architecture_signal'])
    def test_missing_gff_is_unknown_not_intronless(self):
        a=row(chrom='chr2');b=row(chrom='chrX')
        r=compare_retrocopy(a,[b],None)
        self.assertFalse(r['retrocopy_architecture_signal'])
        self.assertEqual(r['target_model']['status'],'unavailable')
    def test_no_overlapping_cds_not_intronless(self):
        r=models_for_loci(self.gff,[row(chrom='chr3')])
        self.assertEqual(r[0]['status'],'no_annotated_cds')
        self.assertIsNone(r[0]['exon_count'])
