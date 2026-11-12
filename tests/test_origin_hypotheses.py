import unittest
from unittest.mock import patch
from evogene.origin_tests import screen,distribute
from evogene.genome_probes import discover_target_locus,assess_synteny
def h(cds='unknown_no_gff',chrom='2L',start=100,end=400):
    return {'query_coverage':.94,'identity_approx':.9,'coding_overlap':cds,
            'chrom':chrom,'start':start,'end':end,'strand':'+'}
class OriginHypothesisTests(unittest.TestCase):
    def test_denovo_needs_two_syntenic_noncoding_homologs(self):
        r={'relatives':{'a':[h('no_cds_overlap_not_proof_of_noncoding')],
                        'b':[h('no_cds_overlap_not_proof_of_noncoding')]}}
        self.assertEqual(screen(r)['de_novo']['support_points'],0)
        r['synteny_evidence']={s:{'status':'flanks_collinear','locus_status':'homolog_no_cds_overlap'} for s in ('a','b')}
        self.assertGreater(screen(r)['de_novo']['support_points'],0)
        r['relatives']['b']=[h('overlaps_cds')]
        self.assertEqual(screen(r)['de_novo']['support_points'],0)
    def test_duplication_needs_extra_locus(self):
        r={'relatives':{'a':[h()]}}
        self.assertEqual(screen(r)['duplication']['support_points'],0)
        r['extra_target_copies']=[h(chrom='3R')]
        self.assertGreater(screen(r)['duplication']['support_points'],0)
    def test_transposition_needs_displacement_and_synteny(self):
        r={'relatives':{'a':[h()]}}
        self.assertEqual(screen(r)['transposition']['support_points'],0)
        r['synteny_evidence']={'a':{'status':'flanks_collinear','locus_status':'homolog_elsewhere'}}
        self.assertGreater(screen(r)['transposition']['support_points'],0)
    def test_conservation_with_cds(self):
        r={'relatives':{'a':[h('overlaps_cds')]}}
        self.assertGreater(screen(r)['ancestral_conservation']['support_points'],0)
    def test_other_is_unassessed(self):
        self.assertEqual(screen({'relatives':{}})['other_mechanism']['status'],'not_assessed')
    def test_100_percent_rounding_and_uncertainty(self):
        for n in range(6):
            r={'relatives':{str(i):[h('overlaps_cds')] for i in range(n)}}
            d=distribute(r)
            self.assertEqual(round(sum(x['support_percent'] for x in d['options']),1),100)
            self.assertFalse(d['is_calibrated_probability'])
    @patch('evogene.genome_probes.map_genome')
    def test_target_multiple_copies(self,m):
        m.return_value=[h(),h(start=200,end=450),h(chrom='3R',start=5000,end=5300)]
        main,extra=discover_target_locus('A'*300,'file')
        self.assertEqual(main['chrom'],'2L')
        self.assertEqual(len(extra),1)
    @patch('evogene.genome_probes.map_genome')
    def test_anchor_synteny(self,m):
        m.side_effect=[[h(start=100,end=900)],[h(start=1700,end=2500)]]
        r=assess_synteny(('A'*4000,'G'*4000),
             [h('overlaps_cds',start=1200,end=1400)],'genome','x')
        self.assertEqual(r['locus_status'],'homolog_overlaps_cds')
