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
    def test_other_mosaic_needs_independent_nonoverlapping_query_segments(self):
        a=h(start=100,end=400);b=h(chrom='3R',start=500000,end=500400)
        a.update(query_coverage=.40,query_start=0,query_end=400,query_length=1000)
        b.update(query_coverage=.40,query_start=600,query_end=1000,query_length=1000)
        r={'relatives':{'a':[a,b]}}
        self.assertGreater(screen(r)['other_mechanism']['support_points'],0)
        b['query_start']=100;b['query_end']=500
        self.assertEqual(screen(r)['other_mechanism']['support_points'],0)
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

class PipelineSummaryIntegrationTests(unittest.TestCase):
    @patch('evogene.ncbi_auto.assess_synteny')
    @patch('evogene.ncbi_auto.focal_flanks')
    @patch('evogene.ncbi_auto.discover_target_locus')
    @patch('evogene.ncbi_auto.annotate')
    @patch('evogene.ncbi_auto.map_genome')
    @patch('evogene.ncbi_auto.download')
    def test_auto_analysis_runs_screens_then_distribution(
            self,download,mapping,annotate,locate,flanks,synteny):
        from pathlib import Path
        from tempfile import TemporaryDirectory
        from evogene.ncbi_auto import analyze
        download.side_effect=lambda name,cache: {
            'species':name,'fasta':name+'.fa','gff':None,
            'accession':'GCF_000000001.1'}
        mapping.side_effect=lambda query,genome,species:[h('overlaps_cds')]
        locate.return_value=(h(),[h(chrom='3R')])
        flanks.return_value=('A'*4000,'G'*4000)
        synteny.return_value={'status':'flanks_collinear',
                              'locus_status':'homolog_overlaps_cds'}
        with TemporaryDirectory() as folder:
            f=Path(folder)/'gene.fasta'
            f.write_text('>gene'+chr(10)+'ATG'*110+chr(10))
            report=analyze(f,'Drosophila simulans','test_gene',count=2)
        self.assertIn('synteny_evidence',report)
        self.assertIn('duplication',report['origin_support']['tests'])
        self.assertGreater(report['origin_support']['tests']['duplication']['support_points'],0)
        self.assertEqual(round(sum(x['support_percent']
                                 for x in report['origin_support']['options']),1),100.0)
