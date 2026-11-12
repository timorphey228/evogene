import unittest
from pathlib import Path
from evogene.core import Evidence, classify
from evogene.io import read_properties, read_evidence
from evogene.orf import find_orfs

class CoreTests(unittest.TestCase):
    def test_denovo_requires_synteny_and_outgroups(self):
        e = Evidence('gene1', paralog_detected=False, distant_coding_homolog=False,
                     syntenic_outgroup_locus=True, outgroup_locus_noncoding=True,
                     homolog_search_completed=True, outgroup_count=3)
        self.assertEqual(classify(e).hypothesis, 'de_novo')
        e.outgroup_count = 1
        self.assertEqual(classify(e).hypothesis, 'undetermined')

    def test_paralog_prevents_denovo(self):
        e = Evidence('gene2', paralog_detected=True, syntenic_outgroup_locus=True,
                     outgroup_locus_noncoding=True, homolog_search_completed=True, outgroup_count=3)
        self.assertEqual(classify(e).hypothesis, 'duplication')
        self.assertTrue(classify(e).conflicting)

    def test_unknown_is_not_negative(self):
        self.assertEqual(classify(Evidence('gene')).hypothesis, 'undetermined')

    def test_ambiguous_positive_evidence(self):
        e = Evidence('gene', paralog_detected=True, distant_coding_homolog=True)
        self.assertEqual(classify(e).hypothesis, 'undetermined')

    def test_orf_strands(self):
        dna = 'ATG' + 'AAA' * 30 + 'TAA'
        self.assertEqual(find_orfs(dna)[0].aa_length, 31)
        self.assertIn('-', [x.strand for x in find_orfs(dna.translate(str.maketrans('ACGT','TGCA'))[::-1])])

    def test_demo_csv(self):
        rows = read_evidence('data/demo_evidence.csv')
        self.assertEqual(len(rows), 3)
        self.assertEqual([classify(x).hypothesis for x in rows], ['de_novo', 'duplication', 'undetermined'])

class SpreadsheetTests(unittest.TestCase):
    def test_uploaded_workbook(self):
        p = Path('/mnt/data/41467_2024_45028_MOESM10_ESM.xlsx')
        if p.exists():
            rows = read_properties(p)
            self.assertIn('FBgn0040784', rows)
            self.assertEqual(rows['FBgn0040784']['De novo or not?'], 'de novo')
            self.assertEqual(rows['FBgn0040784']['de novo gene age'], 'Br1')

if __name__ == '__main__':
    unittest.main()
