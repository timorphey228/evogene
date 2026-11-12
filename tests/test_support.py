import unittest
from evogene.support import distribution
class SupportTests(unittest.TestCase):
    def test_no_evidence(self):
        d=distribution({"relatives":{}})
        self.assertEqual(d["options"][-1]["support_percent"],100)
        self.assertFalse(d["is_calibrated_probability"])
    def test_coding_support(self):
        h={"query_coverage":.95,"identity_approx":.92,"coding_overlap":"overlaps_cds"}
        d=distribution({"relatives":{"dsim":[h],"dyak":[h]}})
        x={z["mechanism"]:z["support_percent"] for z in d["options"]}
        self.assertGreater(x["ancestral_conservation"],0)
        self.assertEqual(x["de_novo"],0)
        self.assertEqual(x["transposition"],0)
        self.assertEqual(sum(x.values()),100)
    def test_failures_increase_uncertainty(self):
        h={"query_coverage":.95,"identity_approx":.92}
        a=distribution({"relatives":{"a":[h]}})
        b=distribution({"relatives":{"a":[h]},"download_or_analysis_errors":[{"species":"b"}]})
        self.assertGreater(b["options"][-1]["support_percent"],a["options"][-1]["support_percent"])
    def test_duplication_only_with_extra(self):
        h={"query_coverage":.8,"identity_approx":.8}
        d=distribution({"relatives":{"a":[h]},"extra_target_copies":[h]})
        self.assertGreater(d["options"][1]["support_percent"],0)
    def test_rounding_sums_100(self):
        for n in range(6):
            h={"query_coverage":.8,"identity_approx":.8,"coding_overlap":"overlaps_cds"}
            d=distribution({"relatives":{str(i):[h] for i in range(n)}})
            self.assertEqual(round(sum(x["support_percent"] for x in d["options"]),1),100)
if __name__=="__main__":unittest.main()
