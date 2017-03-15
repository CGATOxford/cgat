import unittest
import pandas
import CGAT.Counts as Counts


class TestCountsClass(unittest.TestCase):

    def setUp(self):

        self.counts = Counts.Counts(pandas.DataFrame(
            {'sample1': [0, 1, 2],
             'sample2': [2, 4, 3]}))

    def test_removeSamples_removes_low_counts(self):
        self.counts.removeSamples(min_counts_per_sample=3)

        self.assertTrue(
            self.counts.table.equals(
                pandas.DataFrame(
                    {'sample2': [2, 4, 3]})))

    def test_removeSamples_fails_with_string(self):
        self.assertRaises(
            self.counts.removeSamples,
            min_counts_per_sample='3')
