import unittest
import os

import CGAT.GTF as GTF
import CGAT.IOTools as IOTools


class TestIteration(unittest.TestCase):

    filename = os.path.join("data", "hg19.small.gtf.gz")

    def test_number_of_intervals_is_correct(self):

        with IOTools.openFile(self.filename) as inf:
            records = list(GTF.iterator(inf))

        self.assertEqual(len(records),
                         100)


if __name__ == "__main__":
    unittest.main()
