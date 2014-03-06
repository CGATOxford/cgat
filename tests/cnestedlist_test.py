"""test low-level interface to cnestedlist."""

import unittest
import random
import tempfile
import shutil
import os
from CGAT.NCL import *


class TestIntervalDB(unittest.TestCase):

    def setUp(self):
        self.l = [(10, 20, 1),
                  (15, 25, 2),
                  (30, 50, 3),
                  ]

        self.tests = (((10, 15), (1,)),
                      ((10, 30), (1, 2)),
                      ((0, 10), (), ),
                      ((10, 50), (1, 2, 3), ),
                      ((25, 30), ()),
                      )

    def checkIntervals(self, index):

        for a, b in self.tests:
            result = tuple(
                sorted([x[2] for x in index.find_overlap(a[0], a[1])]))
            self.assertEqual(result, b)

    def testBuild(self):
        index = IntervalDB()
        index.fromlist(self.l)

        # test adding an empty interval
        index = IntervalDB()
        l = self.l[:]
        index.fromlist(l)
        l.append((0, 0, 4))
        # AH: disabled, need to check if raising
        # exception is correct.
        # self.assertRaises(ValueError, index.fromlist, l)

    def buildIndex(self, l):
        index = IntervalDB()
        index.fromlist(self.l)
        return index

    def testSortedList(self):
        """test sorted list."""
        l = self.l[:]
        l.sort()
        self.checkIntervals(self.buildIndex(l))

    def testReverseList(self):
        """test reverse sorted list."""
        l = self.l[:]
        l.sort()
        l.reverse()
        self.checkIntervals(self.buildIndex(l))

    def testRandomizedList(self):
        """test randomized list."""
        l = self.l[:]
        random.shuffle(l)
        self.checkIntervals(self.buildIndex(l))

    def testPointIntervals(self):
        bits = [0] * 100
        for start, end, val in self.l:
            bits[start:end] = [1] * (end - start)

        index = self.buildIndex(self.l)

        for x in xrange(0, len(bits)):
            r = len(list(index.find_overlap(x, x + 1)))
            self.assertEqual(
                r > 0, bits[x], "invalid return at position %i (expected %i, got %i)" % (x, bits[x], r))

    def testEmptyIntervals(self):
        bits = [0] * 100

        index = self.buildIndex(self.l)

        for x in xrange(0, len(bits)):
            self.assertRaises(IndexError, index.find_overlap, x, x)


class TestIntervalFileDB(TestIntervalDB):

    def setUp(self):
        TestIntervalDB.setUp(self)
        self.tmpdir = tempfile.mkdtemp()
        self.tmpfile = os.path.join(self.tmpdir, "tmp")

    def buildIndex(self, l):
        tmp = IntervalDB()
        tmp.fromlist(self.l)
        tmp.write_binaries(self.tmpfile)
        index = IntervalFileDB(self.tmpfile)
        return index

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

if __name__ == '__main__':
    unittest.main()
