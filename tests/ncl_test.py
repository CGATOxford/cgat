"""test high-level interface."""

import unittest
import random
import tempfile
import shutil
import os
from CGAT.NCL import *


class TestNCLSimpleNegativeIntervals(unittest.TestCase):

    def setUp(self):
        self.mIndex = NCLSimple

    def testBuild1(self):
        index = self.mIndex()
        self.assertRaises(ValueError, index.add, 0, 0)

    def testBuild2(self):
        index = self.mIndex()
        self.assertRaises(ValueError, index.add, -10, 0)


class TestNCLSimple(unittest.TestCase):

    def setUp(self):
        self.l = [(10, 20, 0),
                  (15, 25, 1),
                  (30, 50, 2),
                  ]

        self.tests = (((10, 15), (0,)),
                      ((10, 30), (0, 1)),
                      ((0, 10), (), ),
                      ((10, 50), (0, 1, 2), ),
                      ((25, 30), ()),
                      )

        self.mIndex = NCLSimple

    def checkIntervals(self, index):

        for a, b in self.tests:
            result = tuple(sorted([x[2] for x in index.find(a[0], a[1])]))
            self.assertEqual(result, b)

    def buildIndex(self, l):
        index = self.mIndex()
        for start, end, value in self.l:
            index.add(start, end)
        return index

    def testBuild(self):
        index = self.buildIndex(self.l)
        self.assertRaises(ValueError, index.add, 0, 0)

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
            r = len(list(index.find(x, x + 1)))
            self.assertEqual(
                r > 0, bits[x], "invalid return at position %i (expected %i, got %i)" % (x, bits[x], r))

    def testEmptyIntervals(self):
        bits = [0] * 100

        index = self.buildIndex(self.l)

        for x in xrange(0, len(bits)):
            self.assertRaises(IndexError, index.find, x, x)


class TestNCL(TestNCLSimple):

    def setUp(self):
        TestNCLSimple.setUp(self)
        self.mIndex = NCL

    def buildIndex(self, l):
        index = self.mIndex()
        for start, end, value in self.l:
            index.add(start, end, value)
        return index

    def testBuild(self):
        index = self.buildIndex(self.l)
        self.assertRaises(ValueError, index.add, 0, 0, 4)


class TestNCLDisk(TestNCLSimple):

    """basic tests on database."""

    def setUp(self):
        TestNCLSimple.setUp(self)
        self.mIndex = NCL
        self.tmpdir = tempfile.mkdtemp()
        self.tmpfile = os.path.join(self.tmpdir, "tmp")

    def buildIndex(self, l):
        index = self.mIndex(filestem=self.tmpfile)
        for start, end, value in self.l:
            index.add(start, end, value)
        del index
        index = self.mIndex(filestem=self.tmpfile)
        return index

    def testBuild(self):
        index = self.buildIndex(self.l)
        self.assertRaises(ValueError, index.add, 0, 0, 4)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

if __name__ == '__main__':
    unittest.main()
