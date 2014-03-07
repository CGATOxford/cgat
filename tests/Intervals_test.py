##########################################################################
#   Gene prediction pipeline
#
#   $Id: Intervals_test.py 2784 2009-09-10 11:41:14Z andreas $
#
#   Copyright (C) 2004 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################
"""unit testing module for the Tree.py class."""

import CGAT.Intervals as Intervals
import unittest


class TruncateCheck(unittest.TestCase):

    def testEmpty(self):
        """test empty input."""
        self.assertEqual(Intervals.truncate([], []), [])

    def testHalfEmpty(self):
        """test empty input."""
        self.assertEqual(Intervals.truncate([], [(0, 5)]), [])
        self.assertEqual(Intervals.truncate([(0, 5)], []), [(0, 5)])

    def testSingle(self):
        """test empty input."""
        self.assertEqual(Intervals.truncate([(0, 5)], [(0, 5)]), [])
        self.assertEqual(Intervals.truncate([(0, 5)], [(0, 3)]), [(3, 5)])
        self.assertEqual(Intervals.truncate([(0, 3)], [(0, 5)]), [])
        self.assertEqual(Intervals.truncate([(0, 5)], [(3, 5)]), [(0, 3)])
        self.assertEqual(Intervals.truncate([(3, 5)], [(0, 5)]), [])
        self.assertEqual(Intervals.truncate([(5, 10)], [(5, 10)]), [])
        self.assertEqual(Intervals.truncate([(5, 10)], [(5, 20)]), [])
        self.assertEqual(Intervals.truncate([(5, 10)], [(0, 10)]), [])
        self.assertEqual(Intervals.truncate([(5, 10)], [(0, 10)]), [])
        self.assertEqual(Intervals.truncate([(5, 10)], [(0, 20)]), [])

    def testMultiple(self):
        """test empty input."""
        self.assertEqual(
            Intervals.truncate([(0, 5), (10, 15)], [(0, 5)]), [(10, 15)])
        self.assertEqual(
            Intervals.truncate([(0, 5), (10, 15)], [(0, 10)]), [(10, 15)])
        self.assertEqual(Intervals.truncate([(0, 5), (10, 15)], [(0, 15)]), [])
        self.assertEqual(Intervals.truncate([(0, 5), (5, 10)], [(0, 10)]), [])
        self.assertEqual(
            Intervals.truncate([(0, 5), (5, 10)], []), [(0, 5), (5, 10)])

    def testNoOverlap(self):
        """test empty input."""
        self.assertEqual(
            Intervals.truncate([(0, 5), (10, 15)], [(5, 10)]), [(0, 5), (10, 15)])
        self.assertEqual(
            Intervals.truncate([(5, 10)], [(0, 5), (10, 15)]), [(5, 10)])
        self.assertEqual(
            Intervals.truncate([(0, 5), (5, 10)], [(10, 15)]), [(0, 5), (5, 10)])


class IntersectCheck(unittest.TestCase):

    def testEmpty(self):
        """test empty input."""
        self.assertEqual(Intervals.intersect([], []), [])

    def testHalfEmpty(self):
        """test empty input."""
        self.assertEqual(Intervals.intersect([(0, 5)], []), [])
        self.assertEqual(Intervals.intersect([], [(0, 5)]), [])

    def testSingle(self):
        """test empty input."""
        self.assertEqual(Intervals.intersect([(0, 5)], [(0, 5)]), [(0, 5)])
        self.assertEqual(Intervals.intersect([(0, 5)], [(0, 3)]), [(0, 3)])
        self.assertEqual(Intervals.intersect([(0, 3)], [(0, 5)]), [(0, 3)])
        self.assertEqual(Intervals.intersect([(0, 5)], [(3, 5)]), [(3, 5)])
        self.assertEqual(Intervals.intersect([(3, 5)], [(0, 5)]), [(3, 5)])
        self.assertEqual(Intervals.intersect([(5, 10)], [(5, 20)]), [(5, 10)])
        self.assertEqual(Intervals.intersect([(5, 10)], [(0, 20)]), [(5, 10)])

    def testMultiple(self):
        """test empty input."""
        self.assertEqual(
            Intervals.intersect([(0, 5), (10, 15)], [(0, 5)]), [(0, 5)])
        self.assertEqual(
            Intervals.intersect([(0, 5), (10, 15)], [(0, 10)]), [(0, 5)])
        self.assertEqual(
            Intervals.intersect([(0, 5), (10, 15)], [(0, 15)]), [(0, 5), (10, 15)])
        self.assertEqual(
            Intervals.intersect([(0, 5), (5, 10)], [(0, 10)]), [(0, 5), (5, 10)])

    def testNoOverlap(self):
        """test empty input."""
        self.assertEqual(
            Intervals.intersect([(0, 5), (10, 15)], [(5, 10)]), [])
        self.assertEqual(
            Intervals.intersect([(5, 10)], [(0, 5), (10, 15)]), [])


class FromArrayCheck(unittest.TestCase):

    def testEmpty(self):
        """test empty input."""
        self.assertEqual(Intervals.fromArray([]), [])

    def testArray1(self):
        """test simple array."""
        a = [1, 1, 1, 0, 0, 0, 1, 1, 1]
        self.assertEqual(Intervals.fromArray(a), [(0, 3), (6, 9)])
        self.assertEqual(Intervals.fromArray([not x for x in a]), [(3, 6)])

    def testArray2(self):
        """test longer array."""
        a = [1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1]
        self.assertEqual(Intervals.fromArray(a), [(0, 3), (6, 9), (12, 15)])
        self.assertEqual(
            Intervals.fromArray([not x for x in a]), [(3, 6), (9, 12)])


if __name__ == "__main__":
    unittest.main()
