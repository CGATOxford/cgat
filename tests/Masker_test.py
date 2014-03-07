##########################################################################
#   Gene prediction pipeline
#
#   $Id: Masker_test.py 2784 2009-09-10 11:41:14Z andreas $
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

import CGAT.Masker as Masker
import unittest


class SegCheck(unittest.TestCase):

    mMasker = Masker.MaskerSeg()

    def testEmpty(self):
        """test empty input."""
        self.assertEqual(self.mMasker(""), "")

    def testProtein(self):
        """test protein input."""
        self.assertEqual(self.mMasker(
            "ACDEFGHIKLWWWWWWWWWWWWWWwwwwwwwwwwwacdefghikl"),
            "ACDEFGHIKLXXXXXXXXXXXXXXxxxxxxxxxxxacdefghikl")

    def testCoding(self):
        """test coding sequence input."""
        self.assertEqual(self.mMasker(
            "GCCTGCGACGAGTTCGGCCACATCAAGCT"
            "GTGGTGGTGGTGGTGGTGGTGGTGGTGGT"
            "GGTGGTGGTGGTGGTGGTGGTGGTGGTGG"
            "tggtggtggtggtggtgggcctgcgacga"
            "gttcggccacatcaagctg"),
            "GCCTGCGACGAGTTCGGCCACATCAAGCT"
            "GNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
            "NNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
            "nnnnnnnnnnnnnnnnnngcctgcgacga"
            "gttcggccacatcaagctg")


class DustMaskerCheck(unittest.TestCase):
    mMasker = Masker.MaskerDustMasker()

if __name__ == "__main__":
    unittest.main()
