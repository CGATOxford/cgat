##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
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
'''
test_tree_strain2species.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python test_tree_strain2species.py --help

Type::

   python test_tree_strain2species.py --help

for command line help.

Documentation
-------------

Code
----

'''
import unittest
import sys
from pprint import pprint

import CGAT.TreeTools as TreeTools
import tree_strain2species


class tree_strain2speciesTest(unittest.TestCase):

    """
    A test class for the feedparser module.
    """

    def setUp(self):
        """
        set up data used in the tests.
        setUp is called before each test function execution.
        """

        class Options:
            separator = "|"
            loglevel = 3
            stdlog = sys.stdout
            stderr = sys.stderr
            stdout = sys.stdout
            pattern_gene = "J%06i"

        self.mTestData = [
            # simple binary tree
            ("((species1a|gene1,species1b|gene2),(species2a|gene1,species2b|gene2));",
             ((('species1a', 'gene1'), ('species1b', 'gene2')),
              (('species2a', 'gene1'), ('species2b', 'gene2'))),
             {'species1a': 'species1',
              'species1b': 'species1',
              'species2a': 'species2',
              'species2b': 'species2'},
             Options()),
            # tree with extra species, trifurcating tree at root
            ("((species1a|gene1,species1b|gene2),(species2a|gene1,species2b|gene2),species3|gene1);",
             ((('species1a', 'gene1'), ('species1b', 'gene2')),
              (('species2a', 'gene1'), ('species2b', 'gene2'))),
             {'species1a': 'species1',
              'species1b': 'species1',
              'species2a': 'species2',
              'species2b': 'species2'},
             Options()),
            # tree with extra species, binary tree
            ("((species1a|gene1,species1b|gene2),((species2a|gene1,species2b|gene2),species3|gene1));",
             ((('species1a', 'gene1'), ('species1b', 'gene2')),
              (('species2a', 'gene1'), ('species2b', 'gene2'))),
             {'species1a': 'species1',
              'species1b': 'species1',
              'species2a': 'species2',
              'species2b': 'species2'},
             Options()),
            # tree with extra species, extra species within synonyms should
            # prevent joining.
            ("((species1a|gene1,species1b|gene2),((species2a|gene1,species3|gene1),species2b|gene2));",
             ((('species1a', 'gene1'), ('species1b', 'gene2')), ),
             {'species1a': 'species1',
              'species1b': 'species1',
              'species2a': 'species2',
              'species2b': 'species2'},
             Options()), ]

    def tearDown(self):
        """
        tear down any data used in tests
        tearDown is called after each test function execution.
        """
        pass

    def testGetMergers(self):
        """
        test.

        TODO: add testing for transcripts
        """
        print "testGetMergers()"

        for lines, reference, map_strain2species, options in self.mTestData:
            nexus = TreeTools.Newick2Nexus(lines)
            mergers = tree_strain2species.getMergers(
                nexus.trees[0], map_strain2species, options)
            for node_id, species, strain_x, gene_x, strain_y, gene_y in mergers:
                key1 = ((strain_x, gene_x), (strain_y, gene_y))
                key2 = ((strain_y, gene_y), (strain_x, gene_x))
                if key1 not in reference and key2 not in reference:
                    self.fail("%s not in reference %s" %
                              (str(key1), str(reference)))

if __name__ == '__main__':
    unittest.main()
