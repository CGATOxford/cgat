import unittest, sys
from pprint import pprint

import TreeTools
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
            ## simple binary tree
            ("((species1a|gene1,species1b|gene2),(species2a|gene1,species2b|gene2));",
             ( ( ('species1a', 'gene1'), ('species1b', 'gene2')),
               ( ('species2a', 'gene1'), ('species2b', 'gene2') ) ), 
            { 'species1a' : 'species1',
              'species1b' : 'species1',
              'species2a' : 'species2',
              'species2b' : 'species2' },
            Options()),
            ## tree with extra species, trifurcating tree at root
            ("((species1a|gene1,species1b|gene2),(species2a|gene1,species2b|gene2),species3|gene1);",
             ( ( ('species1a', 'gene1'), ('species1b', 'gene2')),
               ( ('species2a', 'gene1'), ('species2b', 'gene2') ) ), 
            { 'species1a' : 'species1',
              'species1b' : 'species1',
              'species2a' : 'species2',
              'species2b' : 'species2' },
            Options()),
            ## tree with extra species, binary tree
            ("((species1a|gene1,species1b|gene2),((species2a|gene1,species2b|gene2),species3|gene1));",
             ( ( ('species1a', 'gene1'), ('species1b', 'gene2')),
               ( ('species2a', 'gene1'), ('species2b', 'gene2') ) ), 
            { 'species1a' : 'species1',
              'species1b' : 'species1',
              'species2a' : 'species2',
              'species2b' : 'species2' },
            Options()),             
            ## tree with extra species, extra species within synonyms should prevent joining.
            ("((species1a|gene1,species1b|gene2),((species2a|gene1,species3|gene1),species2b|gene2));",
             ( ( ('species1a', 'gene1'), ('species1b', 'gene2') ), ),
            { 'species1a' : 'species1',
              'species1b' : 'species1',
              'species2a' : 'species2',
              'species2b' : 'species2' },
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
            nexus = TreeTools.Newick2Nexus( lines )
            mergers = tree_strain2species.getMergers( nexus.trees[0], map_strain2species, options )
            for node_id, species, strain_x, gene_x, strain_y, gene_y in mergers:
                key1 = ((strain_x, gene_x), (strain_y, gene_y))
                key2 = ((strain_y, gene_y), (strain_x, gene_x))                
                if key1 not in reference and key2 not in reference:
                    self.fail( "%s not in reference %s" % (str(key1), str(reference)))
                    
if __name__ == '__main__':
    unittest.main()


