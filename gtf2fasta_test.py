################################################################################
#   Gene prediction pipeline 
#
#   $Id: gtf2fasta_test.py 2822 2009-11-24 16:08:44Z andreas $
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
#################################################################################
"""unit testing module for the Tree.py class."""

import sys, os, shutil, optparse, random
import unittest, tempfile
import gtf2fasta
import IndexedFasta

import GTF, GFF

class GeneralTest(unittest.TestCase):

    mDummy = True
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

        self.outfile_genome = os.path.join( self.tmpdir, "genome_in" )
        self.outfile_gtf = os.path.join( self.tmpdir, "exons.gtf" )
        self.outfile_output = os.path.join( self.tmpdir, "output" )

        self.length = 1000
        
        genome = iter( ( ("chr1", "A" * self.length), )  )

        IndexedFasta.createDatabase( self.outfile_genome, genome )
        self.reference = ["g"] * self.length

    def tearDown( self ):
        shutil.rmtree( self.tmpdir )

    def testRun( self ):
        
        if self.mDummy: return

        gtf2fasta.main( [
                "dummy",
                "--genome-file=%s" % self.outfile_genome,
                "--stdin=%s" %  self.outfile_gtf,
                "--output-filename-pattern=%s/%%s" % self.tmpdir,
                "--stdout=%s" % self.outfile_output ] )


        result = open(self.outfile_output).readlines()[1][:-1]
        reference = "".join(self.reference)

        self.assertEqual( len(result), 
                          len(reference), 
                          "unequal length %i != %i" % (len(result), len(reference) ) )
        for x in xrange( len(result)):
            self.assertEqual( result[x], reference[x], "mismatch in characters at postition %i: %s!=%s" %\
                                  (x, result[x], reference[x] ) ) 

    def addIntron( self, start, end, is_positive ):

        if is_positive: i,s = "i", "s"
        else: i,s = "I", "S"

        self.reference[start:end] = [i] * 50
        self.reference[start:start+2] = [s] * 2
        self.reference[end-2:end] = [s] * 2

class TestForwardStrand( GeneralTest ):
    mDummy = False
    def setUp( self ):

        GeneralTest.setUp(self)

        outfile = open( self.outfile_gtf, "w" )

        # one gene of < 200 codons
        gene = "abc" * 200

        entry = GTF.Entry()
        entry.contig = "chr1"
        entry.strand = "+"
        entry.feature = "exon"
        entry.gene_id = entry.transcript_id = "gene1"
        entry.start = 50
        entry.end = 100
        entry.frame = 0
        entry.source = "protein_coding"

        entry.feature = "exon"
        outfile.write( str(entry) + "\n" )
        entry.feature = "CDS"
        outfile.write( str(entry) + "\n" )

        self.reference[entry.start:entry.end] = gene[0:50]
        self.addIntron(entry.end,entry.start+100, True) 

        entry.start += 100
        entry.end += 100
        entry.frame = 1

        entry.feature = "exon"
        outfile.write( str(entry) + "\n" )
        entry.feature = "CDS"
        outfile.write( str(entry) + "\n" )

        self.reference[entry.start:entry.end] = gene[50:100]
        self.addIntron(entry.end,entry.start+100, True) 

        entry.start += 100
        entry.end += 100
        entry.frame = 2
        entry.feature = "exon"
        outfile.write( str(entry) + "\n" )
        entry.feature = "CDS"
        outfile.write( str(entry) + "\n" )

        self.reference[entry.start:entry.end] = gene[100:150]
        self.addIntron(entry.end,entry.start+100, True) 

        entry.start += 100
        entry.end += 100
        entry.frame = 0
        entry.feature = "exon"
        outfile.write( str(entry) + "\n" )
        entry.feature = "CDS"
        outfile.write( str(entry) + "\n" )

        self.reference[entry.start:entry.end] = gene[150:200]

        outfile.close()

class TestReverseStrand( GeneralTest ):
    """still to do"""
    mDummy = False
    def setUp( self ):

        GeneralTest.setUp(self)

        outfile = open( self.outfile_gtf, "w" )

        # one gene of < 200 codons
        gene = "abc" * 200

        entry = GTF.Entry()
        entry.contig = "chr1"
        entry.strand = "+"
        entry.feature = "exon"
        entry.gene_id = entry.transcript_id = "gene1"
        entry.start = 50
        entry.end = 100
        entry.frame = 0
        entry.source = "protein_coding"

        entry.feature = "exon"
        outfile.write( str(entry) + "\n" )
        entry.feature = "CDS"
        outfile.write( str(entry) + "\n" )

        self.reference[entry.start:entry.end] = gene[0:50]
        self.addIntron(entry.end,entry.start+100, True) 

        entry.start += 100
        entry.end += 100
        entry.frame = 1

        entry.feature = "exon"
        outfile.write( str(entry) + "\n" )
        entry.feature = "CDS"
        outfile.write( str(entry) + "\n" )

        self.reference[entry.start:entry.end] = gene[50:100]
        self.addIntron(entry.end,entry.start+100, True) 

        entry.start += 100
        entry.end += 100
        entry.frame = 2
        entry.feature = "exon"
        outfile.write( str(entry) + "\n" )
        entry.feature = "CDS"
        outfile.write( str(entry) + "\n" )

        self.reference[entry.start:entry.end] = gene[100:150]
        self.addIntron(entry.end,entry.start+100, True) 

        entry.start += 100
        entry.end += 100
        entry.frame = 0
        entry.feature = "exon"
        outfile.write( str(entry) + "\n" )
        entry.feature = "CDS"
        outfile.write( str(entry) + "\n" )

        self.reference[entry.start:entry.end] = gene[150:200]

        outfile.close()

if __name__ == "__main__":
    unittest.main()
