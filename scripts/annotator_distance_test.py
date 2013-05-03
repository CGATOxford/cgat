################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: snp2table.py 2861 2010-02-23 17:36:32Z andreas $
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
#################################################################################
'''
annotator_distance_test.py - test script for annotator_distance.py
==================================================================

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

   python <script_name>.py --help

Type::

   python <script_name>.py --help

for command line help.

Documentation
-------------

Code
----

'''


import sys
import os
import shutil
import optparse
import random
import unittest
import tempfile
import annotator_distance

import CGAT.GTF as GTF
import CGAT.GFF as GFF

class AnnotatorDistanceCheck(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.workspace = os.path.join( self.tmpdir, "workspace.gtf" )
        self.segments = os.path.join( self.tmpdir, "segments.gtf" )

        
    def tearDown( self ):
        shutil.rmtree( self.tmpdir )
        
class AnnotatorCheckWorkspaceIntergenic( AnnotatorDistanceCheck ):

    def setUp( self ):
 
        AnnotatorDistanceCheck.setUp(self)
        
        outfile = open( self.workspace, "w" )
        e = GTF.Entry()

        e.contig, e.strand, e.gene_id, e.transcript_id = "chr1", "+", "gene1", "trans1"

        e.start, e.end = 0, 1000
        outfile.write( str(e) + "\n" )

        e.start, e.end = 3000, 4000
        outfile.write( str(e) + "\n" )

        e.contig, e.strand, e.gene_id, e.transcript_id = "chr1", "+", "gene2", "trans1"

        e.start, e.end = 10000, 11000
        outfile.write( str(e) + "\n" )

        e.start, e.end = 13000, 14000
        outfile.write( str(e) + "\n" )

        e.contig, e.strand, e.gene_id, e.transcript_id = "chr1", "-", "gene3", "trans1"

        e.start, e.end = 20000, 21000
        outfile.write( str(e) + "\n" )

        e.start, e.end = 23000, 24000
        outfile.write( str(e) + "\n" )

        outfile.close()

    def testWorkspaceGFF( self ):
        
        workspace = annotator_distance.readWorkspace( open( self.workspace ), 
                                                      workspace_builder = "gff",
                                                      label = "none" )

        self.assertEqual( workspace, {'chr1': [(0, 1000), (3000, 4000), (10000, 11000), (13000, 14000), (20000, 21000), (23000, 24000)]} )

    def testWorkspaceGTFIntergenic1( self ):
        
        workspace = annotator_distance.readWorkspace( open( self.workspace ), 
                                                      workspace_builder = "gtf-intergenic",
                                                      label = "none" )

        self.assertEqual( workspace, {'chr1': [(4000, 10000, ( ('X', ), ('X', ))), (14000, 20000, ( ('X',) , ('X',)))]} )

    def testWorkspaceGTFIntergenic2( self ):
        
        workspace = annotator_distance.readWorkspace( open( self.workspace ), 
                                                      workspace_builder = "gtf-intergenic",
                                                      label = "direction" )
        self.assertEqual( workspace, {'chr1': [(4000, 10000, (('3',), ('5',))), (14000, 20000, (('3',), ('3',)))]} )
 

class AnnotatorCheckWorkspaceIntronic( AnnotatorDistanceCheck ):

    def setUp( self ):
 
        AnnotatorDistanceCheck.setUp(self)

        e = GTF.Entry()
        e.contig, e.strand = "chr1", "+"
        
        outfile = open( self.workspace, "w" )
        start, inc, size = 0, 1000, 100
        for x in range( 0, 2):
            start = x * 10000

            for y in range( 0, 3 ):
                e.gene_id, e.transcript_id = "gene_%i" % x , "transcript_%i" % y
                e.start, e.end = start, start + size
                outfile.write( str(e) + "\n" )
                start += inc

            if e.strand == "+": e.strand = "-"
            else: e.strand = "+"
        outfile.close()

    def testWorkspaceGTFIntronic1( self ):
        
        workspace = annotator_distance.readWorkspace( open( self.workspace ), 
                                                      workspace_builder = "gtf-intronic",
                                                      label = "none" )


        self.assertEqual( workspace, {'chr1': [(100, 1000, (('X',), ('X',))), (1100, 2000, (('X',), ('X',))), (10100, 11000, (('X',), ('X',))), (11100, 12000, (('X',), ('X',)))]} )

    def testWorkspaceGTFIntronic2( self ):
        
        workspace = annotator_distance.readWorkspace( open( self.workspace ), 
                                                      workspace_builder = "gtf-intronic",
                                                      label = "direction" )
        self.assertEqual( workspace, {'chr1': [(100, 1000, (('3',), ('5',))), (1100, 2000, (('3',), ('5',))), (10100, 11000, (('5',), ('3',))), (11100, 12000, (('5',), ('3',)))]} )

 

class AnnotatorDistanceCheckUniform(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.workspace = os.path.join( self.tmpdir, "workspace.gtf" )
        self.segments = os.path.join( self.tmpdir, "segments.gtf" )

        outfile = open( self.workspace, "w" )
        e = GTF.Entry()

        e.contig, e.strand, e.gene_id, e.transcript_id = "chr1", "+", "gene1", "trans1"

        # 10kb genes every 100000kb for 10Mb
        for x in range( 100000, 10000000, 100000 ):
            e.gene_id, e.transcript_id = "gene%i" % x, "trans1"
            e.start, e.end = x, x + 10000
            outfile.write( str(e) + "\n" )
        outfile.close()

        # segments: uniformly distributed every 1kb with random length
        outfile = open( self.segments, "w" )
        e = GTF.Entry()
        e.contig, e.strand  = "chr1", "+"
        for x in range( 0, 10000000, 1000 ):
            e.gene_id, e.transcript_id = "gene%i" % x, "trans1"
            e.start, e.end = x, x+random.randint( 50, 150 )
            outfile.write( str(e) + "\n" )
            
        outfile.close()

    def testMain( self ):
        pass
        #annotator_distance.main( ["--workspace=%s" % self.workspace, "--segments=%s" %  self.segments,
        #                          "--workspace-labels=direction", 
        #                          "--workspace-builder=gtf-intergenic" ])

        
class AnnotatorDistanceCheck3End(AnnotatorDistanceCheck):
    """equally sized intergenic space, over-representation of 3' end."""
    
    def setUp(self):
        AnnotatorDistanceCheck.setUp(self)

        outfile = open( self.workspace, "w" )
        e = GTF.Entry()

        e.contig, e.strand, e.gene_id, e.transcript_id = "chr1", "+", "gene1", "trans1"

        # 10kb genes every 100000kb for 10Mb
        for x in range( 100000, 10000000, 100000 ):
            e.gene_id, e.transcript_id = "gene%i" % x, "trans1"
            e.start, e.end = x, x + 10000
            outfile.write( str(e) + "\n" )
        outfile.close()

        # segments: concentrated at 5' end
        outfile = open( self.segments, "w" )
        e = GTF.Entry()
        e.contig, e.strand  = "chr1", "+"
        for x in range( 110000, 10000000, 100000 ):
            y = x
            inc = 200
            while y < x + 100000:
                e.gene_id, e.transcript_id = "gene%i" % (y), "trans1"
                e.start, e.end = y, y+random.randint( 50, 150 )
                outfile.write( str(e) + "\n" )
                y += inc
                inc += random.randint( 0, 100 )

        outfile.close()

    def testMain( self ):
        return
        annotator_distance.main( ["--workspace=%s" % self.workspace, "--segments=%s" %  self.segments,
                                  "--workspace-labels=direction", 
                                  "--workspace-builder=gtf-intergenic" ])

class AnnotatorDistanceCheck3End(AnnotatorDistanceCheck):
    """unequally sized intergenic space, over-representation of 3' end."""
    
    def setUp(self):
        AnnotatorDistanceCheck.setUp(self)

        work_outfile = open( self.workspace, "w" )
        segs_outfile = open( self.segments, "w" )

        w = GTF.Entry()
        w.contig, w.strand, w.gene_id, w.transcript_id = "chr1", "+", "gene1", "trans1"
        e = GTF.Entry()
        e.contig, e.strand  = "chr1", "+"

        # 10kb genes, size of intergenic space grows by random increment
        x, y = 0, 0
        w_inc = 0
        while x < 10000000:
            
            w.gene_id, w.transcript_id = "gene%i" % x, "trans1"
            e.start, e.end = x, x + 10000

            work_outfile.write( str(e) + "\n" )
            
            x += 10000
            w_inc += random.randint( 0, 10000)
            end = x + w_inc 
            y = x
            s_inc = 0
            while y < end:
                e.gene_id, e.transcript_id = "gene%i" % (y), "trans1"
                e.start, e.end = y, y+random.randint( 50, 150 )
                segs_outfile.write( str(e) + "\n" )
                y += s_inc
                s_inc += random.randint( 0, 100 )

            x = end

        work_outfile.close()
        segs_outfile.close()

    def testMain( self ):
        annotator_distance.main( ["--workspace=%s" % self.workspace, 
                                  "--segments=%s" %  self.segments,
                                  "--workspace-labels=direction", 
                                  "--workspace-builder=gtf-intergenic" ])


if __name__ == "__main__":
    unittest.main()
