################################################################################
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
#################################################################################
'''
gff_ensembl2gbrowser.py - 
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

   python gff_ensembl2gbrowser.py --help

Type::

   python gff_ensembl2gbrowser.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import getopt
import time
import sets
import optparse
import math
import tempfile

import CGAT.Experiment as E

USAGE="""

convert an ensembl gff file into a gbrowser gff file for upload.

Input:
X       ensembl CDS     60708731        60709127        .       -       0        gene_id "ENSMODG00000009707"; transcript_id "ENSMODT00000038508"; exon_number "4"; protein_id "ENSMODP00000036912";
Output:
chr2    gpipe   gene    211249748       211254644       .       -       .       Id=ENSMODG00000009707
chr3    gpipe   mRNA    362769650       362870268       45.24   -       .       Id=ENSMODG00000009707.ENSMODT00000038508 ; Parent=ENSMODG00000009707
chr8    gpipe   CDS     202417740       202417805       59.0909090909   -       0       ID=ENSMODG00000009707.ENSMODT00000038508.4 ; Parent=ENSMODG00000009707.ENSMODT00000038508

"""

parser = E.OptionParser( version = "%prog version: $Id: gff_ensembl2gbrowser.py 2781 2009-09-10 11:33:14Z andreas $")

if __name__ == "__main__":

    parser.add_option( "-s", "--source-id", dest="source_id", type="string",
                       help="source_id to use." )
    
    parser.add_option( "-f", "--feature", dest="feature", type="string",
                       help="name of feature (CDS, exon, ...)." )

    parser.set_defaults(
        source_id = "ensembl",
        feature = "CDS",
        )

    (options, args) = E.Start( parser, add_csv_options  = True)

    input_gffs = GTF.readFromFile( sys.stdin, separator = " " )

    ## sort by genes
    input_gffs.sort( lambda x,y: cmp( (x["gene_id"], x["transcript_id"], x.start), (y["gene_id"], y["transcript_id"], y.start) ) )

    new_gffs = []

    gene_from = None
    gene_to = None
    last_g = None
    mrna_from = None
    mrna_to = None    
    last_mrna = None

    for g in input_gffs:

        if g.feature != options.feature: continue

        if last_g:
            if last_g["transcript_id"] != g["transcript_id"]:
                n = GTF.Entry()
                n.Fill(last_g)
                n.clearAttributes()                
                n.feature = "mRNA"
                n.start = mrna_from
                n.end = mrna_to
                n.frame = "."
                n.addAttribute( "mRNA", "%s.%s" % (last_g["gene_id"], last_g["transcript_id"]), separator = " " )
                n.addAttribute( "Gene", "%s" % (last_g["gene_id"]), separator = " ")
                new_gffs.append( n )
                mrna_from, mrna_to = g.start, g.end
                
            if last_g["gene_id"] != g["gene_id"]:
                n = GTF.Entry()
                n.Fill(last_g)
                n.clearAttributes()
                n.feature = "gene"                
                n.start = gene_from
                n.end = gene_to
                n.frame = "."
                new_gffs.append( n )
                gene_from, gene_to = g.start, g.end
                n.addAttribute( "Gene", "%s" % (last_g["gene_id"]), separator = " " )
                
        else:
            gene_from, gene_to = g.start, g.end
            mrna_from, mrna_to = gene_from, gene_to
                
        n = GTF.Entry()            

        n.Fill( g )
        n.clearAttributes()
        
        n.addAttribute( "mRNA", "%s.%s" % (g["gene_id"], g["transcript_id"] ), separator = " ")

        new_gffs.append( n )

        gene_from, gene_to = min( g.start, gene_from), max(g.end, gene_to)
        mrna_from, mrna_to = min( g.start, mrna_from), max(g.end, mrna_to)        

        last_g = g
        
    ## combine gffs to genes and dump:
    
    for n in new_gffs:
        n.source = options.source_id
        options.stdout.write(str(n) + "\n")
        
    E.Stop()
