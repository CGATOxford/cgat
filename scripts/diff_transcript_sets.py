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
diff_transcript_sets.py - 
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

   python diff_transcript_sets.py --help

Type::

   python diff_transcript_sets.py --help

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

""" program $Id: diff_transcript_sets.py 2781 2009-09-10 11:33:14Z andreas $

show differences between two transcript sets.
"""
import CGAT.Experiment as E
import CGAT.IOTools as IOTools

def parseIdentifier( id, options ):
    data = id.split(options.separator)

    if len(data) == 4:
        species, transcript, gene = data[:3]
    elif len(data) == 2:
        species, gene = data
        transcript = gene
    else:
        options.stderr.write("# parsing error for '%s'\n" % id )
        raise ValueError
    
    return species, transcript, gene

def countGenesTranscripts( inlist, options ):
    """count number of genes/transcripts in list."""

    genes = {}
    transcripts = {}
    
    for x in inlist:

        try:
            species, transcript, gene = parseIdentifier( x, options )
        except ValueError:
            continue
        
        if species not in genes:
            genes[species] = set()
            transcripts[species] = set()

        transcripts[species].add(transcript)
        genes[species].add(gene)
        
    return genes, transcripts

def getTranscriptsForGenes( genes, transcripts, options ):
    """get transcripts for list of genes."""

    result = []
    for x in transcripts:
        try:
            species, transcript, gene = parseIdentifier( x, options )
        except ValueError:
            continue
        
        if gene in genes:
            result.append( x )
    return result

##--------------------------------------------------------------------------
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: diff_transcript_sets.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-p", "--add-percent", dest="add_percent", action="store_true",
                      help="add percent columns" )
    parser.add_option("-d", "--dump-sets", dest="dump_sets", action="append", type="choice",
                      choices=("rest_genes1", "rest_genes2", "intersection", "union"),
                      help="dump sets of transcripts/genes" )
    parser.add_option("-o", "--output-pattern", dest="output_pattern", type="string",
                      help="output pattern to use for dumped sets. Should contain one %s." )
    
    parser.set_defaults(
        separator = "|",
        add_percent = "False",
        dump_sets = [],
        output_pattern = "%s",
        )

    (options, args) = E.Start( parser, add_psql_options = True )

    options.filename1, options.filename2 = args

    ids1, nerrors1 = IOTools.ReadList( open(options.filename1, "r") )
    ids2, nerrors2  = IOTools.ReadList( open(options.filename2, "r") )

    genes1, transcripts1 = countGenesTranscripts( ids1, options )
    genes2, transcripts2 = countGenesTranscripts( ids2, options )    

    options.stdout.write( "species\tngenes1\tntranscripts1\tngenes2\tntranscripts2\ttr_inter\ttr_union\ttr_rest1\ttr_rest2\ttr_inter\tg_union\tg_rest1\tg_rest2" )
    options.stdout.write( "\ttr_rest1\ttr_rest2\tg_rest1\tg_rest2")
    
    options.stdout.write ("\n" )
    
    for species in set(genes1.keys()).union(set(genes2.keys())):
        nt1, nt2, ng1, ng2 = "na", "na", "na", "na"
        
        if species in genes1:
            g1 = genes1[species]
            t1 = transcripts1[species]
            nt1 = "%i" % len(transcripts1[species])
            ng1 = "%i" % len(genes1[species])
        else:
            t1, g1 = None, None
        
        if species in genes2:
            g2 = genes2[species]
            t2 = transcripts2[species]            
            nt2 = "%i" % len(transcripts2[species])
            ng2 = "%i" % len(genes2[species])
        else:
            t2, g2 = None, None
            
        if species in transcripts1 and transcripts2:
            ct = "%i" % len(t1.intersection(t2))
            ut = "%i" % len(t2.union(t1))
            rt1 = "%i" % len(t1.difference(t2))
            rt2 = "%i" % len(t2.difference(t1))
        else:
            ct, ut, rt1, rt2 = ["na"] * 4

        if species in genes1 and genes2:
            cg = "%i" % len(g1.intersection(g2))
            ug = "%i" % len(g2.union(g1))
            rg1 = "%i" % len(g1.difference(g2))
            rg2 = "%i" % len(g2.difference(g1))
        else:
            cg, ug, rg1, rg2 = ["na"] * 4
            
        options.stdout.write( "\t".join( (species, nt1, ng1, nt2, ng2)) )
        options.stdout.write( "\t" )
        options.stdout.write( "\t".join( (ct, ut, rt1, rt2)) )        
        options.stdout.write( "\t" )        
        options.stdout.write( "\t".join( (cg, ug, rg1, rg2)) )

        if options.add_percent:
            if species in genes1 and genes2:
                rg1 = "%5.2f" % (100.0 * len(g1.difference(g2)) / len(g1))
                rg2 = "%5.2f" % (100.0 * len(g2.difference(g1)) / len(g2))
            if species in transcripts1 and transcripts2:
                rt1 = "%5.2f" % (100.0 * len(t1.difference(t2)) / len(t1))
                rt2 = "%5.2f" % (100.0 * len(t2.difference(t1)) / len(t2))
            options.stdout.write( "\t" )                        
            options.stdout.write( "\t".join( (rt1, rt2, rg1, rg2)) )
        
        options.stdout.write( "\n" )

        for choice in options.dump_sets:

            output_set = None
            
            if choice == "rest_genes1" and g1 and g2:
                output_set = getTranscriptsForGenes(g1.difference(g2), ids1, options )
                
            elif choice == "rest_genes2" and g1 and g2:
                output_set = getTranscriptsForGenes(g2.difference(g1), ids2, options )
                
            if output_set:
                outfile = open( options.output_pattern % (choice), "w" )
                for x in output_set:
                    outfile.write( "%s\n" % (x,) )
                outfile.close()
                    
    E.Stop()
