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
align_pairs.py - align nucleotide sequence pairs
================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Read a file of unaligned nucleotide sequence pairs and align them.

Usage
-----

Example::

   python align_pairs.py --help

Type::

   python align_pairs.py --help

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
import tempfile
import optparse
import CGAT.Experiment as E
import CGAT.AlignedPairs as AlignedPairs
import CGAT.Blat as Blat
import CGAT.FastaIterator as FastaIterator
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IOTools as IOTools

def iterate_double_fasta ( fn1, fn2 ):
    iterator = FastaIterator.iterate_together( fn1, fn2 )
    for seq1, seq2 in iterator:
        yield AlignedPairs.UnalignedPair( 
            token1 = seq1.title,
            sequence1 = seq1.sequence,
            token2 = seq2.title,
            sequence2 = seq2.sequence )

def iterate_single_fasta ( fn1 ):
    iterator = FastaIterator.FastaIterator( fn1 )
    while 1:
        seq1, seq2 = iterator.next(), iterator.next()
        yield AlignedPairs.UnalignedPair( 
            token1 = seq1.title,
            sequence1 = seq1.sequence,
            token2 = seq2.title,
            sequence2 = seq2.sequence )

def iterate_list ( infile, idx1, idx2 = None):

    fasta1 = IndexedFasta.IndexedFasta( idx1 )
    if idx2 == None:
        fasta2 = fasta1
    else:
        fasta2 = IndexedFasta.IndexedFasta( idx2 )

    first = True
    for line in infile:
        if line[0] == "#": continue
        id1, id2 = line[:-1].split("\t")[:2]

        try:
            yield AlignedPairs.UnalignedPair( 
                token1 = id1,
                sequence1 = fasta1.getSequence( id1 ),
                token2 = id2,
                sequence2 = fasta2.getSequence( id2 ) )
        except KeyError, msg:
            if first:
                first = False
                continue
            raise KeyError(msg)
    
def iterate_table( infile ):
    for line in infile:
        if line[0] == "#": continue
        
        unaligned_pair = AlignedPairs.UnalignedPair()
        unaligned_pair.Read(line)
        yield unaligned_pair


def getFile( section, options ):
    if options.output_filename_pattern:
        if "%s" in options.output_filename_pattern:
            return open( options.output_filename_pattern % section, "w" )
        else:
            return open( options.output_filename_pattern, "w" )
    return options.stdout
        
##------------------------------------------------------------
if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: align_pairs.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("--skip-statistics", dest="skip_stats", action="store_true",
                      help="do not compute alignment statistics [%default]."  )

    parser.add_option("--method", dest="methods", type="choice", action="append",
                      choices=("dialign", "clustal", "blastz", "nw", "sw", "dba", "dialignlgs" ),
                      help="alignment method [%default]."  )

    parser.add_option("--anchor-alignment", dest="anchor_alignment", type="int",
                      help="anchor alignmet with xxx residues [%default]."  )

    parser.add_option("--output-format", dest="output_formats", type="choice", action="append",
                      choices=("fasta", "stats", "psl" ),
                      help="anchor alignment with xxx residues [%default]."  )

    parser.add_option("--input-format", dest="input_format", type="choice", 
                      choices=("fasta", "list" ),
                      help="input format of stdin [%default]."  )

    parser.add_option("--output-filename-pattern", dest="output_filename_pattern", type="string",
                      help="output pattern for multiple files [%default]."  )

    parser.add_option("--filename-sequences1", dest="filename_sequences1", type="string",
                      help="first indexed input filename with sequences [%default]."  )

    parser.add_option("--filename-sequences2", dest="filename_sequences2", type="string",
                      help="second indexed input filename with sequences [%default]."  )

    parser.add_option("--options-blastz", dest="options_blastz", type="string",
                      help="command line options for blastz [%default]."  )

    parser.set_defaults( 
        skip_stats = False,
        methods = [],
        output_formats = [],
        input_format = "fasta",
        output_filename_pattern = None,
        filename_sequences1 = None,
        filename_sequences2 = None,
        anchor_alignment = 0,
        options_blastz = "C=2 B=1 T=0 W=6 K=2200" )

    (options, args) = E.Start( parser, add_pipe_options = True )

    if len(options.methods) == 0:
        print USAGE
        print "please specify an alignment method."
        sys.exit(1)

    if len(options.output_formats) == 0:
        print USAGE
        print "please specify at least one output format."
        sys.exit(1)

    if len(args) == 2:
        iterator = iterate_double_fasta( args[0], args[1] )        
    elif options.filename_sequences1 and options.filename_sequences2:
        if len(args) == 0 or (len(args) == 1 and args[0] == "-"):
            infile = options.stdin
        elif len(args) == 1:
            infile = open( args[0], "r") 
                
        iterator = iterate_list( infile, options.filename_sequences1, options.filename_sequences2 )
    else:
        iterator = iterate_single_fasta( options.stdin )
        

    npairs, ntoken_pairs = 0, 0
    ninput, nskipped, nerrors = 0, 0, 0

    outfile_table = None
    outfile_fasta = None
    outfile_psl = None
    if "table" in options.output_formats:
        outfile_table = getFile( "table ", options )
        outfile_table.write( """# CATEGORY:       category [intron|exon]
# METHOD:         alignment method
# TOKEN:          name
# ID:             segment id
# TOTAL:          number of segments
# LEN:            length of segment
# NALIGNED:       number of aligned positions
# PALIGNED:       percentage of aligned positions
# IDENT:          number of identical positions
# TRANSIT:        number of transitions
# TRANSVERS:      number of transversion
# MATCHES:        number of matching positions
# PIDENT:         percentage of identical positions
# PTRANSIT:       precentage of transitions
# PTRANSVERS:     precentage of transversion
# BLOCKSIZES:     alignment, length of blocks
# GAPS:           gap sizes in sequence 1/2
CATEGORY\tMETHOD\tTOKEN1\tID1\tTOTAL1\tLEN1\tTOKEN2\tID2\tTOTAL2\tLEN2\tNALIGNED\tPALIGNED\tIDENT\tTRANSIT\tTRANSVER\tMATCHES\tPIDENT\tPTRANSVIT\tPTRANVER\tBLOCKSIZES\tGAPSIZES\tGAPSIZES\tTYPE1\tTYPE2\n""")

    if "fasta" in options.output_formats:
        outfile_fasta = getFile( "fasta", options )

    if "psl" in options.output_formats:
        outfile_psl = getFile( "psl", options )

    ## setup alignment objects
    for unaligned_pair in iterator:

        ninput += 1
        
        for method in options.methods:

            pair = AlignedPairs.AlignedPair( unaligned_pair )
            pair.mOptionsBlastZ = options.options_blastz

            try:
                pair.Align( method, anchor = options.anchor_alignment )
            except AlignedPairs.AlignmentError, msg:
                
                if options.loglevel >= 1:
                    options.stdlog.write( "# %s - %s: %s\n" % (msg, unaligned_pair.mToken1, unaligned_pair.mToken2))
                    if options.loglevel >= 2:
                        options.stdlog.write( "# input=%s\n" % (str(unaligned_pair)))

                nskipped += 1
                continue

            if outfile_table:
                outfile_table.write( str(pair) + "\n" )
            
            if outfile_fasta:
                outfile_fasta.write( ">%s\n%s\n>%s\n%s\n" % (pair.mToken1, pair.mAlignedSequence1, pair.mToken2, pair.mAlignedSequence2 ) )

            if outfile_psl:
                entry = Blat.Match()
                entry.mQueryId, entry.mSbjctId = pair.mToken1, pair.mToken2
                entry.strand = pair.strand
                entry.fromMap( pair.mAlignment )
                outfile_psl.write( str(entry) + "\n" )

            npairs += 1
     
    E.info( "input=%i, skipped=%i, nerrors=%i, transcripts=%i, pairs=%i" % (ninput, nskipped,
                                                                            nerrors,
                                                                            ntoken_pairs, npairs ) )
    E.Stop()
    sys.exit(0)

