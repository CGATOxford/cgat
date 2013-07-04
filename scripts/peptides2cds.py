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
peptides2cds.py - map peptide sequences onto nucleic acid sequences
===================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

build a map between peptides to cdnas. This script deals with
frameshifts but not with larger indels. 

This script requires the :mod:`alignlib` library and :file:`muscle`
to be installed.

Usage
-----

Example::

   python peptides2cds.py --help

Type::

   python peptides2cds.py --help

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
import tempfile
import time
import optparse
import CGAT.Experiment as E
import CGAT.BlastAlignments as BlastAlignments
import alignlib
import CGAT.Genomics as Genomics
import CGAT.FastaIterator as FastaIterator
import CGAT.IOTools as IOTools

def AlignExhaustive( seq_wobble, seq_cds, seq_peptide, map_p2c, options, diag_width = 2 ):
    """Align two sequences.

    Align in chunks to keep memory low. Both sequences are roughly the same,
    thus align only in diagonal.
    """

    gop, gep = -1.0, -1.0
    matrix = alignlib.makeSubstitutionMatrixBackTranslation( 1, -10, 1, alignlib.getDefaultEncoder() )
    alignlib.setDefaultSubstitutionMatrix( matrix )

    if seq_wobble.getLength() < 10000:
        if options.loglevel >= 6:
            options.stdlog.write( "# using full dynamic programing matrix.\n" )
            options.stdlog.flush()
        # do not penalize gaps at the end, because sometimes the last codon might be missing
        alignator = alignlib.makeAlignatorDPFull( alignlib.ALIGNMENT_GLOBAL, gop, gep, 1, 1 )
    else:
        diag_width = abs(seq_wobble.getLength() - seq_cds.getLength()) + 1
        if options.loglevel >= 6:
            options.stdlog.write( "# using dot alignment with diagonal %i\n" % diag_width )
            options.stdlog.flush()

        dots = alignlib.makeAlignmentMatrixRow()
            
        for x in range(0, seq_wobble.getLength()):
            xr = seq_wobble.asResidue( x )
            for y in range( max(0, x - diag_width), min( seq_cds.getLength(), x + diag_width)):
                s = matrix.getValue( xr, seq_cds.asResidue(y) )
                if s >= 0:
                    dots.addPair( x, y, float(s) )

        if options.loglevel >= 6:
            options.stdlog.write( "# finished adding %i dots" % dots.getLength() )
            options.stdlog.flush()

        alignator_dummy = alignlib.makeAlignatorPrebuilt( dots )
        
        alignator = alignlib.makeAlignatorDots( alignator_dummy, gop, gep )

    alignator.align( map_p2c, seq_wobble, seq_cds )

def FindMatch( matrix, seq1, x, seq2, y, options, max_advance = 10, min_match = 3 ):
    """advance on sequence2 from position y until you find a minimum match of
    min_match residues.

    returns the offset in y. If no match is found, an offset of 0 is returned.
    """

    #for z in range(0,10):
    #print seq1.asChar(x+z), seq2.asChar(y+z)

    min_match = min( seq1.getLength() - x, min_match )
    l2 = seq2.getLength()
    # get list of reference residues
    rr = []
    for xx in range(x, x + min_match ):
        rr.append( seq1.asResidue(xx) )

    max_advance = min( max_advance, l2 - y)
    for yy in range(y, y + max_advance):
        m = True
        for i in range( 0, min_match ):
            if matrix.getValue( rr[i], seq2.asResidue(yy + i) ) < 0:
                m = False
                break
        if m: return yy - y

    return 0

def AlignCodonBased( seq_wobble, seq_cds, seq_peptide, map_p2c, options,
                     diag_width = 2, max_advance = 2 ):
    """advance in codons in seq_wobble and match to nucleotides in seq_cds.

    Due to alinglib this is all in one-based coordinates.
    Takes care of frameshifts.
    """
    
    map_p2c.clear()

    gop, gep = -1.0, -1.0
    matrix = alignlib.makeSubstitutionMatrixBackTranslation( 1, -10, 1, alignlib.getDefaultEncoder() )

    pep_seq = seq_peptide.asString()
    cds_seq = seq_cds.asString()
    wobble_seq = seq_wobble.asString()
    
    lcds = seq_cds.getLength()
    lwobble = seq_wobble.getLength()
    y = 0
    x = 0

    last_start = None

    while x < lwobble and y < lcds:

        xr = seq_wobble.asResidue( x )
        # skip over masked chars in wobble - these are gaps
        if seq_wobble.asChar(x) == "X": 
            x += 1
            continue

        # skip over masked chars in wobble - these are from
        # masked chars in the peptide sequence
        # Note to self: do not see all implications of this change
        # check later.
        if seq_wobble.asChar(x) == "N": 
            x += 1
            continue

        # skip over gaps in wobble 
        if seq_wobble.asChar(x) == "-": 
            x += 1
            continue

        s = matrix.getValue( xr, seq_cds.asResidue(y) )

        if options.loglevel >= 6:
            if (x % 3 == 0):
                c = seq_cds.asChar(y) + seq_cds.asChar(y+1) + seq_cds.asChar(y+2)
                options.stdlog.write( "# c=%s, x=%i, y=%i, aa=%s target=%s\n" % (c, x, y,
                                                                                 Genomics.MapCodon2AA( c ),
                                                                                 pep_seq[int(x/3)]) )
                                      
            options.stdlog.write( "# x=%i\twob=%s\ty=%i\tcds=%s\txr=%s\tcds=%i\tscore=%s\n" % \
                                      (x, seq_wobble.asChar(x), y, seq_cds.asChar(y), xr, seq_cds.asResidue(y), str(s) ))
            
        # deal with mismatches
        if s <= 0:

            tmp_map_p2c = alignlib.makeAlignmentVector()

            ## backtrack to previous three codons and align
            ## three codons for double frameshifts that span two codons and
            ## produce two X's and six WWWWWW.

            ## number of nucleotides to extend (should be multiple of 3)
            ## less than 12 caused failure for some peptides.
            d = 15
            
            # extend by amound dx
            dx = (x % 3) + d
            
            x_start = max(0, x - dx )
            # map to ensure that no ambiguous residue mappings
            # exist after re-alignment
            y_start = max(0, map_p2c.mapRowToCol( x_start, alignlib.RIGHT ))

            if (x_start, y_start) == last_start:
                raise ValueError( "infinite loop detected" )

            last_start = (x_start, y_start)

            x_end = min(x_start + 2 * d, len(wobble_seq) )
            y_end = min(y_start + 2 * d, len(cds_seq) )

            wobble_fragment = alignlib.makeSequence(wobble_seq[x_start:x_end])
            cds_fragment = alignlib.makeSequence(cds_seq[y_start:y_end])
            
            AlignExhaustive( wobble_fragment, cds_fragment, "", tmp_map_p2c, options )

            if options.loglevel >= 10:
                 options.stdlog.write("# fragmented alignment from %i-%i, %i-%i:\n%s\n" % (x_start, x_end,
                                                                                           y_start, y_end,
                                                                                           str(alignlib.AlignmentFormatExplicit( tmp_map_p2c,
                                                                                                                                 wobble_fragment, 
                                                                                                                                 cds_fragment ))))
                 
                 options.stdlog.flush()

            ## clear alignment
            map_p2c.removeRowRegion( x_start, x_end )
            ngap = 0
            last_x, last_y = None, None
            for xxx in range( tmp_map_p2c.getRowFrom(), tmp_map_p2c.getRowTo() ):
                yyy = tmp_map_p2c.mapRowToCol(xxx)

                if yyy >= 0:
                    x = xxx + x_start
                    y = yyy + y_start
                    xr = seq_wobble.asResidue(x)
                    s = matrix.getValue( seq_wobble.asResidue(x), seq_cds.asResidue(y) )
                    if s < 0:
                        raise ValueError("mismatched residue wobble: %i (%s), cds: %i (%s)" % (x, seq_wobble.asChar(x), y, seq_cds.asChar(y)))
                    
                    map_p2c.addPair( x, y, s)
                    last_x, last_y = x, y
                    if options.loglevel >= 6:
                        options.stdlog.write( "# reset: x=%i\twob=%s\ty=%i\tcds=%s\txr=%s\tcds=%i\tscore=%i\n" % \
                                              (x, seq_wobble.asChar(x), y, seq_cds.asChar(y), xr, seq_cds.asResidue(y), s ))
                        options.stdlog.flush()
                    ngap = 0
                else:
                    ngap += 1

                # treat special case of double frameshifts. They might cause a petide/wobble residue
                # to be eliminated and thus the translated sequences will differ.
                # simply delete the last residue between x and y and move to next codon.
                if ngap == 3:
                    map_p2c.removeRowRegion( last_x, last_x + 1 )

                    last_x += 1
                    map_p2c.addPair( last_x, last_y )
                    if options.loglevel >= 6:
                        options.stdlog.write( "# double: x=%i\twob=%s\ty=%i\tcds=%s\txr=%s\tcds=%i\tscore=%i\n" % \
                                              (last_x, seq_wobble.asChar(last_x), last_y, seq_cds.asChar(last_y), xr, seq_cds.asResidue(last_y), s ))
                        options.stdlog.flush()                    
                    ngap = 0
                    
            ## exit condition if alignment is shorter than problematic residue
            ## need to catch this to avoid infinite loop.
            if tmp_map_p2c.getRowTo() < d:
                if lwobble - x <= 4:
                    ## only last codon is missing, so ok
                    break
                else:
                    raise ValueError("failure to align in designated window.")
                    
            s = 0
            
        s = matrix.getValue( xr, seq_cds.asResidue(y) )

        if s < 0:
            raise ValueError("mis-matching residues.")
        
        map_p2c.addPair( x, y, float(s) )
        
        # advance to next residues
        x += 1
        y += 1

    # sanity checks
    assert( map_p2c.getRowTo() <= seq_wobble.getLength() )
    assert( map_p2c.getColTo() <= seq_cds.getLength() )

def PrintPrettyAlignment( seq_wobble, seq_cds, seq_pep, map_p2c, options ):
    """print a pretty alignment."""

    f = alignlib.AlignmentFormatExplicit( map_p2c, seq_wobble, seq_cds )
    wobble_ali, cds_ali = f.mRowAlignment, f.mColAlignment
    
    wi, ci, pi = 0, 0, 0
    frags_w, frags_c, frags_p = [], [], []
    for x in range( 0, len( wobble_ali )):

        if wi % 3 == 0:
            if pi < len(seq_pep):
                frags_p.append( "  %s " % seq_pep[pi])
            frags_w.append (" ")
            frags_c.append (" ")
            pi += 1

        frags_w.append( wobble_ali[x] )
        frags_c.append( cds_ali[x] )
        if wobble_ali[x] != "-":
            wi += 1
            
        if len(frags_w) > 120 and len(frags_w) % 3 == 0:
            options.stdlog.write( "#" + "".join(frags_w) + "\n" )
            options.stdlog.write( "#" + "".join(frags_p) + "\n" )
            options.stdlog.write( "#" + "".join(frags_c) + "\n" )                
            options.stdlog.write( "#\n" )
            frags_w, frags_c, frags_p = [], [], []


    options.stdlog.write( "#" + "".join(frags_w) + "\n" )
    options.stdlog.write( "#" + "".join(frags_p) + "\n" )
    options.stdlog.write( "#" + "".join(frags_c) + "\n" )                
    options.stdlog.write( "#\n" )
    

##--------------------------------------------------------------------------
##--------------------------------------------------------------------------
##--------------------------------------------------------------------------
##--------------------------------------------------------------------------        
def getMapPeptide2Cds( peptide_sequence, cds_sequence, options ):
    """get map between peptide sequence and cds sequence.
    
    The returned alignment is in nucleotides.

    """
    
    ## remove whitespaces form protein sequence
    p = re.sub(" ", "", peptide_sequence )

    ## remove gaps and whitespaces from cds
    c = re.sub("[ .-]", "", cds_sequence )

    w = Genomics.Protein2Wobble( p.upper() )

    if options.loglevel >= 6:
        options.stdlog.write( "# peptide original (%5i): %s\n" % (len(p), p) )
        options.stdlog.write( "# cds original     (%5i): %s\n" % (len(c), c) )
        options.stdlog.write( "# wobble sequence  (%5i): %s\n" % (len(w), w) )
        options.stdlog.flush()

    seq_wobble = alignlib.makeSequence( w )
    seq_cds = alignlib.makeSequence( string.upper(c) )
    seq_peptide = alignlib.makeSequence( p )

    map_p2c = alignlib.makeAlignmentVector()

    try:
        AlignCodonBased( seq_wobble, seq_cds, seq_peptide, map_p2c, options = options )
    except ValueError, msg:
        raise ValueError( "mapping error for sequence: %s" % (msg) )

    ## if there are more than five frameshifts - do exhaustive alignment
    max_gaps = 5
    num_peptide_gaps = len( re.sub("[^-]", "", p ) )
    ngaps = map_p2c.getNumGaps() - (num_peptide_gaps *  3) - abs(len(w)-len(c))
    
    if options.loglevel >= 6:
        options.stdlog.write("# alignment between wobble and cds: ngaps=%i, npeptide_gaps=%i\n" % (ngaps, num_peptide_gaps) )
        PrintPrettyAlignment( seq_wobble, seq_cds, p, map_p2c, options )

    if ngaps > max_gaps:
        if options.loglevel >= 2:
            options.stdlog.write("# too many gaps (%i>%i), realigning exhaustively.\n" % (ngaps, max_gaps ) )
            options.stdlog.flush()
        full_map_p2c = alignlib.makeAlignmentVector()
        
        AlignExhaustive( seq_wobble, seq_cds, seq_peptide, full_map_p2c, options )
        if options.loglevel >= 6:
            options.stdlog.write("# full alignment between wobble and cds:\n" )
            options.stdlog.flush()
            PrintPrettyAlignment( seq_wobble, seq_cds, p, full_map_p2c, options )

        map_p2c = full_map_p2c
        
    ## remove incomplete codons
    x = 0
    while x < len(p) * 3:
        if (map_p2c.mapRowToCol( x ) < 0 or \
            map_p2c.mapRowToCol( x+1 ) < 0 or \
            map_p2c.mapRowToCol( x+2 ) < 0 ):
            map_p2c.removeRowRegion( x, x + 3)
        x += 3

    if map_p2c.getLength() == 0:
        if options.loglevel >= 1:
            options.stdlog.write("# WARNING: empty alignment\n" )
            if options.loglevel >= 6:
                options.stdlog.write( "# peptide original: %s\n" % p )
                options.stdlog.write( "# cds original    : %s\n" % c )
                options.stdlog.write( "# wobble sequence : %s\n" % w )

        raise ValueError( "empty alignment" )

    assert( map_p2c.getRowTo() <= seq_wobble.getLength() )
    assert( map_p2c.getColTo() <= seq_cds.getLength() )

    return map_p2c
    
        
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: peptides2cds.py 2890 2010-04-07 08:58:54Z andreas $")

    parser.add_option("-p", "--peptides", dest="filename_peptides", type="string",
                      help="filename with peptide sequences [%default]."  )
    
    parser.add_option("-c", "--cds", "--cdnas", dest="filename_cdna", type="string",
                      help="filename with cdna sequences [%default]."  )

    parser.add_option("-m", "--map", dest="filename_map", type="string",
                      help="filename with map of peptide identifiers to cdna identifiers [%default]."  )

    parser.add_option( "--output-identifier", dest="output_identifier", type="choice",
                       choices=("cdna", "peptide"),
                       help="output identifier to use [%default]."  )


    parser.add_option("-f", "--output-format=", dest="output_format", type="choice",
                      choices=("alignment", "fasta"),
                      help="output format.")
    
    parser.set_defaults(
        peptides=None,
        filename_cdna = None,
        output_format="alignment",
        filename_map = None,
        stop_codons = ("TAG", "TAA", "TGA"),
        output_identifier = "peptide",
        )

    (options, args) = E.Start( parser, add_pipe_options = True )

    if not options.filename_cdna:
        raise ValueError("please supply filename with cds sequences.")

    if options.filename_peptides:
        infile = open(options.filename_peptides, "r") 
        E.info("reading from %s" % options.filename_peptides)
    else:
        E.info("reading from stdin")
        infile = sys.stdin

    if options.filename_map:
        E.info( "reading map" )
        map_peptide2cds = IOTools.readMap( IOTools.openFile( options.filename_map, "r" ) )
        E.info( "read map for %i identifiers" % len(map_peptide2cds) )
    else:
        map_peptide2cds = {}

    E.info( "reading cds sequences" )
        
    cds_sequences = Genomics.ReadPeptideSequences( IOTools.openFile(options.filename_cdna, "r") )    

    E.info( "read %i cds sequences" % len(cds_sequences))

    ninput, noutput = 0, 0
    nskipped, nnosequence = 0, 0

    # iterate over peptide sequences
    iterator = FastaIterator.FastaIterator( infile )

    use_cds_id = options.output_identifier == "cds"

    for cur_record in iterator:

        ninput += 1
        
        peptide_identifier = re.split("\s+", cur_record.title)[0]
        cds_identifier = map_peptide2cds.get( peptide_identifier, peptide_identifier )

        if cds_identifier not in cds_sequences:
            nnosequence += 1
            continue

        p = cur_record.sequence
        c = cds_sequences[cds_identifier]
        
        E.debug("processing %s: laa=%i (without gaps=%i), lna=%i" % (peptide_identifier, len(p), len(re.sub("-", "", p)), len(c)))

        try:
            map_p2c = getMapPeptide2Cds( p, c, options )
        except ValueError:
            nskipped += 1
            continue
            
        if use_cds_id:
            identifier = cds_identifier
        else:
            identifier = peptide_identifier

        if options.output_format =="alignment":
            options.stdout.write("\t".join( map(str, (identifier, alignlib.AlignmentFormatEmissions( map_p2c ),
                                                      len(cur_record.sequence), len(cds_sequences[identifier])) ) )+"\n")
            
        elif options.output_format == "fasta":

            map_p2c.switchRowCol()

            alignatum = alignlib.makeAlignatum( c )
            
            alignatum.mapOnAlignment( map_p2c, len(p) * 3 )

            s = alignatum.getString()
            if len(s) != len(p) * 3:
                raise ValueError ("incomplete aligned string for %s: %s, cds=%s" % (cur_record.title, s, c ))
            
            options.stdout.write( ">%s\n%s\n" % (identifier, s ))

        noutput += 1
        sys.stdout.flush()

    E.info( "ninput=%i, noutput=%i, nnosequence=%i, nskipped=%i" % (ninput, noutput, nnosequence, nskipped) )
        
    E.Stop()
            




            
