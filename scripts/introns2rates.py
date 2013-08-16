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
introns2rates.py - compute rates of aligned intron sequences
============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Read a list of aligned introns and compute pairs.

Options::

    -h, --help                      print this message.
    -v, --verbose=                  loglevel.
    --is-compressed                 alignments are in compressed format (unaligend sequence has to
                                    be provided)
    --do-gblocks                    apply gblocks to alignment before calculation.
    --skip-distance                 skip distance calculation
    --skip-alistats                 skip alignment characteristics calculation
    --echo-unaligned                echo unaligned transcripts
    --method=                       alignment method [dbaligned|clusaligned|dialigned|dialignedlgs]
    --fixed-alpha=                  fix alpha at value
    --anchor-alignment=             anchor alignment with x residues of A at each side

Usage
-----

Example::

   python introns2rates.py --help

Type::

   python introns2rates.py --help

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

import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.Exons as Exons
import CGAT.WrapperBaseML as WrapperBaseML
import CGAT.WrapperGblocks as WrapperGblocks
import alignlib

import CGAT.WrapperDialign as WrapperDialign
import CGAT.WrapperDBA as WrapperDBA
# clustalw wrapper not up-to-date
# import CGAT.WrapperClustal as WrapperClustal

param_long_options=["verbose=", "help",
                    "write-exons=", "write-introns=",
                    "extend-introns=",
                    "is-compressed", "do-gblocks",
                    "skip-distance", "skip-alistats", "echo-unaligned",
                    "method=", "fixed-alpha=",
                    "anchor-alignment=", "version"]

param_short_options="v:hg:"

param_loglevel = 2

param_echo_unaligned = 0
param_is_compressed = 0

param_do_gblocks = 0

param_do_distance = 1

param_intron_types = ( ( "U2-GT/AG", "GT", "AG"),
                       ( "U2-nc-GC/AG", "GC", "AG"),
                       ( "U12-AT/AC", "AT", "AC") )

param_method = None

param_do_alignment = False
param_do_alistats = True

param_fixed_alpha = None

param_anchor_alignment = 0

class IntronPair:
    def __init__(self):
        self.mCategory = ""
        self.mMethod = ""
        self.mToken1 = ""
        self.mIntronId1 = ""
        self.mNumIntrons1 = 0
        self.mLenIntron1 = 0
        self.mToken2 = ""
        self.mIntronId2 = ""
        self.mNumIntrons2 = 0
        self.mLenIntron2 = 0
        self.mIdentical = 0
        self.mTransitions = 0
        self.mTransversions = 0
        self.mNumGaps = 0
        self.mLength = 0
        self.mAligned = 0        
        self.mMatches = 1
        self.mFrom1 = 0
        self.mTo1 = 0
        self.mAlignedSequence1 = ""
        self.mFrom2 = 0
        self.mTo2 = 0
        self.mAlignedSequence2 = ""        
        self.mUnaligned = 0
        self.mBlockSizes = []
        self.mGapSizes1 = []
        self.mGapSizes2 = []
        self.mDistance = 0
        self.mKappa = 0.0
        self.mAlpha = 0.0
        self.mType1 = ""
        self.mType2 = ""
        self.mIntronFrom1 = 0
        self.mIntronTo1 = 0
        self.mIntronFrom2 = 0
        self.mIntronTo2 = 0
        
    def Read( self, line ):
        """read data from tab-separated line."""

        data = string.split( line[:-1], "\t")
        
        if len(data) == 19:        
            (self.mCategory, self.mMethod,
             self.mToken1, self.mIntronId1, self.mNumIntrons1, self.mLenIntron1,
             self.mToken2, self.mIntronId2, self.mNumIntrons2, self.mLenIntron2,
             self.mNumGaps, self.mLength, self.mAligned,
             self.mFrom1, self.mTo1, self.mAlignedSequence1, 
             self.mFrom2, self.mTo2, self.mAlignedSequence2 ) = data
        elif len(data) == 12:
            (self.mCategory, self.mMethod,
             self.mToken1, self.mIntronId1, 
             self.mToken2, self.mIntronId2, 
             self.mFrom1, self.mTo1, self.mAlignedSequence1, 
             self.mFrom2, self.mTo2, self.mAlignedSequence2 ) = data
            
            (self.mNumIntrons1, self.mLenIntron1, self.mNumIntrons2, self.mLenIntron2,
             self.mNumGaps, self.mLength, self.mAligned) = [0] * 7
        elif len(data) == 23:
            (self.mCategory, self.mMethod,
             self.mToken1, self.mIntronId1, self.mNumIntrons1, self.mLenIntron1,
             self.mToken2, self.mIntronId2, self.mNumIntrons2, self.mLenIntron2,
             self.mNumGaps, self.mLength, self.mAligned,
             self.mFrom1, self.mTo1, self.mAlignedSequence1, 
             self.mFrom2, self.mTo2, self.mAlignedSequence2,
             self.mIntronFrom1, self.mIntronTo1,
             self.mIntronFrom2, self.mIntronTo2 ) = data

        ( self.mIntronId1, self.mNumIntrons1, self.mLenIntron1,
          self.mIntronId2, self.mNumIntrons2, self.mLenIntron2,
          self.mNumGaps, self.mLength, self.mAligned,
          self.mFrom1, self.mTo1, self.mFrom2, self.mTo2,
          self.mIntronFrom1, self.mIntronTo1,
          self.mIntronFrom2, self.mIntronTo2) = \
          map(int, ( self.mIntronId1, self.mNumIntrons1, self.mLenIntron1,
                     self.mIntronId2, self.mNumIntrons2, self.mLenIntron2,
                     self.mNumGaps, self.mLength, self.mAligned,
                     self.mFrom1, self.mTo1, self.mFrom2, self.mTo2,
                     self.mIntronFrom1, self.mIntronTo1,
                     self.mIntronFrom2, self.mIntronTo2, ))

    def __str__( self ):
        return string.join( map(str,\
                                ( self.mCategory, self.mMethod,
                                  self.mToken1, self.mIntronId1, self.mNumIntrons1, self.mLenIntron1,
                                  self.mToken2, self.mIntronId2, self.mNumIntrons2, self.mLenIntron2,
                                  self.mNumGaps, self.mLength, self.mAligned, 
                                  self.mFrom1, self.mTo1, self.mAlignedSequence1,
                                  self.mFrom2, self.mTo2, self.mAlignedSequence2)), "\t")


def CalculateBlockSizes( pair, gap_char = "-" ):
    """return the block-sizes in the alignment.
    """
    block_sizes = []
    left_gap_sizes = []
    right_gap_sizes = []    

    was_left_gap = pair.mAlignedSequence1[0] == gap_char
    was_block = pair.mAlignedSequence1[0] != gap_char and pair.mAlignedSequence2[0] != gap_char
    size = 1
    
    for x in range(1,min(len(pair.mAlignedSequence1), len(pair.mAlignedSequence2))):
        is_left_gap  = pair.mAlignedSequence1[x] == gap_char
        is_right_gap = pair.mAlignedSequence2[x] == gap_char
        is_block = not is_left_gap and not is_right_gap
        if is_block and not was_block:
            if was_left_gap:
                left_gap_sizes.append(size)
                size = 0
            else:
                right_gap_sizes.append(size)
                size = 0
        elif not is_block and was_block:
            block_sizes.append(size)
            size = 0
        elif is_left_gap and not was_left_gap:
            right_gap_sizes.append(size)
            size = 0
        elif is_right_gap and was_left_gap:
            left_gap_sizes.append(size)
            size = 0
        elif is_left_gap and is_right_gap:
            raise "double gap"

        was_block = is_block
        was_left_gap = is_left_gap
            
        size += 1
            
    if was_block:
        block_sizes.append(size)
    else:
        if was_left_gap:
            left_gap_sizes.append(size)
        else:
            right_gap_sizes.append(size)
        
    return block_sizes, left_gap_sizes, right_gap_sizes
    
def CalculatePercentIdentity( pair, gap_char = "-" ):
    """return number of idential and transitions/transversions substitutions
    in the alignment.
    """
    transitions   = ("AG", "GA", "CT", "TC")
    transversions = ("AT", "TA", "GT", "TG", "GC", "CG", "AC", "CA" )    

    nidentical = 0
    naligned = 0
    ndifferent = 0
    ntransitions = 0
    ntransversions = 0
    nunaligned = 0

    for x in range(min(len(pair.mAlignedSequence1), len(pair.mAlignedSequence2))):
        if pair.mAlignedSequence1[x] != gap_char and \
           pair.mAlignedSequence2[x] != gap_char:
            naligned += 1
            if pair.mAlignedSequence1[x] == pair.mAlignedSequence2[x]:
                nidentical += 1
            else:
                ndifferent += 1

            if (pair.mAlignedSequence1[x] + pair.mAlignedSequence2[x]) in transitions:
                ntransitions += 1
            if (pair.mAlignedSequence1[x] + pair.mAlignedSequence2[x]) in transversions:
                ntransversions += 1
        else:
            nunaligned += 1
            
    return nidentical, ntransitions, ntransversions, naligned, nunaligned

##------------------------------------------------------------
def GetIntronType( sequence ):
    """return intron type for an intronic sequence."""

    for name, prime5, prime3 in param_intron_types:
        if sequence[:len(prime5)].upper() == prime5 and \
           sequence[-len(prime3):].upper() == prime3:
            return name
    else:
        return "unknown-" + sequence[:5] + "-" + sequence[-5:]

##------------------------------------------------------------
def AlignPair( pair, anchor = 0 ):
    """align a pair of introns."""

    map_intron_a2b = alignlib.makeAlignmentVector()

    if param_loglevel >= 1:
        print "# aligning %s-%i with %s-%i: lengths %i and %i" % (pair.mToken1, pair.mIntronId1,
                                                                  pair.mToken2, pair.mIntronId2,
                                                                  len(pair.mAlignedSequence1),
                                                                  len(pair.mAlignedSequence2))
        sys.stdout.flush()

    s1 = "A" * anchor + pair.mAlignedSequence1 + "A" * anchor
    s2 = "A" * anchor + pair.mAlignedSequence2 + "A" * anchor

    if param_method == "dialigned":
        dialign.Align( s1, s2, map_intron_a2b )
    elif param_method == "dialignedlgs":
        dialignlgs.Align( s1, s2, map_intron_a2b ) 
    elif param_method == "dbaligned":
        dba.Align( s1, s2, map_intron_a2b )
    elif param_method == "clusaligned":
        raise NotImplementedError("clustalw wrapper not up-to-date")
        clustal.Align( s1, s2, map_intron_a2b )

    if anchor:
        map_intron_a2b.removeRowRegion( anchor + len(pair.mAlignedSequence1) + 1, map_intron_a2b.getRowTo() )
        map_intron_a2b.removeRowRegion( 1, anchor)        
        map_intron_a2b.removeColRegion( anchor + len(pair.mAlignedSequence2) + 1, map_intron_a2b.getColTo() )        
        map_intron_a2b.removeColRegion( 1, anchor)
        map_intron_a2b.moveAlignment( -anchor, -anchor )

    if map_intron_a2b.getLength() == 0:
        if param_loglevel >= 1:
            print "# Error: empty intron alignment"
        return False


    seq1 = alignlib.makeSequence( pair.mAlignedSequence1 )
    seq2 = alignlib.makeSequence( pair.mAlignedSequence2 )
    
    data = alignlib.AlignmentFormatExplicit( map_intron_a2b, seq1, seq2 )

    pair.mFrom1, pair.mAlignedSequence1, pair.mTo1 = data.mRowFrom, data.mRowAlignment, data.mRowTo
    pair.mFrom2, pair.mAlignedSequence2, pair.mTo2 = data.mColFrom, data.mColAlignment, data.mColTo
    pair.mMethod = param_method

    pair.mNumGaps, pair.mLength = map_intron_a2b.getNumGaps(), map_intron_a2b.getLength()
    pair.mAligned = pair.mLength - pair.mNumGaps

    if param_loglevel >= 2:
        print "# alignment success", pair.mAlignedSequence1, pair.mAlignedSequence2

    return True

##------------------------------------------------------------
if __name__ == '__main__':

    try:
        optlist, args = getopt.getopt(sys.argv[1:], param_short_options, param_long_options)
    except getopt.error, msg:
        print globals()["__doc__"], msg
        sys.exit(2)

    for o,a in optlist:
        if o in ( "-v", "--verbose" ):
            param_loglevel = int(a)
        elif o in ( "-h", "--help" ):
            print globals()["__doc__"]
            sys.exit(0)
        elif o in ( "--version", ):
            print "version="
            sys.exit(0)
        elif o == "--is-compressed":
            param_is_compressed = 1
        elif o == "--do-gblocks":
            param_do_gblocks = 1
        elif o == "--skip-distance":
            param_do_distance = 0
        elif o == "--skip-alistats":
            param_do_alistats = 0
        elif o == "--echo-unaligned":
            param_echo_unaligned = 1
        elif o == "--method":
            param_method = a
            param_do_alignment = True
        elif o == "--fixed-alpha":
            param_fixed_alpha = float(a)
        elif o == "--anchor-alignment":
            param_anchor_alignment = int(a)

    if len(args) > 0:
        print globals()["__doc__"], "no arguments required."
        sys.exit(2)

    baseml = WrapperBaseML.BaseML()
    if param_fixed_alpha != None:
        baseml.SetOption( "alpha", param_fixed_alpha )
        baseml.SetOption( "fix_alpha", 1 )
        
    gblocks = WrapperGblocks.Gblocks()
    
    print E.GetHeader()
    print E.GetParams()
    print baseml.GetOptions()
    sys.stdout.flush()
    
    last_token1, last_token2 = None, None
    
    nintron_pairs, ntoken_pairs = 0, 0
    ninput, nskipped, nerrors = 0, 0, 0

    if param_loglevel >= 3:
        dump_result = 1
    else:
        dump_result = 0

    unaligned_pair = None

    ## setup alignment objects
    if param_do_alignment:
        dialign = WrapperDialign.Dialign( "-n" )
        dialignlgs = WrapperDialign.Dialign( "-n -it -thr 2 -lmax 30 -smin 8" )    
        dba = WrapperDBA.DBA()    
        clustal = WrapperClustal.Clustal()

    print """# CATEGORY:       category [intron|exon]
# METHOD:         alignment method
# TOKEN:          name
# ID:             id of exon/intron
# NINTRONS:       number of introns
# LINTRON:        length of intron
# NALIGNED:       number of aligned positions
# PALIGNED:       percentage of aligned positions
# DISTANCE:       baseML distance
# IDENT:          number of identical positions
# TRANSIT:        number of transitions
# TRANSVERS:      number of transversion
# MATCHES:        number of matching positions
# PIDENT:         percentage of identical positions
# PTRANSIT:       precentage of transitions
# PTRANSVERS:     precentage of transversion
# BLOCKSIZES:     alignment, length of blocks
# GAPS:           gap sizes in sequence 1/2
# ALPHA:          estimated alpha (from BASEML)
# KAPPA:          estimated kappa (from BASEML)
# CATEGORY\tMETHOD\tTOKEN1\tID1\tNINTRONS1\tNLINTRON1\tTOKEN2\tID2\tNINTRONS2\tNLINTRON2\tNALIGNED\tPALIGNED\tDISTANCE\tIDENT\tTRANSIT\tTRANSVER\tMATCHES\tPIDENT\tPTRANSVIT\tPTRANVER\tBLOCKSIZES\tGAPSIZES\tGAPSIZES\tTYPE1\tTYPE2\tALPHA\tKAPPA"""

    for line in sys.stdin:

        if line[0] == "#": continue
        ninput += 1
        data = string.split( line[:-1], "\t")
        
        if data[0] != "intron": continue

        pair = IntronPair()
        pair.Read(line)

        if param_do_alignment:

            pair.mType1 = GetIntronType( string.replace( pair.mAlignedSequence1, "-", "" ) )
            pair.mType2 = GetIntronType( string.replace( pair.mAlignedSequence1, "-", "" ) )
            
            is_ok = AlignPair( pair, anchor = param_anchor_alignment )
            
            if not is_ok:
                nerrors += 1
                continue
            
        if len(pair.mAlignedSequence1) == 0 or len(pair.mAlignedSequence2) == 0:

            if param_loglevel >= 1:
                print "# skipped entry:", str(pair)
            nskipped += 1
            continue

        if pair.mMethod == "unaligned":
            unaligned_pair = pair
            pair.mType1 = GetIntronType( unaligned_pair.mAlignedSequence1 )
            pair.mType2 = GetIntronType( unaligned_pair.mAlignedSequence2 )
            do_print = param_echo_unaligned
        else:
            do_print = 1
            if param_is_compressed:
                if unaligned_pair and \
                       unaligned_pair.mToken1 == pair.mToken1 and \
                       unaligned_pair.mToken2 == pair.mToken2 and \
                       unaligned_pair.mIntronId1 == pair.mIntronId1:

                    map_a2b = alignlib.makeAlignmentVector()
                    f = AlignmentFormatEmissions( 
                        pair.mFrom1, 
                        pair.mAlignedSequence1,
                        pair.mFrom2, 
                        pair.mAlignedSequence2).copy( map_a2b )
                    map_a2b.moveAlignment( -unaligned_pair.mFrom1 + 1, -unaligned_pair.mFrom2 + 1 )            

                    data = alignlib.AlignmentFormatExplicit( map_a2b,
                                                             alignlib.makeSequence( unaligned_pair.mAlignedSequence1),
                                                             alignlib.makeSequence( unaligned_pair.mAlignedSequence2) )

                    from1, ali1, to1 = data.mRowFrom, data.mRowAlignment, data.mRowTo
                    from2, ali2, to2 = data.mColFrom, data.mColAlignment, data.mColTo

                    pair.mAlignedSequence1 = ali1
                    pair.mAlignedSequence2 = ali2

                else:
                    raise "sequence not found for pair %s" % str(pair)

                    
            if param_do_gblocks:
                if param_loglevel >= 4:
                    print "# length before: %i %i" % (len(pair.mAlignedSequence1), pair.mAligned)
                pair.mAlignedSequence1, pair.mAlignedSequence2 = gblocks.GetBlocks( pair.mAlignedSequence1, pair.mAlignedSequence2)
                if param_loglevel >= 4:
                    print "# length after: %i" % len(pair.mAlignedSequence1)

                if len(pair.mAlignedSequence1) == 0:
                    nerrors += 1
                    continue

            if param_loglevel >= 5:
                print ">seq1\n%s\n>seq2\n%s" % (pair.mAlignedSequence1, pair.mAlignedSequence2 )
                
            if param_do_distance:
                try:
                    result = baseml.Run( ( ("seq1", pair.mAlignedSequence1),
                                           ("seq2", pair.mAlignedSequence2)),
                                         dump_result = dump_result )
                    pair.mDistance = result.mMatrix["seq1"]["seq2"]
                    pair.mAlpha    = result.mAlpha
                    pair.mKappa    = result.mKappa
                except WrapperBaseML.ParsingError, x:
                    if param_loglevel >= 2:
                        print "# parsing error in BaseML output", x
                        
                    nerrors += 1
                    pair.mDistance = 999.0

                if pair.mDistance > 500:
                    nerrors += 1
                    print "# error", str(pair)
                    continue
            else:
                pair.mDistance = 0.0

            #if param_do_alignment:
            pair.mIdentical, pair.mTransitions, pair.mTransversions, pair.mMatches, pair.mUnaligned = \
                             CalculatePercentIdentity( pair )

            pair.mBlockSizes, pair.mGapSizes1, pair.mGapSizes2 = CalculateBlockSizes( pair )

        if do_print:
            if not pair.mLength: pair.mLength = 1
            if not pair.mMatches: pair.mMatches = 1
            print string.join( map(str, ( \
                pair.mCategory, pair.mMethod,
                pair.mToken1, pair.mIntronId1, pair.mNumIntrons1, pair.mLenIntron1,
                pair.mToken2, pair.mIntronId2, pair.mNumIntrons2, pair.mLenIntron2,
                pair.mAligned,
                "%5.4f" % (float(pair.mAligned) / float(pair.mLength)),
                pair.mDistance,
                pair.mIdentical,
                pair.mTransitions,
                pair.mTransversions,
                pair.mMatches,
                "%5.4f" % (float(pair.mIdentical) / float(pair.mMatches)),
                "%5.4f" % (float(pair.mTransitions) / float(pair.mMatches)),
                "%5.4f" % (float(pair.mTransversions) / float(pair.mMatches)),
                string.join(map(str, pair.mBlockSizes), ","),
                string.join(map(str, pair.mGapSizes1), ","),
                string.join(map(str, pair.mGapSizes2), ","),
                pair.mType1, pair.mType2,
                pair.mAlpha, pair.mKappa,
                )), "\t")
        
            nintron_pairs += 1
            
        if last_token1 != pair.mToken1 and\
           last_token2 != pair.mToken2:
            ntoken_pairs += 1
            last_token1 = pair.mToken1
            last_token2 = pair.mToken2

    print "# input=%i, skipped=%i, nerrors=%i, transcripts=%i, introns=%i" % (ninput, nskipped, nerrors,
                                                                              ntoken_pairs, nintron_pairs )
    print E.GetFooter()
