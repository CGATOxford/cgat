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
mali2mali.py - operate on a multiple alignments
===============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

perform operations (masking, renaming, conversion, ...) on a multiple alignment

Available edit operations are:

remove-unaligned-ends
   removes unaligned ends (lowercase characters)

remove-end-gaps
   removes gapped ends.

shift-alignment
   add offset to locations. This requires one parameter, which
   should be filename of a file containing id-offset in tabular form.

mark-transitions
   marks transitions in the alignment. The parameter denotes either

      1 a filename to a file containing id-transition pairs in
        tabular form. There can be multiple transitions for the same id, but
        they should be on separate lines. Use 'mali' as an id for multiple
        alignment coordinates.

      2 a colon (':') -separated list of transitions

filter-even-transitions
   same parameter as above, but will remove all even segments (0, 2, ...).

keep-odd-segments
   synonym to filter-even-transitions

filter-odd-transitions
   same parameter as above, but will remove all odd segments (1, 3, ...).

keep-even-segments
   synonym to filter-odd-transitions

lower
   convert the multiple alignment to lowercase characters.

upper
   convert the multiple alignment to uppercase characters.

remove-all-gaps
   remove columns that are completely gaps.

remove-any-gaps
   remove columns that contain at least one gap.

remove-some-gaps
   remove columns that contain some gaps. The amount is specified in
   the parameter list as minimum_gaps. A value of
   --parameter=3 will remove all columns with three or more gaps.

add-annotation
   add annotation to mali.

recount
   recount alignment ranges.

remove-empty-sequences
   removes all empty sequences from mali.

remove-unaligned-pairs
   removes unaligned pairs in the multiple alignment. 

filter-3rd
   return alignment with only third codon positions.

filter-4d               
   return alignment with only the four-fold degenerate site.

mask-random-proportion  
   mask random proportion of colums

exclude-with-stop
   remove all sequences that contain a stop codon

exclude-with-frameshift
   remove all sequences that contain an indel that is not
   a multiple of three. Indels are determined with respect to the first
   sequence.

Parameters are given to the option parameters in a comma-separated list in the order
that the edit operations are called upon.

Usage
-----

Example::

   python mali2mali.py --help

Type::

   python mali2mali.py --help

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
import optparse
import math
import time

import CGAT.Genomics as Genomics
import CGAT.Masker as Masker
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Mali as Mali

##------------------------------------------------------------
##------------------------------------------------------------
##------------------------------------------------------------
def removeUnalignedPairs( mali, options ):
    """remove unaligned pairs in multiple alignment.

    This procedure works by finding the column that has most sequences
    aligned. Only keep sequences, that are aligned in this region.
    """
    ids = mali.getIdentifiers()

    nids = len(ids)
    columns = []
    max_size = 0

    ## compute columns with maximum occupancy
    for x in range(mali.getWidth()):
        aligned = set()
        for y in range(nids):
            if mali[ids[y]][x] not in options.gap_chars:
                aligned.add( y )
        if len(aligned) > max_size:
            columns = []
            max_size = len(aligned)
        if len(aligned) == max_size:
            columns.append( (x, aligned) )

    ## compute sequences contributing most to occupancy
    counts = {}
    for col, s in columns:
        ss = list(s)
        ss.sort()
        ss = tuple(ss)
        if ss not in counts: counts[ss] = 0
        counts[ss] += 1

    ## compute sequences contributing most to occupancy
    cc = []
    for s, c in counts.items(): cc.append( (c, s) )
    cc.sort()

    take_rows = set(cc[-1][1])

    delete_rows = [ ids[x] for x in set(range(len(ids))).difference( set(take_rows) )]
    
    for id in delete_rows:
        mali.deleteEntry( id )

    E.info( "ninput=%i, nremoved=%i, noutput=%i" % (nids, len(delete_rows), mali.getLength()))

##------------------------------------------------------------
##------------------------------------------------------------
##------------------------------------------------------------
def maskColumns( mali, annotation_type, annotation_file ):
    """mask columns in mali according to annotation type."""

    if annotation_type.upper() == "SLR":
    
        try:
            import WrapperSlr
        except ImportError:
            print "# could not add Slr annotation: module WrapperSlr not found."
            return
        
        slr = WrapperSlr.Slr()
        
        result = slr.parseOutput( open(annotation_file, "r").readlines() )

        columns = []

        for site in result.mSites:
            if site.isPositive():
                columns += list( range( 3 * (site.mResidue - 1), 3 * (site.mResidue - 1) + 3) )

    mali.maskColumns( columns )

##------------------------------------------------------------
##------------------------------------------------------------
##------------------------------------------------------------
def addAnnotation(mali, annotation_type, annotation_file ):
    """add annotation to a multiple alignment."""

    if annotation_type.upper() in ("SLR-CODON", "SLR-AA"):
        
        try:
            import WrapperSlr
        except ImportError:
            print "# could not add Slr annotation: module WrapperSlr not found."
            return
        
        slr = WrapperSlr.Slr()
        
        result = slr.parseOutput( open(annotation_file, "r").readlines() )

        codes = []

        if annotation_type.upper() == "SLR-CODON":
            l = 3
        else:
            l = 1
        
        for site in result.mSites:
            if site.mResult:
                if site.mResult == "0":
                    c = "@@@"[:l]
                else:
                    c = site.mResult[:l]
            else:
                c = ""
                
            codes.append( c + "_" * (l - len(c)) )

        mali.addAnnotation( "slr", "".join(codes))

    elif annotation_type.upper() in ("CODEML-CODON", "CODEML-AA"):

        try:
            import WrapperCodeML
        except ImportError:
            print "# could not add Codeml annotation: module WrapperCodeML not found."
            return
        
        codeml = WrapperCodeML.CodeMLSites()
        
        result = codeml.parseOutput( open(annotation_file, "r").readlines() )

        codes = []

        if annotation_type.upper() == "CODEML-CODON":
            l = 3
        else:
            l = 1

        threshold = 0.95
        if options.loglevel >= 2:
            options.stdlog.write("# Warning: only positive sites of > %f significance for Bayes Empirical Bayes\n" % threshold)
            
        for x in ("1", "0", "3", "2", "7", "8"):

            sites = result.mSites[x].mBEB.mPositiveSites
            if len(sites) == 0: continue
            
            codes = ["-"] * mali.getWidth()
            found = False
            for site in sites:
                if site.mProbability < threshold: continue
                found = True
                pos = (site.mResidue - 1) * l
                codes[pos:pos+l] = ["+"] * l
            if not found: continue
            mali.addAnnotation( "codeml_%s_%s" % (x, "beb"), "".join(codes))

##------------------------------------------------------------
def maskMali( mali, method="seg" ):
    """mask multiple alignment according to an external masker.
    """
    
    if method == "seg":
        masker = Masker.MaskerSeg()
    elif method == "bias":
        masker = Masker.MaskerBias()
    elif method == "random":
        masker = Masker.MaskerRandom()
        
    if mali.getAlphabet() == "na" and method in ("seg", "bias"):
        for id, s in mali.items():
            ss = Genomics.TranslateDNA2Protein( s.mString )
            mss = masker( ss )
            columns = []
            for x in range( 0, len(mss) ):
                if mss[x] in ("X", "x"):
                    columns += range( x, x+3 )
            mali.getEntry(id).maskColumns( columns )
    else:
        for id, s in mali.items():
            mali[id].mString = masker( s.mString )

##------------------------------------------------------------
def filterMali( mali, method="3rd" ):
    """build a new multiple alignment based on a filter.

    valid methods are
    3rd:        only third positions
    4d:         only four-fold degenerate sites
    """

    if method not in ("3rd", "4d"):
        raise "unknown method %s" % method

    if method == "3rd":
        columns = range( 2, mali.getWidth(), 3 )
        
    elif method == "4d":
        ## translate
        trans_mali = Mali.Mali()
        for id, seq in mali.items():
            s = []
            sequence = seq.mString
            l = len(sequence)
            for codon in [ sequence[x:x+3] for x in range(0, l, 3) ]:
                aa = Genomics.MapCodon2AA( codon )
                s.append( aa )

            trans_mali.addSequence( id, 0, l, "".join( s ) )

        ## get four-fold (or higher) degenerate amino acids            
        aa_columns = trans_mali.getColumns()
        columns = []
        for c in range(len(aa_columns)):
            chars = set(aa_columns[c])
            chars = chars.difference(set(mali.mGapChars))
            if len(chars) == 1:
                char = list(chars)[0].upper()
                try:
                    deg = Genomics.DegeneracyAA[char]
                except KeyError:
                    continue
                if deg >= 4:
                    columns.append( c * 3 )
                    
    mali.takeColumns( columns )
        
##------------------------------------------------------------
def main( argv = sys.argv ):

    parser = E.OptionParser( version = "%prog version: $Id: mali2mali.py 2782 2009-09-10 11:40:29Z andreas $", 
                                    usage = globals()["__doc__"])

    parser.add_option("-i", "--input-format", dest="input_format", type="choice",
                      choices=("plain", "fasta", "clustal", "stockholm", "phylip" ),
                      help="input format of multiple alignment [default=%default]."  )
    
    parser.add_option("-o", "--output-format", dest="output_format", type="choice",
                      choices=("plain", "fasta", "stockholm", "phylip", "nexus", "plain-fasta" ),
                      help="output format of multiple alignment [default=%default]."  )

    parser.add_option("--with-ranges", dest="with_ranges", action="store_true",
                      help="output alignment ranges (suffix /from-to after identifier) [default=%default]."  )

    parser.add_option( "--without-ranges", dest="with_ranges", action="store_false",
                      help="do not output alignment ranges (suffix /from-to after identifier) [default=%default]."  )

    parser.add_option( "-u", "--allow-duplicates", dest="allow_duplicates", action="store_true",
                      help="permit duplicate entries [default=%default]."  )

    parser.add_option("-m", "--method", dest="methods", type="string",
                      help="""methods to apply. Several methods can be specified in a ','-separated list [default=%default]."""  )

    parser.add_option("-p", "--parameters", dest="parameters", type="string",
                      help="parameter stack for methods that require one [default=%default]."  )

    parser.add_option("-a", "--mask-char", dest="mask_char", type="string",
                      help="character to identify/set masked characters [default=%default]."  )

    parser.set_defaults(
        input_format="fasta",
        output_format="fasta",
        methods = "",
        parameters = "",
        mask_char = "x",
        gap_chars = "-.nN",
        with_ranges = True,
        allow_duplicates = False,
        )

    (options, args) = E.Start( parser )

    options.methods = options.methods.split(",")
    options.parameters = options.parameters.split(",")    
    
    ## 1. read multiple alignment in various formats
    if options.allow_duplicates:
        mali = Mali.SequenceCollection()
    else:
        mali = Mali.Mali()

    t1 = time.time()
    
    mali.readFromFile( options.stdin, format = options.input_format )

    E.info( "read mali with %i entries in %i seconds." % (len(mali), time.time() - t1))
        
    if len(mali) == 0:
        raise ValueError("empty multiple alignment")

    for method in options.methods:

        t1 = time.time()
        
        if method == "remove-unaligned-ends":
            mali.removeUnalignedEnds()
        elif method == "remove-end-gaps":
            mali.removeEndGaps()
        elif method == "remove-all-gaps":
            mali.removeGaps( minimum_gaps = len(mali) )
        elif method == "remove-any-gaps":
            mali.removeGaps( minimum_gaps = 1 )
        elif method == "remove-some-gaps":
            minimum_gaps = int(options.parameters[0])
            del options.parameters[0]
            mali.removeGaps( minimum_gaps = minimum_gaps )
        elif method == "remove-empty-sequences":
            mali.removeEmptySequences()
        elif method == "upper":
            mali.upperCase()            
        elif method == "lower":
            mali.lowerCase()
        elif method == "mark-codons":
            mali.markCodons()
        elif method == "remove-stops":
            mali.removePattern( lambda x: x.upper() in ("TAG", "TAA", "TGA"),
                                allowed_matches = 0,
                                minimum_matches = 1,
                                delete_frame = 3,
                                search_frame = 3)
        elif method == "shift-alignment":
            map_id2offset = IOTools.ReadMap( open( options.parameters[0], "r"),
                                             map_functions=(str,int) )
            del options.parameters[0]
            mali.shiftAlignment( map_id2offset )
        elif method == "propagate-masks":
            mali.propagateMasks( mask_char = options.mask_char )

        elif method == "recount":
            mali.recount()

        elif method in ("mark-transitions", "filter-odd-transitions", "filter-even-transitions",
                        "keep-even-segments", "keep-odd-segments" ):

            if os.path.exists(options.parameters[0]):
                map_id2transitions = IOTools.readMultiMap( open( options.parameters[0], "r"),
                                                           map_functions=(str,int) )
            else:
                map_id2transitions={}
                r = map(int, options.parameters[0].split(':'))
                r.sort()
                map_id2transitions["mali"] = r

            del options.parameters[0]
            if method == "mark-transitions":
                mali.markTransitions( map_id2transitions )            
            elif method in ("filter-odd-transitions", "keep-even-segments"):
                mali.markTransitions( map_id2transitions, mode = "keep-odd" )
            elif method in ("filter-even-transitions", "keep-odd-segments"):
                mali.markTransitions( map_id2transitions, mode = "keep-even" )                                            

        elif method == "propagate-transitions":
            mali.propagateTransitions()
                
        elif method == "map-annotation":
            ## map annotations in one mali (stockholm-format) to the annotations in another.
            ## Note: the first two sequence identifiers must be shared and the sequence of the
            ## same length
            other_mali = Mali.Mali()
            other_mali.readFromFile( open(options.parameters[0], "r"), format="stockholm")
            del options.parameters[0]            
            mali.copyAnnotations( other_mali )

        elif method == "add-annotation":
            annotation_type, annotation_file = options.parameters[:2]
            del options.parameters[:2]
            AddAnnotation( mali, annotation_type, annotation_file )

        elif method == "mask-columns":
            annotation_type, annotation_file = options.parameters[:2]
            del options.parameters[:2]
            maskColumns( mali, annotation_type, annotation_file )

        elif method == "remove-unaligned-pairs":
            removeUnalignedPairs( mali, options )

        elif method == "filter-3rd":
            filterMali( mali, "3rd" )

        elif method == "filter-4d":
            filterMali( mali, "4d" )

        elif method in ("mask-seg", "mask-bias" ):
            a, b = method.split("-")
            maskMali( mali, b )
            
        elif method == "exclude-with-stop":
            mali.filter( method = "with-stop" )

        elif method == "exclude-with-stop":
            mali.filter( method = "with-frameshift" )

        E.info( "applied method %s in %i seconds." % (method, time.time() - t1))

            
    mali.writeToFile( options.stdout, 
                      format = options.output_format, 
                      write_ranges = options.with_ranges )
        
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
    
