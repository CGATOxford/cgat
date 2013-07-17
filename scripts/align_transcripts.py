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
align_transcripts.py - multiply align transcripts
=================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script will multiply align transcripts. Transcripts
are aligned as amino acid sequences and later back-translated
to codon sequences.

The script is aware of gene structures will correctly align
exons in sequence. 

Usage
-----

Example::

   python align_transcripts.py --help

Type::

   python align_transcripts.py --help

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

USAGE="""python %s [OPTIONS] < transcripts > filtered

Align peptide/DNA sequences.

If there are multiple transcripts for a gene, concatenate the
exons into one pseudo-transcript and submit that to the multiple
alignment in order to avoid incompatible exons to be aligned.
After alignment, the pseudo-transcripts are split into separate
parts again.

If the peptide sequence contains at least one character ``U``, the
sequence is assumed to be selonoprotein. 

""" % sys.argv[0]

import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.Mali as Mali
import CGAT.Exons as Exons
import CGAT.WrapperMuscle as WrapperMuscle
import CGAT.Intervals as Intervals
from peptides2cds import getMapPeptide2Cds
import alignlib
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools

def buildGeneMap( identifiers, separator = "|" ):
    """build map of predictions to genes.

    Use an identifier syntax of species|transcript|gene.
    If none is given, all transcripts are assumed to be 
    from their own gene.
    """

    map_id2gene, map_gene2ids = {}, {}

    for id in identifiers:
        f = id.split(separator)
        if len(f) < 3:
            gene = id
        else:
            gene = f[0] + separator + f[2]
        map_id2gene[id] = gene
        if gene not in map_gene2ids: map_gene2ids[gene] = []
        map_gene2ids[gene].append(id)

    return map_id2gene, map_gene2ids

class Segment:
    def __init__(self):
        pass

#---------------------------------------------------------
#---------------------------------------------------------
#---------------------------------------------------------
## Note: maps are in one-based coordinates, but residue
## return values are in zero-based coordinates
## The calling functions have to deal with this
def MapC2PRight( map_p2c, residue ):
    """mapping function."""
    if residue == map_p2c.getColTo():
        return map_p2c.getRowTo()

    while residue < map_p2c.getColTo():
        x = map_p2c.mapColToRow( residue )
        if x >= 0: return x
        residue += 1

    return -1

def MapP2CRight( map_p2c, residue ):
    """mapping function."""

    if residue == map_p2c.getRowTo():
        return map_p2c.getColTo()

    while residue < map_p2c.getRowTo():
        x = map_p2c.mapRowToCol( residue )
        if x >= 0: return x
        residue += 1

    return -1

def MapP2CLeft( map_p2c, residue ):
    """mapping function."""
    while residue >= 0:
        x = map_p2c.mapRowToCol( residue )
        if x >= 0: return x
        residue -= 1

    return -1

#---------------------------------------------------------
#---------------------------------------------------------
#---------------------------------------------------------

def buildFragments( exons, input, pseudogenes, options, coordinate_factor=3,
                    map_peptide2cds = {},
                    cds_sequences = {},
                    ):
    """build fragments out of a sorted list of overlapping exons.

    The exon list has to be sorted by genome_from, length.

    The algorithm works in the following way::
    
       start                                                                   |end
       |
       =======================================================================
          ====================================================================
              =============================================
          ----                                             -------------------
       ---    ---------------------------------------------

    There are three coordinates for each exon involved:

    -> genomic exon coordinates: fields genome_from, genome_to
    -> coordinates in peptide sequence: fields exon_from, exon_to
    -> coordinates in cds sequence: correspond to exon_from, exon_to,
       unless the transcript contains frameshifts. cds coordinates
       are obtained using the maps in map_peptide2cds.

    If map_peptide2cds is given, it is used for mapping exon boundaries.
    If cds_sequences is given, nucleotide alignments are produces as
    well as aa alignments.
    """


    # intervals
    exon_intervals = map(lambda x: (x.mGenomeFrom, x.mGenomeTo), exons)

    # build intersection set of intervals
    small_intervals = Intervals.getIntersections( exon_intervals )

    nomap = lambda x: x

    if options.loglevel >= 5:
        for e in exons:
            options.stdlog.write( "# exon: %s\n" % str(e))

    if coordinate_factor != 1:

        frames = []
        # collect frames for each interval and assign exons to it
        # both outgoing and in-going frames are checked
        # Frame is the number of residues at start of exon that
        # belong to previous codon.
        
        for start, end in small_intervals:
            ff = {}
            for e in exons:
                if min(end, e.mGenomeTo) - max(start, e.mGenomeFrom) <= 0:
                    continue

                ## collect the three sets of coordinates for this exon
                genome_from, genome_to = e.mGenomeFrom, e.mGenomeTo
                peptide_from, peptide_to = e.mPeptideFrom, e.mPeptideTo
                try:
                    cds_from, cds_to = e.mCdsFrom, e.mCdsTo
                except AttributeError:
                    continue
                    raise "no cds boundaries for %s" % (str(e))
                
                # get cds coordinates of this fragment within exon
                cds_start = cds_from + (start - genome_from)
                cds_end   = cds_to - (genome_to - end) 

                # get map to adjust frame-shifted coordinates
                if e.mQueryToken in map_peptide2cds:
                    map_p2c = lambda x: MapP2CRight( map_peptide2cds[e.mQueryToken], x)
                    map_c2p = lambda x: MapC2PRight( map_peptide2cds[e.mQueryToken], x)                    
                else:
                    map_p2c = nomap
                    map_c2p = nomap
                
                # get peptide coordinates of this fragment within exon
                peptide_start = map_c2p(cds_start)
                peptide_end = map_c2p(cds_end) 
                
                # compute frame based on peptide coordinates
                if peptide_start % coordinate_factor != 0:
                    start_frame = coordinate_factor - peptide_start % coordinate_factor
                else:
                    start_frame = 0

                end_frame = peptide_end % coordinate_factor

                # adjust frame: frame of fragment - original frame
                key = (start_frame, end_frame)
                if key not in ff: ff[key] = []
                ff[key].append(e)
                
            frames.append(ff)

        segments = []
        
        ## build up segments according to exons
        ## iteration over intervals
        for x in range(len(small_intervals)):

            start,end = small_intervals[x]
            l = end - start

            ## iteration over frames within interval
            for key, exons_in_intervall_and_frame in frames[x].items():

                start_frame, end_frame = key
                
                if options.loglevel >= 4:
                    options.stdlog.write("# interval: %i-%i in frame %i:%i\n" % (start, end, start_frame,end_frame) )
                
                segment = Segment()
                segment.mGenomeFrom = start
                segment.mGenomeTo = end
                segment.frame = start_frame
                segment.mMembers = []
                segment.mFragments = []
                segment.mCds = []

                ## iteration over exons with the same frame within an interval
                for e in exons_in_intervall_and_frame:

                    # sanity check - why was this added?
                    if min( end, e.mGenomeTo) - max(start, e.mGenomeFrom) <= 0:
                        raise "within sanity check."
                        continue
                        
                    assert( e.mQueryToken not in segment.mMembers )

                    ## collect the three sets of coordinates for this exon
                    genome_from, genome_to = e.mGenomeFrom, e.mGenomeTo
                    peptide_from, peptide_to = e.mPeptideFrom, e.mPeptideTo
                    cds_from, cds_to = e.mCdsFrom, e.mCdsTo

                    # get cds coordinates of this fragment within exon
                    cds_start = cds_from + (start - genome_from)
                    cds_end   = cds_to - (genome_to - end)

                    # get map to adjust frame-shifted coordinates
                    if e.mQueryToken in map_peptide2cds:
                        map_p2c = lambda x: MapP2CRight( map_peptide2cds[e.mQueryToken], x)
                        map_c2p = lambda x: MapC2PRight( map_peptide2cds[e.mQueryToken], x)                   
                    else:
                        map_p2c = nomap
                        map_c2p = nomap

                    # get peptide coordinates of this fragment within exon
                    peptide_start = map_c2p(cds_start)
                    peptide_end = map_c2p(cds_end)

                    # compute frame based on peptide coordinates
                    #
                    # if frame is not 0, include residue from previous exon. 
                    # Thus: split codon residues are always assigned to the next exon.
                    #
                    # IAW: Frame is the number of residues at start of exon that
                    # belong to previous codon.

                    if start_frame:
                        dfrom = 3 - start_frame
                    else:
                        dfrom, dto = 0, 0

                    dto = (peptide_end - peptide_start + dfrom) % 3

                    x = int( (peptide_start - dfrom) / coordinate_factor)
                    y = int( (peptide_end - dto) / coordinate_factor)

                    if options.loglevel >= 5:
                        options.stdlog.write("# %s: start_frame=%i end_frame=%i dfrom=%i dto=%i x=%i y=%i peptide:%4i-%4i cds:%4i-%4i len=%i len=%i seq=%s\n" % \
                                                 (e.mQueryToken, start_frame, end_frame,
                                                  dfrom, dto, x, y,
                                                  peptide_start, peptide_end,
                                                  cds_start, cds_end,                                                  
                                                  len(input[e.mQueryToken]),
                                                  y-x, 
                                                  input[e.mQueryToken][x:y]))

                    ## assertion only work for genes without frameshifts
                    assert( (peptide_start - dfrom) % coordinate_factor == 0)
                    if (peptide_start - dfrom) % coordinate_factor != 0:
                        pseudogenes.add( e.mQueryToken )

                    assert( (peptide_end - dto) % coordinate_factor == 0)

                    if (peptide_end - dto) % coordinate_factor != 0 :
                        pseudogenes.add( e.mQueryToken )

                    # assert that boundaries work, but ignore errors in last exon or intervals that are smaller than a codon
                    if x < len(input[e.mQueryToken]) and end - start > 3:
                        assert x <= y, "sanity check x <= y failed: x=%i, y=%i, l=%i" % (x, y, len(input[e.mQueryToken]) )

                    # only add meaningfull fragments with non-zero length.
                    if x < y and x < len(input[e.mQueryToken]) and x >= 0:
                        segment.mMembers.append( e.mQueryToken )
                        pep_fragment = input[e.mQueryToken][x:y] 
                        segment.mFragments.append( pep_fragment )

                        ## add nucleotides for the peptide sequence whilst
                        ## skipping over codons containing frameshifts.
                        ## adding "NNN" for skipped codons (identified by
                        ## gaps that are multiples of 3)
                        ## walk along the peptide sequence
                        if cds_sequences and e.mQueryToken in cds_sequences:

                            my_map = map_peptide2cds[e.mQueryToken]
                            # keep case of sequence
                            cds_sequence = cds_sequences[e.mQueryToken]
                            s = []
                            last_i = x * coordinate_factor
                            last_m = my_map.mapRowToCol(x)
                            
                            for peptide_pos in range(x, y):
                                
                                pos = [ my_map.mapRowToCol( i ) for i in range( peptide_pos * coordinate_factor,
                                                                                (peptide_pos + 1) * coordinate_factor) ]
                                
                                ## do not add codons that are incomplete in the peptide sequence
                                if -1 in pos: 
                                    s.append( "NNN" )
                                    continue
                                
                                s += [cds_sequence[x] for x in pos ]

#                             for i in range( y * coordinate_factor, my_map.getRowTo()):
#                                 m = my_map.mapRowToCol(i)
#                                 if m >= 0:
#                                     break

#                             if i - last_i >= 3:
#                                 s.append( "NNN" * int(math.floor( (i-last_i) / 3.0)))
                                
#                             # fill up so that cds sequence is in-frame with nucleotide sequence
                            cds_fragment = "".join(s)

#                             if len(s) % 3 != 0:
#                                  s += "N" * (3 - (len(s) % 3) )

                            segment.mCds.append( cds_fragment )                                

                assert( len(segment.mCds) == len(segment.mFragments) )
                
                ## sanity check: are all sequence segments the same?
                # do not compare boundary residues, as they might be different
                # due to alternative splicing.
                xset = set()
                for x in segment.mFragments:
                    xset.add(x[1:-1])

                ## sequence segments might not be identical
                ## due to different locations of frameshifts within the
                ## sequences. In this case resolve by aligning
                ## the segments and providing the consensus alignment as sequence
                if len(xset) > 1:

                    muscle = WrapperMuscle.Muscle()
                    mali = Mali.Mali()
                    for x in range(len(segment.mFragments)):
                        mali.addSequence( segment.mMembers[x], 0, 0, segment.mFragments[x] )

                    aligned = muscle.Run( mali )

                    pep_consensus = aligned.getConsensus()
                    if options.loglevel >= 6:
                        options.stdlog.write("# consensus peptide alignment:\n")
                        aligned.writeToFile( options.stdlog )
                        options.stdlog.write( pep_consensus + "\n" )
                        
                    ## substitute for consensus
                    segment.mMali = aligned
                    segment.mFragments = [pep_consensus]

                    ## do the same for the nucleotide fragments
                    ## thread each sequence through the peptide alignment
                    if segment.mCds:
                        
##                         muscle = WrapperMuscle.Muscle()
##                         mali = Mali.Mali()
##                         for x in range(len(segment.mCds)):
##                             mali.addSequence( segment.mMembers[x], 0, 0, segment.mCds[x] )
                            
##                         aligned = muscle.Run( mali )

##                         cds_consensus = aligned.getConsensus( mark_with_gaps = True )
                        
##                         if options.loglevel >= 6:
##                             options.stdlog.write("# consensus cds alignment:\n")
##                             aligned.writeToFile( options.stdlog )
##                             options.stdlog.write( cds_consensus + "\n" )
##                             options.stdlog.flush()

                        cds_consensus = [ "-" ] * len(pep_consensus) * 3

                        ## mask inconsistent positions in the consensus peptide
                        ## these are due to frameshifts within an exon
                        columns = segment.mMali.getColumns()
                        for c in range(len(columns)):
                            s = columns[c]
                            counts = [(a, s.count(a)) for a in set(list(s)).difference( set("-") ) ]
                            if len(counts) > 1:
                                cds_consensus[c*3:c*3+3] = ["N"] * 3

                        for c in range(len(segment.mCds)):
                        
                            if options.loglevel >=2:
                                options.stdlog.write("# building map between consensus peptide and %s.\n" % segment.mMembers[c] )
                                options.stdlog.flush()
                            
                            # build map of consensus cds to peptide cds
                            cons_map_p2c = alignlib.makeAlignmentVector()
                            this_cds = segment.mCds[c]
                            try:
                                cons_map_p2c = getMapPeptide2Cds( pep_consensus,
                                                                  this_cds,
                                                                  options )

                            except ValueError, msg:
                                if options.loglevel >= 2:
                                    options.stdlog.write("# Warning: sequence %s not mappable: %s\n" % (segment.mMembers[c],msg) )
                                continue
                            
                            for x in range( cons_map_p2c.getRowFrom(), cons_map_p2c.getRowTo()):
                                y = cons_map_p2c.mapRowToCol( x )
                                if y < 0: continue
                                if cds_consensus[x] == "-":
                                    cds_consensus[x] = this_cds[y]

                        cds_consensus = "".join(cds_consensus)
                        if options.loglevel >= 6:
                            options.stdlog.write("# consensus cds alignment %i, %i:\n" % (len(cds_consensus), len(pep_consensus) ))
                            options.stdlog.write( cds_consensus + "\n" )
                            options.stdlog.flush()
                                
                        ## substitute for consensus
                        segment.mCds = [cds_consensus]

                else:
                    segment.mMali = None
                    
                if len(segment.mFragments) == 0:
                    continue
                
                segment.mSequence = segment.mFragments[0]
                segment.mCdsSequence = segment.mCds[0]

                # assert(len(segment.mSequence) == len(segment.mCdsSequence) / 3 )
                
                segments.append( segment)
                
    return segments

##------------------------------------------------------------
def writeToFile( mali, section, options, is_aligned = True ):
    """write mali to file."""

    outfile = open(options.output_filename_pattern % section, "w" )

    mali.writeToFile(outfile, format=options.output_format )
        
    outfile.close()
        
    if is_aligned and not mali.checkLength():
        raise "mali in file %s has entries of different lengths" % (options.output_filename_pattern % section)
        
##------------------------------------------------------------
if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: align_transcripts.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-m", "--master", dest="master", type="string",
                      help="master sequence."  )

    parser.add_option("-p", "--master-pattern", dest="master_pattern", type="string",
                      help="master pattern."  )

    parser.add_option( "--master-species", dest="master_species", type="string",
                      help="species to use as master sequences."  )

    parser.add_option("-t", "--translate", dest="filename_translation", type="string",
                      help="filename on where to store translated sequences."  )

    parser.add_option("-e", "--exons", dest="filename_exons", type="string",
                      help="filename on where to exon information."  )

    parser.add_option("-g", "--gtf", dest="filename_gtf", type="string",
                      help="filename with exon information in gtf format."  )

    parser.add_option("-c", "--mark-codons", dest="mark_codons", action="store_true",
                      help="mark codons."  )

    parser.add_option( "--remove-stops", dest="remove_stops", action="store_true",
                      help="remove stop codons."  )

    parser.add_option("--mask-stops", dest="mask_stops", action="store_true",
                      help="mask stop codons."  )

    parser.add_option("--mask-char", dest="mask_char", type="string",
                      help="masking character to use."  )

    parser.add_option("-f", "--remove-frameshifts", dest="remove_frameshifts", action="store_true",
                      help="remove columns corresponding to frameshifts."  )

    parser.add_option("-s", "--split-exons", dest="split_exons", action="store_true",
                      help="split columns aligned to different exons in the same gene."  )

    parser.add_option("-a", "--target", dest="target", type="choice",
                      choices=("paml", ),
                      help="perform cleaning up for certain targets." )
    
    parser.add_option("--stop-at", dest="stop_at", type="choice",
                      choices=("aligned", "unaligned", "unpacked"),
                      help="stop at intermediate stage and dump out that alignment.")

    parser.add_option("--force-map", dest="force_map", action="store_true",
                      help="force mapping of sequences that have changed to previous sequence." )

    parser.add_option( "--cds", dest="filename_cds", type="string",
                      help="""filename with cds - useful if you expect pseudogenes in your set. The peptide
sequences will be aligned to the cds sequences. This produces better coordinates."""  )

    parser.add_option("--output", dest="output", type="choice", action="append",
                      choices=("final_aa", "final_na", "aligned_aa", "aligned_na", "all", "unaligned_aa", "unaligned_na", "coords" ),
                      help="which alignment to output: aligned=aligned sequences, but before untangling exons; "
                           " final=final multiple alignment; unaligned=unaligned sequences; "
                           " coords=genomic coordinates (corresponding to 'final_na')." )
    
    parser.add_option("--output-filename-pattern", dest="output_filename_pattern", type="string",
                      help="filename pattern for multiple alignment output files. "
                           " If no --output option is given, stdout is used." )

    parser.add_option("--output-filename-coords", dest="output_filename_coords", type="string",
                      help="filename to output coordinates to." )

    parser.add_option("--output-format", dest="output_format", type="choice",
                      choices=("fasta", "stockholm", "clustal", "plain-fasta") ,
                      help="output format of multiple alignments." )

    parser.add_option("--strict", dest="strict", action="store_true",
                      help="apply strict checking." )
    
    parser.set_defaults(
        gap_char = "-",
        mask_char = "x",
        gap_chars = "-.",
        separator = "|",
        master = None,
        master_species = None,
        filename_translation = None,
        filename_exons = None,
        filename_gtf = None,
        master_pattern = None,
        remove_stops = False,
        mark_codons = False,
        mask_unaligned = False,
        split_exons = False,
        remove_frameshifts = False,
        min_segment_length = 5,
        mask_stops = False,
        target = None,
        sequence_separator = "",
        stop_at = None,
        filename_cds = None,
        force_map = False,
        output_filename_pattern="%s.fasta",
        output_filename_coords="coords.tsv.gz",
        output = [],
        output_format = "fasta",
        max_percent_gaps = 0.1,
        max_gaps = 40,
        strict= False,
        )

    (options, args) = E.Start( parser )

    ########################################################
    if "all" in options.output:
        options.output = ["final_aa", "final_na", "aligned_aa", "aligned_na", "unaligned_aa", "unaligned_na", "coords"]

    ########################################################
    ########################################################
    ########################################################    
    # read unaligned sequences
    input = Mali.Mali()
    input.readFromFile( sys.stdin, format="fasta" )

    all_identifiers = input.getIdentifiers()

    alphabet = input.getAlphabet()

    E.info( "sequences are of alphabet: %s" % alphabet )
        
    if alphabet == "aa":
        coordinate_factor = 3.0

        # build list of selenoproteins
        selenoproteins = set( [x.mId for x in input.values() if "U" in x.mString.upper() ] )

    else:
        coordinate_factor = 1.0
        selenoproteins = ()

    ########################################################
    ########################################################
    ########################################################    
    # read cds sequences
    if options.filename_cds:
        cds_sequences = Genomics.ReadPeptideSequences( open(options.filename_cds, "r") )
        map_peptide2cds = {}

        nskipped = 0
        for x in all_identifiers:

            E.debug( "building map_peptide2cds for sequence %s" % (x) )
                
            if x not in cds_sequences:
                nskipped += 1
                continue
            
            try:
                map_p2c = getMapPeptide2Cds( input[x],
                                             cds_sequences[x],
                                             options )

            except ValueError, msg:
                E.warn( "sequence %s not mappable: %s" % (x,msg) )
                nskipped += 1
                continue

            num_peptide_gaps = len( re.sub("[^-]", "", input[x] ) )
            ngaps = map_p2c.getNumGaps() - (num_peptide_gaps *  3) - abs(len(input[x]*3)-len(cds_sequences[x]))
            
            l = map_p2c.getLength()
            if float(ngaps) / float(l) > options.max_percent_gaps or ngaps > options.max_gaps:
                E.warn( "map between cds to peptide has too many gaps: %s: %i gaps out of length %i" % (x, ngaps, l) )
                nskipped += 1
                continue
            
            map_peptide2cds[x] = map_p2c

            if options.loglevel >= 3:
                f = alignlib.AlignmentFormatEmissions( map_p2c )
                options.stdlog.write("# p2c: " + "\t".join( map(str, (x, str(f),
                                                                      len(input[x]), len(cds_sequences[x])) ) )+"\n")
            

        E.info("built %i maps between peptides and cds - %i skipped" % \
                   (len(map_peptide2cds),
                    nskipped ))
    else:
        map_peptide2cds = {}

    map_id2gene, map_gene2ids = buildGeneMap( all_identifiers )

    ########################################################
    ########################################################
    ########################################################
    ## read exon informmation
    if options.filename_exons or options.filename_gtf:

        # read exon boundaries and keep forward coordinates
        if options.filename_exons:
            exons = Exons.ReadExonBoundaries( open(options.filename_exons, "r"),
                                              filter=set(all_identifiers),
                                              from_zero = True)
        elif options.filename_gtf:
            exons = Exons.ReadExonBoundaries( open(options.filename_gtf, "r"),
                                              filter=set(all_identifiers),
                                              format = "gtf",
                                              gtf_extract_id = re.compile('transcript_id \"(\S+)\";'),
                                              from_zero = True)

            
        ########################################################
        ########################################################        
        ########################################################
        ## Exon boundaries contain the fields mPeptideFrom and mPeptideTo.
        ## These fields correspond to the cds coordinates, not peptide coordinates.
        ## This is no problem, if there are no frameshifts. However, in the presence
        ## of frameshifts, these two coordinates will be different. Thus, if cds
        ## sequences are given, the exon boundaries are corrected, such that
        ## mPeptideFrom and mPeptideTo contain the peptide coordinates and
        ## mCdsFrom and mCdsTo contain the cds coordinates.
        ########################################################
        ########################################################
        ########################################################
        
        ## adjust exon boundaries
        if map_peptide2cds:
            
            # if map_peptide2cds is given, pseudogenes can be treated properly
            E.info("checking exon boundaries." )
                
            nmissing, ndifferences, nstop_codons, ndeleted_empty, nunmappable = 0, 0, 0, 0, 0
            
            # minimum genomic coordinates and strands for a gene
            genome_starts = {}
            # maps of cds sequence to genomic coordinates. These are
            # zeroed and increasing for both forward and reverse strand
            map_cds2genome = {}
            
            for key, ee in exons.items():

                if key not in map_peptide2cds:
                    nmissing += 1
                    continue

                map_p2c = map_peptide2cds[key]

                ## a patch to eliminate empty last exons
                if ee[-1].mPeptideTo == ee[-1].mPeptideFrom and \
                   ee[-1].mGenomeTo == ee[-1].mGenomeFrom:
                    E.warn("%s of length %i: deleting empty last exon: %s." % \
                               (key, len(input[key]), str(ee[-1])))
                    del ee[-1]
                    ndeleted_empty += 1
                    
                if ee[-1].mPeptideTo != map_p2c.getColTo():

                    E.debug( "%s" % str(ee[-1]) )
                    E.warn( "%s of length %i: peptide and exon do not correspond: %i <> %i" %\
                                (key, len(input[key]), ee[-1].mPeptideTo, map_peptide2cds[key].getColTo()) )
                    ndifferences += 1
                    d = ee[-1].mPeptideTo - map_peptide2cds[key].getColTo()
                    if d == 3:
                        E.warn("%s: assuming difference is stop codon - exon shortened." )
                        nstop_codons += 1
                    elif d > 0:
                        E.warn("%s: fixing difference of %i nucleotides - incomplete stop-codon?" % (key,d) )
                    else:
                        # if the exon information extends to higher residues than the map_peptide2cds
                        raise ValueError("not able to fix exon boundary incongruence: %i != %i." %\
                                             (ee[-1].mPeptideTo, map_p2c.getColTo()))
                    
                    old_peptide_end, old_genome_end = ee[-1].mPeptideTo, ee[-1].mGenomeTo 
                    ee[-1].mPeptideTo -= d
                    ee[-1].mGenomeTo -= d

                    d = ee[-1].mPeptideTo - ee[-1].mPeptideFrom
                    if d <= 0:
                        del ee[-1]
                        ee[-1].mPeptideTo += d
                        ee[-1].mGenomeTo += d

                    E.debug( "%s: fixed exon end from %i to %i (%i to %i)" % \
                                 (key, old_peptide_end, ee[-1].mPeptideTo,
                                  old_genome_end, ee[-1].mGenomeTo ) )

                is_negative_strand = ee[0].mSbjctStrand == "-"
                # note that exon coordinates are already inverted
                # and negative strand coordinates are negative
                # Thus the following works for both forward and reverse strand
                genome_start = ee[0].mGenomeFrom
                genome_starts[key] = (is_negative_strand, genome_start)

                map_c2g = alignlib.makeAlignmentBlocks()

                for e in ee:
                    # map boundaries
                    # note: map_p2c is in 1 based coordinates
                    peptide_from = MapC2PRight(map_p2c, e.mPeptideFrom )
                    peptide_to = MapC2PRight(map_p2c, e.mPeptideTo)

                    if peptide_from < 0 or peptide_to < 0:
                        E.debug( "%s" % str(e) )
                        E.warn( "%s of length %i: exon boundary could not be mapped: from %i->%i to %i->%i" %\
                                    (key, len(input[key]),
                                     e.mPeptideFrom, peptide_from,
                                     e.mPeptideTo, peptide_to))
                        E.debug("%s" %str(alignlib.AlignmentFormatEmissions( map_p2c )) )
                        nunmappable += 1
                        continue

                    e.mCdsFrom = e.mPeptideFrom
                    e.mCdsTo = e.mPeptideTo
                    e.mPeptideFrom = peptide_from
                    e.mPeptideTo = peptide_to

                    # build map of cds to genomic sequence
                    map_c2g.addDiagonal( e.mCdsFrom, e.mCdsTo, e.mGenomeFrom - genome_start - e.mCdsFrom ) 
                    map_cds2genome[key] = map_c2g

            E.info("checked exon boundaries against cds: missing=%i, differences=%i, fixed_stops=%i, deleted_empty=%i, nunmappable=%i" %\
                                     (nmissing, ndifferences, nstop_codons, ndeleted_empty, nunmappable))

        else:
            E.info("no checking of exon boundaries - assumed to be correct.")
            E.info("removing stop codons.")
            #
            # remove stop codons from exons
            #
            # Tests whether the length of the peptide sequence and
            # the exon correspond. If there is a difference, truncate.
            # Note: this fails for pseudogenes. Thus use peptide2cds
            # to map these
            for key, ee in exons.items():

                reference_length = 3 * len(input[key])

                if ee[-1].mPeptideTo > reference_length:
                    ee[-1].mPeptideTo -= 3
                    ee[-1].mGenomeTo -= 3

                if ee[-1].mPeptideTo - ee[-1].mPeptideFrom == 0:
                    del ee[-1]

                for e in ee:
                    e.mCdsFrom = e.mPeptideFrom
                    e.mCdsTo = e.mPeptideTo

            map_peptide2cds = {}

        E.debug( "read exons for %i sequences." % len(exons) )

    else:
        # if no exons are given, all transcripts are assumed to be single exon
        # genes
        exons = {}
        for gene, identifiers in map_gene2ids.items():
            for id in identifiers:
                e = Exons.Exon()
                e.mQueryToken = id
                e.mPeptideFrom = map_peptide2cds[id].getRowFrom()
                e.mPeptideTo = map_peptide2cds[id].getRowTo()
                e.mGenomeFrom = map_peptide2cds[id].getColFrom()
                e.mGenomeTo = map_peptide2cds[id].getColTo()
                e.mCdsFrom = e.mGenomeFrom
                e.mCdsTo = e.mGenomeTo
                exons[id] = [e]
                
    ##########################################################################
    ##########################################################################
    ##########################################################################
    ## Build the packed alignment
    ##########################################################################        
    unaligned = Mali.Mali()
    unaligned_cds = Mali.Mali()
    
    map_gene2fragments = {}

    pseudogenes = set()
    nwarnings_length, nwarnings_sequence = 0, 0
    
    nwarnings_empty = 0
    
    for gene, identifiers in map_gene2ids.items():
        
        ## collect all exons and sort them by genomic location
        exon_list = []
        is_seleno = False
        for id in identifiers:
            exon_list += exons[id]
            is_seleno |= id in selenoproteins

        # sort such that the larger exons come first
        # if they start at the same residue
        exon_list.sort( lambda x,y: cmp( (x.mGenomeFrom, -x.mGenomeTo), (y.mGenomeFrom, -y.mGenomeTo)) )

        e = exon_list[0]

        max_to = e.mGenomeTo
        overlapping_exons = [ e ]
        fragments = []
        for e in exon_list[1:]:

            # ignore empty exons
            # Example for empty exon: ENSP00000344726, where the last exon
            # just consists of the stop codon. Should be filtered out, but
            # you never know.
            if e.mPeptideFrom == e.mPeptideTo: continue

            if max_to <= e.mGenomeFrom:
                # no overlap, process chunks
                fragments += buildFragments( overlapping_exons, input, 
                                             pseudogenes, options,
                                             map_peptide2cds = map_peptide2cds,
                                             cds_sequences = cds_sequences )
                overlapping_exons = []

            overlapping_exons.append( e )
            max_to = max( e.mGenomeTo, max_to )

        fragments += buildFragments( overlapping_exons, input, pseudogenes, options,
                                     map_peptide2cds = map_peptide2cds,
                                     cds_sequences = cds_sequences )
        
        map_gene2fragments[gene] = fragments

        # build unaligned sequence
        sequence = options.sequence_separator.join( map( lambda x: x.mSequence, fragments) )

        # sanity check if cds sequence is given: does it correspond to peptide sequence?
        if cds_sequences:
            cds_sequence = options.sequence_separator.join( map( lambda x: x.mCdsSequence, fragments) )
            mapped_sequence = Genomics.translate( cds_sequence, is_seleno = is_seleno )
            if mapped_sequence != sequence:
                if len(mapped_sequence) == len(sequence):
                    ncounts = 0
                    for x in range(len(mapped_sequence)):
                        if mapped_sequence[x].upper() != sequence[x].upper():
                            E.info("gene %s: amino acid position %i has changed from %s to %s" % (gene, x, sequence[x], mapped_sequence[x]))
                            ncounts += 1
                    if ncounts > 0:
                        E.warn( "gene %s: %i amino acid positions have changed" % (gene, ncounts))
                        nwarnings_sequence += 1
                else:
                    E.warn("gene %s: back-translated sequence is different from original" % (gene))
                    E.warn("      Original : %s" % (sequence))
                    E.warn("      Backtrans: %s" % (mapped_sequence))
                    E.warn("assumed to be pseudogenes and ignored." )
                    nwarnings_sequence += 1

                    if options.strict:
                        raise ValueError("gene %s: back-translated sequence has different length from original\n" % (gene))

        unaligned.addSequence( gene, 0, len(sequence), sequence )
        unaligned_cds.addSequence( gene, 0, len(cds_sequence), cds_sequence )
        
    if "unaligned_aa" in options.output:
        writeToFile( unaligned, "unaligned_aa", options, is_aligned = False )

    if "unaligned_na" in options.output:
        writeToFile( unaligned_cds, "unaligned_na", options, is_aligned = False )

    if options.stop_at == "unaligned":
        unaligned.writeToFile( options.stdout, format=options.output_format )
        E.Stop()
        sys.exit(0)

    ##########################################################################
    ##########################################################################
    ##########################################################################
    ## Perform the multiple alignment
    ##########################################################################        
    muscle = WrapperMuscle.Muscle()

    aligned = muscle.Run( unaligned )

    # substitute U for X in selenoproteins
    # note that this is a greedy substitution - genuine stop-codons
    # will be overridden as well.
    if selenoproteins:
        for gene_id, s in aligned.items():
            for pid in map_gene2ids[gene_id]:
                if pid in selenoproteins:
                    s.mString = s.mString.replace( "X", "U" ).replace( "x", "u")
                    break

    if "aligned_aa" in options.output:
        writeToFile( aligned, "aligned_aa", options)

    # perform sanity checks
    # all aligned sequences should have the same length
    # all sequences should be identical to the unaligned sequences
    width = aligned.getWidth()
    for key,val in aligned.items():
        if len(val.mString) != width:
            raise ValueError("incompatible lenghts in %s: %i should be %i." % (key, len(val.mString), width ))
        try:
            s1 = re.sub( "[%s]" % options.gap_chars, "", unaligned[key] )
        except KeyError:
            continue
        s2 = re.sub( "[%s]" % options.gap_chars, "", val.mString )
        if s1.upper() != s2.upper():
            raise ValueError("sequence changed in %s:\nold=%s\nnew=%s" % (key, s1, s2) )

    if options.stop_at == "aligned":
        aligned.writeToFile( options.stdout, format=options.output_format)
        E.Stop()
        sys.exit(0)

    ##########################################################################
    ##########################################################################
    ##########################################################################
    ## output the packed alignment as nucleotides 
    ## The output does not contain any frameshifts that might have been
    ## present in the original sequences
    if "aligned_na" in options.output:

        aligned_cds = Mali.Mali()

        for id in aligned.getIdentifiers():
            entry = aligned.getEntry(id)

            s = []
            c = 0
            gc = aligned.mGapChar
            cds_sequence = unaligned_cds[id]
            for x in range(len(entry.mString)):
                if entry.mString[x] == gc:
                    s.append( gc * 3 )
                else:
                    s.append( cds_sequence[c:c+3] )
                    c += 3

            aligned_cds.addSequence( id, entry.mFrom * 3, entry.mTo * 3, "".join(s) )

        writeToFile( aligned_cds, "aligned_na", options )
        
    ##########################################################################
    ##########################################################################
    ##########################################################################
    ## unpack the protein level alignment
    ## columns stay the same, but alignments with multiple transcripts per gene
    ## are unpacked.
    ##########################################################################        
    unpacked = Mali.Mali()

    for id in aligned.getIdentifiers():
        
        gene = id
        sequence = aligned[id]
        entry = aligned.getEntry(id)

        fragments = map_gene2fragments[gene]
        ## split aligned sequence in its constituent parts
        transitions = []
        c = 0
        for x in fragments:
            c += len(x.mSequence)
            transitions.append(c)

        segments = entry.getSegments( transitions )

        if len(segments) != len(transitions):

            for x in range(min(len(segments), len(transitions))):

                segment = sequence[segments[x][0]:segments[x][1]]

                if options.loglevel >= 8:
                    options.stdlog.write("# %s: transition=%i segment=%i-%i members=%s\n# %s: sequence=%s segment=%s\n" % (id,
                                                                                                                           transitions[x],
                                                                                                                           segments[x][0],
                                                                                                                           segments[x][1],
                                                                                                                           fragments[x].mMembers,
                                                                                                                           id,
                                                                                                                           fragments[x].mSequence,
                                                                                                                           segment ) )



            print "# %s: nfragments=%i" % (id, len(fragments))
            print "# %s: nsegments=%i segments=%s lsequence=%i" % (id, len(segments), str(segments), len(sequence))
            print "# %s: ntransitions=%i transitions=%s lgene=%i" % (id, len(transitions), transitions, len(unaligned[gene]))
            print sequence
            print unaligned[gene]

        assert (len(segments) == len(transitions))

        sequences = {}          
        for member in map_gene2ids[gene]:
            sequences[member] = []

        fragments = map_gene2fragments[id]

        for x in range(len(segments)):

            segment = sequence[segments[x][0]:segments[x][1]]

            if options.loglevel >= 8:
                options.stdlog.write("# %s: transition=%i segment=%i-%i members=%s\n# %s: sequence=%s segment=%s\n" % (id,
                                                                                                                       transitions[x],
                                                                                                                       segments[x][0],
                                                                                                                       segments[x][1],
                                                                                                                       fragments[x].mMembers,
                                                                                                                       id,
                                                                                                                       fragments[x].mSequence,
                                                                                                                       segment ) )


            added = set()

            ## if it was a consensus string, deconvolute with multiple alignment
            if fragments[x].mMali:
                for member in fragments[x].mMembers:
                    s1 = fragments[x].mMali.getSequence( member ).mString
                    s = []
                    i = 0
                    for ch in segment:
                        if ch in options.gap_chars:
                            s.append( ch )
                        else:
                            s.append(s1[i])
                            i += 1
                    sequences[member].append( "".join(s) )
                    added.add( member )
            else:
            ## simply add the segment
                for member in fragments[x].mMembers:
                    sequences[member].append( segment )
                    added.add( member )

            segment = options.gap_char * (len(segment))

            for member in map_gene2ids[gene]:
                if member not in added:
                    sequences[member].append( segment )

        for member in map_gene2ids[gene]:
            s = "".join(sequences[member])

            if len(re.sub("[%s]" % options.gap_chars, "", s )) == 0:
                E.warn("empty sequence for %s." % member )
                nwarnings_empty += 1
                continue    

            unpacked.addSequence( member, 0, -1, s )

    if options.stop_at == "unpacked":
        aligned.writeToFile( options.stdout, format=options.output_format)
        E.Stop()
        sys.exit(0)

    # perform sanity checks
    # all aligned sequences should have the same length
    # all sequences should be identical to the unaligned sequences
    #   except: sometimes at split codons exon boundaries might
    #   be chosen differently. Thus: allow single residue changes
    
    width = unpacked.getWidth()
    nsubstitutions = 0
    
    for key,val in unpacked.items():
        
        if len(val.mString) != width:
            raise ValueError("incompatible mali lengths in %s: %i should be %i." % (key, len(val.mString), width ))
        
        try:
            sold = re.sub( "[%s]" % options.gap_chars, "", input[key] )
        except KeyError:
            continue
        
        snew = re.sub( "[%s]" % options.gap_chars, "", val.mString )

        if len(sold) != len(snew):
            if key in pseudogenes:
                E.warn("pseudogenic sequence %s changed in length" % (key))
                E.warn("%s:\nold=%s\nnew=%s" % (key, sold, snew))
                nwarnings_length += 1
            elif (len(snew) - len(sold) == 1 and sold == snew[:-1]) or \
                 (len(sold) - len(snew) == 1 and snew == sold[:-1]):
                E.warn("sequence %s changed in length - last residue lost" % (key))
                E.warn("%s:\nold=%s\nnew=%s" % (key, sold, snew))
                nwarnings_length += 1
            elif not options.force_map:
                raise ValueError("sequence changed in length %s:\nold=(%i)\n%s\nnew=(%i)\n%s" % (key, len(sold), sold, len(snew), snew))

        else:
            ## for those with equal length: fix single residue substitutions
            nchanged = 0

            if sold != snew:
                ## fix single character changes
                c = 0
                s = list(val.mString)
                for x in range(len(s)):
                    if s[x] in options.gap_chars:
                        continue
                    if s[x].upper() != sold[c].upper():
                        if options.loglevel >= 8:
                            options.stdlog.write("# %s: fixing split codon: %i->%i %s:%s:%s -> %s:%s:%s\n" %\
                                                     (key, x, c, 
                                                      sold[c-5:c], sold[c], sold[c+1:c+6],
                                                      "".join(s[x-5:x]), s[x], "".join(s[x+1:x+6:])) )
                        nchanged += 1
                        s[x] = sold[c]

                    c += 1

                val.mString = "".join(s)

            if nchanged > 0:
                E.info( "%s: fixed %i split codons." % (key, nchanged) )

            # if more substitutions made than exon boundaries
            # it is likely a programming error
            # should be >=, but was too strict for one fly sequence, so I relaxed it.
            if nchanged > len( exons[key] ):
                raise ValueError("more codons fixed than intron-exon boundaries in sequence %s: %i>=%i" % (key, nchanged, len(exons[key])))
            
        snew = re.sub( "[%s]" % options.gap_chars, "", val.mString )
        
        if sold.upper() != snew.upper():
            if key in pseudogenes:
                E.warn("pseudogenic sequence %s changed." % (key))
                E.warn("%s:\n# old=%s\n# new=%s" % (key, sold, snew))
                nwarnings_sequence += 1
                
            if options.force_map:
                map_old2new = alignlib.makeAlignmentVector()
                alignator = alignlib.makeAlignatorDPFull( alignlib.ALIGNMENT_GLOBAL, -10.0, -1.0 )
                s1 = alignlib.makeSequence( sold )
                s2 = alignlib.makeSequence( snew )
                alignator.align( map_old2new, s1, s2 )

                val.threadSequence( snew, map_old2new )
                E.warn("sequence threaded onto original sequence: %s" % (key))
                nsubstitutions += 1

            else:    
                raise ValueError("# sequence changed in %s:\n# old=%s\n# new=%s" % (key,
                                                                                    sold,
                                                                                    snew))

    if not options.output: 
        unpacked.writeToFile( options.stdout, format=options.output_format )
    else:
        if "final_aa" in options.output:
            writeToFile( unpacked, "final_aa", options )

        if "final_na" in options.output:

            ## map alignment to original cds
            ## frameshifts in the original cdna are ignored
            ## output table with coordinates
            unpacked_cds = Mali.Mali()
            
            info, all_coords = [], []
            
            for key in unpacked.getIdentifiers():

                entry = unpacked.getEntry(key)
                
                p = unpacked[key]
                c = cds_sequences[key]
                
                ###############################
                # build map of cds to alignment position (map_cds2pos)
                try:
                    E.debug("building map for %s" % key)
                    map_pos2cds = getMapPeptide2Cds( p, c,
                                                     options )

                except ValueError, msg:
                    E.warn("final alignment not mappable for %s -skipped" % key )
                    nskipped += 1
                    continue


                #################################
                # map alignment pos to genome
                map_c2g = map_cds2genome[key]
                map_pos2genome = alignlib.makeAlignmentBlocks()
                alignlib.combineAlignment( map_pos2genome, map_pos2cds, map_c2g, alignlib.CR )

                coords = []            
                is_negative_strand, start = genome_starts[key]
                for x in range( len(p) * 3 ):
                    y = map_pos2genome.mapRowToCol(x)
                    if y < 0:
                        coords.append( "" )
                    else:
                        if is_negative_strand:
                            coords.append( 0 - (start + y) )
                        else:
                            coords.append( start + y )
                all_coords.append( coords )
                info.append( ( key, exons[key][0].mSbjctToken, exons[key][0].mSbjctStrand ) )

                ################################
                # create threaded alignment string
                map_pos2cds.switchRowCol()
                map_cds2pos = map_pos2cds

                alignatum = alignlib.makeAlignatum( cds_sequences[key] )

                alignatum.mapOnAlignment( map_cds2pos, len(p) * 3 )
                s = alignatum.getString()
                if len(s) != len(p) * 3:
                    raise ValueError("incomplete aligned string for %s: %s, cds=%s" % (key, s, c ))
                
                unpacked_cds.addSequence( key, 0, len(c), s )

            writeToFile( unpacked_cds, "final_na", options )

            if "coords" in options.output:
                E.info("output genomic coordinates to %s" % options.output_filename_coords)
                with IOTools.openFile( options.output_filename_coords, "w") as outfile:
                    outfile.write( "position\t%s\n" % "\t".join( [ "|".join( map(str,i) ) for i in info ] ) )
                    all_coords = zip( *all_coords )
                    for x, c in enumerate( all_coords ):
                        outfile.write("%i\t%s\n" % (x, "\t".join(map(str, c) ) ) )
                    
    E.info( "ninput=%i, noutput=%i, nskipped=%i, nwarnings_length=%i, nwarnings_sequences=%i, npseudogenes=%i, nsubstitutions=%i" % \
                (input.getLength(),
                 unpacked.getLength(),
                 nskipped,
                 nwarnings_length,
                 nwarnings_sequence,
                 len(pseudogenes),
                 nsubstitutions))
    
    E.Stop()

    

