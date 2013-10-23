################################################################################
#   Gene prediction pipeline 
#
#   $Id: Exons.py 2881 2010-04-07 08:45:38Z andreas $
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
"""
Exons.py - A library to read/write/manage exons.
=====================================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

"""

import re, sys, string
try: import alignlib_lite
except ImportError: pass

class Exon:
    """class for exons.

    contains info about the genomic location of an exon
    and its location within a peptide sequence.

    The field mAlignment is set optionally.
    """
    def __init__(self ):
        self.mQueryToken = ""
        self.mSbjctToken = ""
        self.mSbjctStrand = ""
        self.frame = 0
        self.mRank = 0
        self.mPeptideFrom = 0
        self.mPeptideTo = 0
        self.mGenomeFrom = 0
        self.mGenomeTo = 0
        self.mAlignment = "U"
        
    def Read( self, line, contig_sizes = {}, format = "exons", extract_id = None, converter = None ):
        """read exon from tab-separated line.

        extract_id is a regular expression object to extract the identifier from
        the identifier column.

        if converter is given, it is used to convert to
        zero-based open-closed both strand coordinates.
        """

        if format == "exons":
            (self.mQueryToken,
             self.mSbjctToken,
             self.mSbjctStrand,
             self.frame,
             self.mRank,
             self.mPeptideFrom,
             self.mPeptideTo,
             self.mGenomeFrom,
             self.mGenomeTo) = line[:-1].split("\t")
        elif format == "gtf":
            (self.mSbjctToken,
             a, b,
             self.mGenomeFrom,
             self.mGenomeTo,
             c, 
             self.mSbjctStrand,
             self.frame,
             self.mQueryToken ) = line[:-1].split("\t")
            
            self.mRank = 0
            self.mPeptideFrom, self.mPeptideTo = 0, 0
            self.mGenomeFrom = int(self.mGenomeFrom) - 1

        if extract_id:
            self.mQueryToken = extract_id.search( self.mQueryToken).groups()[0]
        
        if self.frame == ".": self.frame = 0
            
        (self.mRank, self.frame,
         self.mPeptideFrom, self.mPeptideTo,
         self.mGenomeFrom, self.mGenomeTo ) = map( int,\
                                                   (self.mRank, self.frame,
                                                    self.mPeptideFrom, self.mPeptideTo,
                                                    self.mGenomeFrom, self.mGenomeTo ))

        if self.mSbjctStrand == "1": self.mSbjctStrand = "+"
        elif self.mSbjctStrand in ("0", "-1"): self.mSbjctStrand = "-"

        if converter and self.mSbjctToken in contig_sizes:
            self.mGenomeFrom, self.mGenomeTo = converter( self.mGenomeFrom,
                                                          self.mGenomeTo,
                                                          self.mSbjctStrand == "+",
                                                          contig_sizes[self.mSbjctToken] )
        
        if self.mSbjctStrand == "-" and self.mSbjctToken in contig_sizes:
            sbjct_length = contig_sizes[self.mSbjctToken]                    
            self.mGenomeFrom, self.mGenomeTo = sbjct_length - self.mGenomeTo, sbjct_length - self.mGenomeFrom

        if self.mGenomeFrom > self.mGenomeTo:
            self.mGenomeFrom, self.mGenomeTo = self.mGenomeTo, self.mGenomeFrom

        return 1

    def __str__( self ):
        return string.join( map(str, (
            self.mQueryToken,
            self.mSbjctToken,
            self.mSbjctStrand,
            self.frame,
            self.mRank,
            self.mPeptideFrom,
            self.mPeptideTo,
            self.mGenomeFrom,
            self.mGenomeTo)), "\t" )

    def Merge( self, other ):
        """Merge this exon with another (adjacent and preceeding) exon.

        Do not merge if the distance between exons is not divisible by 3.
        Merging of two exons invalidated peptide coordinates for all following
        exons. These need to be updated.
        """
        
        if other.mSbjctStrand != self.mSbjctStrand or \
           other.mSbjctToken != self.mSbjctToken:
            raise ValueError, "exons not on the same contig and strand."

        if other.mGenomeTo >= self.mGenomeFrom:
            raise ValueError, "other exon not preceding this exon"

        difference_sbjct = self.mGenomeFrom - other.mGenomeTo
        difference_query = 0
        
        if difference_sbjct % 3 != 0:
            raise ValueError, "can not merge exons with incompatible phases."
        
        self.mGenomeFrom = other.mGenomeFrom
        self.mPeptideFrom = other.mPeptideFrom
        self.mPeptideTo += difference_sbjct / 3
        self.frame = other.frame
        
        ## remove split codons and add to difference
        ## introduce a gap in the peptide sequence as well
        if other.mAlignment[-1][0] == "S":
            del other.mAlignment[-1]
            del self.mAlignment[0]
            difference_sbjct += 3
            difference_query += 3

        extra = [("G", 0, difference_sbjct)]
        if difference_query:
            extra.append( ("G", difference_query, 0) )

        self.mAlignment = other.mAlignment + extra + self.mAlignment

    def InvertGenomicCoordinates( self, lgenome):
        """invert genomic alignment on sequence.

        Negative strand is calculated from the other end.
        """
        if self.mSbjctStrand == "-":
            x = min(self.mGenomeFrom, self.mGenomeTo)
            y = max(self.mGenomeFrom, self.mGenomeTo)
                
            self.mGenomeFrom = lgenome - y
            self.mGenomeTo   = lgenome - x
            
    def GetCopy( self ):
        
        e = Exon()

        (e.mQueryToken,
         e.mSbjctToken,
         e.mSbjctStrand,
         e.frame,
         e.mRank,
         e.mPeptideFrom,
         e.mPeptideTo,
         e.mGenomeFrom,
         e.mGenomeTo) = \
         (self.mQueryToken,
          self.mSbjctToken,
          self.mSbjctStrand,
          self.frame,
          self.mRank,
          self.mPeptideFrom,
          self.mPeptideTo,
          self.mGenomeFrom,
          self.mGenomeTo)
        return e

##----------------------------------------------------------
def UpdatePeptideCoordinates( exons ):
    """updates peptides coordinates for a list of exons.

    Exons have to  be sorted.
    """
    last_e = exons[0]
    
    for e in exons[1:]:
        if e.mPeptideFrom != last_e.mPeptideTo:
            offset = last_e.mPeptideTo - e.mPeptideFrom
            e.mPeptideFrom += offset
            e.mPeptideTo += offset
        last_e = e
        
##----------------------------------------------------------
def PostProcessExons( all_exons,
                      do_invert = None,
                      remove_utr = None,
                      filter = None,
                      reset = False,
                      require_increase = False,
                      no_invert = False,
                      contig_sizes = {},
                      from_zero = False,
                      delete_missing = False,
                      set_peptide_coordinates = False,
                      set_rank = False):
    """do post-processing of exons

    exons is a dictionary of lists of exons.

    Exons are sorted by mPeptideFrom.
    
    Operations include:
    
    -invert: sort out forward/reverse strand coordinates

    -set_peptide_coordinates: sets the peptide coordinates of
        exons.

    -set-rank: set rank of exons
    
    -remove_utr: remove any utr (needs peptide coordinates)

    -delete_missing: if set set true, exons on contigs not in contig_sizes
        will be deleted.
        
    -from_zero: exon genomic coordinates start at 0

    -reset: exon genomic coordinates start 0
    """
    
    for k in all_exons.keys():

        exons = all_exons[k]
        if delete_missing and contig_sizes and exons[0].mSbjctToken not in contig_sizes:
            del all_exons[k]
            continue
        
        invert = False

        ## ENSEMBL PATCH: exons on different strands, for example CG32491!
        # if exons[0].mGenomeFrom > exons[-1].mGenomeFrom and \
        #       exons[0].mSbjctStrand == exons[-1].mSbjctStrand:
        #    invert = True

        ## convert to forward/reverse strand coordinates if so desired
        if (invert or do_invert or from_zero) and not no_invert:

            if contig_sizes:
                if exons[0].mSbjctToken in contig_sizes:
                    l = contig_sizes[exons[0].mSbjctToken]
                elif "dummy" in contig_sizes:
                    l = contig_sizes["dummy"]
                else:
                    continue
                
            elif from_zero:
                l = 0
            else:
                max_from = max( map(lambda x: x.mGenomeFrom, exons))
                max_to = max( map(lambda x: x.mGenomeTo, exons))
                l = max(max_from, max_to)

            for exon in exons:
                exon.InvertGenomicCoordinates( l )

        ## set peptide coordinates. Exon genomic coordinates are now forward/reverse strand
        ## coordinates so sorting is straight-forward
        if set_peptide_coordinates:
            exons.sort( lambda x,y: cmp( x.mGenomeFrom, y.mGenomeFrom) )
            start = 0
            for exon in exons:
                exon.mPeptideFrom = start
                start += exon.mGenomeTo - exon.mGenomeFrom
                exon.mPeptideTo = start
        else:
            exons.sort( lambda x,y: cmp( x.mPeptideFrom, y.mPeptideFrom) )

        ## set rank
        if set_rank:
            rank = 1
            for exon in exons:
                exon.mRank = rank
                rank += 1

        if remove_utr:
            # remove 5' UTR
            exon = exons[0]
            exon.mGenomeFrom = exon.mGenomeTo - (exon.mPeptideTo - exon.mPeptideFrom)
            # remove 3' UTR
            exon = exons[-1]
            exon.mGenomeTo = exon.mGenomeFrom + (exon.mPeptideTo - exon.mPeptideFrom)

        ## patch for last exon
        ## Occasionally I found an incomplete codon
        ## at the end of a gene. The peptide sequence
        ## was ok (apart from including the stop codon?)
        ## Example: ENSGALP00000006767
        ## NB: happened only in chicken, why?
        ## Check if last codon is complete, otherwise
        ## add appropriate residues to genome_to
        ## dangerous: might be due to pseudogene!
##         e = exons[-1]
##         if e.mPeptideTo % 3 != 0:
##             d = 3 - e.mPeptideTo % 3
##             e.mPeptideTo += d
##             e.mGenomeTo += d

        if reset:
            offset = exons[0].mGenomeFrom
            for e in exons:
                e.mGenomeFrom -= offset
                e.mGenomeTo -= offset
                
    return all_exons

##----------------------------------------------------------
def GetExonBoundariesFromTable( dbhandle,
                                table_name_predictions = "predictions",
                                table_name_exons = "exons",
                                only_good = False,
                                do_invert = None,
                                remove_utr = None,
                                filter = None,
                                reset = False,
                                require_increase = False,
                                contig_sizes = {},
                                prediction_ids = None,
                                table_name_quality = "quality",
                                table_name_redundant = "redundant",
                                non_redundant_filter = False,
                                schema = None,
                                quality_filter = None,
                                from_zero = False,
                                delete_missing = False ):
    """get exon boundaries from table."""

    extra = ""
    extra_tables = ""
    
    if only_good:
        extra += " AND e.is_ok = TRUE"

    if schema:
        table_name_predictions = "%s.%s" % (schema, table_name_predictions)
        table_name_exons = "%s.%s" % (schema, table_name_exons)                
        table_name_redundant = "%s.%s" % (schema, table_name_redundant)
        table_name_quality = "%s.%s" % (schema,table_name_quality)
        
    if prediction_ids:
        extra += " AND p.prediction_id IN ('%s')" % ("','".join(map(str,prediction_ids)))

    if quality_filter:
        extra_tables += ", %s AS q" % table_name_quality
        extra += " AND q.prediction_id = p.prediction_id AND q.class IN ('%s')" % ("','".join(quality_filter))

    if non_redundant_filter:
        extra_tables += ", %s AS r" % table_name_redundant
        extra += " AND r.rep_prediction_id = r.mem_prediction_id AND r.rep_prediction_id = p.prediction_id"
                                                  
    statement = """
    SELECT DISTINCT p.prediction_id, p.sbjct_token, p.sbjct_strand,
    e.exon_frame, 0, e.exon_from, e.exon_to, e.genome_exon_from, e.genome_exon_to
    FROM %s AS p, %s AS e %s
    WHERE p.prediction_id = e.prediction_id AND e.exon_to > 0
    %s
    ORDER BY p.prediction_id, e.genome_exon_from
    """ % (table_name_predictions, table_name_exons, extra_tables, extra )
    cc = dbhandle.cursor()
    cc.execute(statement)
    result = cc.fetchall()
    cc.close()
    
    all_exons = {}
    
    for r in result:
        e = Exon()
        e.Read( "\t".join( map(str, r)) + "\n", contig_sizes)

        if filter and not e.mQueryToken in filter: continue

        if not all_exons.has_key( e.mQueryToken ): all_exons[e.mQueryToken] = []
        all_exons[e.mQueryToken].append( e )
        
    return PostProcessExons( all_exons,
                             do_invert,
                             remove_utr,
                             filter,
                             reset,
                             require_increase,
                             contig_sizes,
                             from_zero = from_zero,
                             delete_missing = delete_missing )


##----------------------------------------------------------
def CountNumExons( exons ):
    """return hash with number of exons per entry."""

    nexonst = {}
    for k, ee in exons.items():
        nexons[k] = len(ee)
    return nexons
        
##----------------------------------------------------------
def SetRankToPositionFlag( exons ):
    """set rank for all exons.

    Set rank to
    1 : if it is first exon,
    -1: if it is the last exon (single exon genes are -1)
    0 : if it is an internal exon.
    """

    for k, ee in exons.items():
        if len(ee) == 1:
            ee[0].mRank = -1
            continue
        ee.sort( lambda x,y: cmp( x.mPeptideFrom, y.mPeptideFrom ) )
        ee[0].mRank = 1
        ee[-1].mRank = -1
        for x in range(1, len(ee) - 1):
            ee[x].mRank = 0

##----------------------------------------------------------
def ReadExonBoundaries( file,
                        do_invert = None,
                        remove_utr = None,
                        filter = None,
                        reset = False,
                        require_increase = False,
                        no_invert = False,
                        contig_sizes = {},
                        converter = None,
                        from_zero = False,
                        delete_missing = False,
                        format = "exons",
                        gtf_extract_id = None ):
    """read exons boundaries from tab separated file.

    if remove_utr is set, the UTR of the first/last exon is removed.

    if reset is set, then the genomic part is moved so that it starts at 1.
    if require_increase is set, then exons are sorted in increasing order.

    If do_invert is set: negative strand coordinates are converted to positive
    strand coordinates

    if no_invert is set: coordinates are kept as they are.

    if from_zero is set: coordinates are mapped from 0. Thus reverse
    strand coordinates will be negative.

    if delete_missing is True and sbjct-token is not in contig_sizes
    but the exon needs to be inverted: delete transcript.

    The exon file format is tab-separated and can be of the two formats:

    format="exons":
    id, contig, strand, frame, rank, peptide_from, peptide_to, genome_from, genome_to

    format = "gtf":
    contig, ignored, ignored, genome_from, genome_to, ignored, strand, frame, id

    if converter is given, use it to convert to forward/reverse strand coordinates.

    gtg_extract_id: regular expression object to extract id from id column.
    """

    all_exons = {}

    l = None
    for line in file:
        if line[0] == "#": continue
        if line[0] == "": continue
        e = Exon()

        e.Read(line,
               contig_sizes = contig_sizes,
               format = format,
               converter = converter,
               extract_id = gtf_extract_id )
        
        l = e
        if filter and not e.mQueryToken in filter: continue
        if not all_exons.has_key( e.mQueryToken ): all_exons[e.mQueryToken] = []
        all_exons[e.mQueryToken].append( e )

    set_rank = False
    set_peptide_coordinates = False
    if format == "gtf":
        set_rank = True
        set_peptide_coordinates = True

    return PostProcessExons( all_exons,
                             do_invert = do_invert,
                             remove_utr = remove_utr,
                             filter = filter,
                             reset = reset,
                             require_increase = require_increase,
                             no_invert = no_invert,
                             contig_sizes = contig_sizes,
                             from_zero = from_zero,
                             delete_missing = delete_missing,
                             set_peptide_coordinates = set_peptide_coordinates,
                             set_rank = set_rank )

    
##------------------------------------------------------------
def Alignment2Exons( alignment, query_from = 0, sbjct_from = 0, add_stop_codon = 1):
    """convert a Peptide2DNA alignment to exon boundaries.
    """

    exons = []

    ## count in nucleotides for query
    query_from *= 3
    query_pos = query_from
    sbjct_pos = sbjct_from
    frame = 0
    ali = []
    
    for state, l_query, l_sbjct in alignment:

        if state in ("S", "M", "G", "F"):
            ali.append( (state, l_query, l_sbjct))

        ## count as nucleotides
        l_query *= 3
                
        if state == "M":
            query_pos += l_query

        elif state == "G":
            query_pos += l_query
            
        elif state == "P":
            ## state P is counted as an intron
            frame = query_from % 3
            if frame != 0: frame = 3 - frame
            exon = Exon()
            exon.mPeptideFrom = query_from
            exon.mPeptideTo = query_pos
            exon.mSbjctFrame = frame
            exon.mGenomeFrom = sbjct_from
            exon.mGenomeTo = sbjct_pos
            exon.mAlignment = ali
            exons.append( exon )
            ali = []
            query_pos += l_query
            query_from = query_pos
            sbjct_from = sbjct_pos + l_sbjct

        elif state == "S":
            query_pos += l_sbjct

        elif state == "I":
            pass
        
        elif state == "5":
            frame = query_from % 3
            if frame != 0: frame = 3 - frame
            exon = Exon()
            exon.mPeptideFrom = query_from
            exon.mPeptideTo = query_pos
            exon.frame = frame
            exon.mGenomeFrom = sbjct_from
            exon.mGenomeTo = sbjct_pos
            exon.mAlignment = ali
            exons.append( exon )
            ali = []
            query_pos += l_query
            
        elif state == "3":
            query_pos += l_query
            query_from = query_pos
            sbjct_from = sbjct_pos + l_sbjct
            
        sbjct_pos += l_sbjct
                
    ## add three for the stop codon:
    if add_stop_codon:
        query_pos += 3
        sbjct_pos += 3
        
    frame = query_from % 3
    if frame != 0: frame = 3 - frame

    exon = Exon()
    exon.mPeptideFrom = query_from
    exon.mPeptideTo = query_pos
    exon.frame = frame
    exon.mGenomeFrom = sbjct_from
    exon.mGenomeTo = sbjct_pos
    exon.mAlignment = ali
    
    exons.append( exon )

    return exons

##-------------------------------------------------------------------
def Exons2Alignment( exons ):
    """build alignment string from a (sorted) list of exons.
    """
    alignment = []
    last_e = exons[0]

    for e in exons[1:]:

        alignment += last_e.mAlignment

        difference = e.mGenomeFrom - last_e.mGenomeTo
        
        alignment.append( ("5", 0, 2) )
        alignment.append( ("I", 0, difference - 4 ) )
        alignment.append( ("3", 0, 2) )
        last_e = e

    alignment += last_e.mAlignment 
    return alignment

##------------------------------------------------------------
class ComparisonResult:

    class Link:
        def __init__(self):
            self.mId1 = 0
            self.mId2 = 0
            self.mPercentIdentity = 0
            self.mPercentSimilarity = 0
            self.mCoverage = 0
            self.mOverlap = 0
            self.mIsSameFrame = True
            self.mIsGoodExon = True

        def __str__( self ):

            return ";".join( map(str, (self.mId1, self.mId2, self.mPercentIdentity, self.mPercentSimilarity,
                                       self.mCoverage, self.mOverlap,
                                       self.mIsSameFrame,
                                       self.mIsGoodExon )))
                                 
    def __init__(self):

        self.mNumDifferenceExons = 0
        self.mSumBoundaryDifferences = 0
        self.mMaxBoundaryDifferences = 0
        self.mNumIdenticalExons  = 0
        self.mNumMissedCmpBoundaries = 0
        self.mNumMissedRefBoundaries = 0
        self.mNumCmpBoundaries = 0
        self.mNumRefBoundaries = 0         
        self.mNumDubiousExons    = 0
        self.mNumDeletedExons    = 0
        self.mNumInsertedExons   = 0
        self.mNumDeletedIntrons  = 0
        self.mNumInsertedIntrons = 0
        self.mNumTruncatedNExons = 0
        self.mNumTruncatedCExons = 0
        self.mNumExtendedNExons = 0
        self.mNumExtendedCExons = 0
        self.mNumDeletedNExons   = 0
        self.mNumDeletedCExons   = 0
        self.mNumInsertedNExons  = 0
        self.mNumInsertedCExons  = 0
        self.mNumSkippedExons    = 0
        self.mEquivalences = []

    def __str__(self):
        a = string.join( map(str, ( \
            self.mNumDifferenceExons,
            self.mNumIdenticalExons,
            self.mNumSkippedExons,
            self.mNumMissedRefBoundaries,                        
            self.mNumMissedCmpBoundaries,            
            self.mNumCmpBoundaries,
            self.mNumRefBoundaries,            
            self.mSumBoundaryDifferences,
            self.mMaxBoundaryDifferences,
            self.mNumDubiousExons    ,
            self.mNumDeletedExons    ,
            self.mNumInsertedExons   ,
            self.mNumDeletedIntrons  ,
            self.mNumInsertedIntrons ,
            self.mNumTruncatedNExons ,
            self.mNumTruncatedCExons ,
            self.mNumExtendedNExons ,
            self.mNumExtendedCExons ,
            self.mNumDeletedNExons   ,
            self.mNumDeletedCExons   ,
            self.mNumInsertedNExons  ,
            self.mNumInsertedCExons )), "\t")
        
        b = " ".join(map(str, self.mEquivalences))
        
        return a + "\t" + str(b)

    def Pretty( self, prefix = "# " ):
        lines = []
        for k,v in self.__dict__.items():
            if k[0] == 'm':
                if k == "mEquivalences":
                    for vv in v:
                        lines.append( "%s%-44s: %s" % (prefix, k, str(vv)) )
                else:
                    lines.append( "%s%-40s: %s" % (prefix, k, str(v)) )
        return string.join(lines, "\n")

    def GetHeader( self ):
        """return header line for tab-separated column format."""
        return string.join( (
            "NumDifferenceExons",
            "NumIdenticalExons",
            "NumSkippedExons",
            "NumMissedRefBoundaries",                        
            "NumMissedCmpBoundaries",            
            "NumCmpBoundaries",
            "NumRefBoundaries",            
            "SumBoundaryDifferences",
            "MaxBoundaryDifferences",
            "NumDubiousExons",
            "NumDeletedExons",
            "NumInsertedExons",
            "NumDeletedIntrons",
            "NumInsertedIntrons",
            "NumTruncatedNExons",
            "NumTruncatedCExons",
            "NumExtendedNExons",
            "NumExtendedCExons",
            "NumDeletedNExons",
            "NumDeletedCExons",
            "NumInsertedNExons",
            "NumInsertedCExons"), "\t")
        
def RemoveRedundantEntries( l ):
    """remove redundant entries (and 0s) from list.

    One liner?
    """

    if len(l) == 0:
        return l
    
    l.sort()
    last = l[0]
    n = [last]
    for x in l[1:]:
        if x != last and x > 0:
            n.append( x )
        last = x
    return n

def CompareGeneStructures( xcmp_exons, ref_exons, 
                           map_ref2cmp = None,
                           cmp_sequence = None,
                           ref_sequence = None,
                           threshold_min_pide = 0,
                           threshold_slipping_exon_boundary = 9,
                           map_cmp2ref = None,
                           threshold_terminal_exon = 15):
    
    """Compare two gene structures.

    This function is useful for comparing the exon boundaries of
    a predicted peptide with the exon boundaries of the query peptide.

    cmp_exons are exons for the gene to test.
    ref_exons are exons from the reference.
    
    Exon boundaries are already mapped to the peptide for the
    reference.
    
    map_ref2cmp: Alignment of protein sequences for cmp and ref.
    map_cmp2ref: Alignment of cmp to ref. If given, mapping is done from cmp to ref.
    Invalid exon boundaries can be set to -1.

    threshold_terminal_exon:
        Disregard terminal exons for counting missed boundaries, if they are
        maximum x nucleotides long.
    
    
    """

    if ref_sequence and cmp_sequence:
        cmp_seq = alignlib_lite.py_makeSequence( cmp_sequence)
        ref_seq = alignlib_lite.py_makeSequence( ref_sequence)
    else:
        ref_seq, cmp_seq = None, None
        
    result = ComparisonResult()
    # unmappable exons have to be kept for counting purpses.
    if map_cmp2ref:
        cmp_exons = MapExons( xcmp_exons, map_cmp2ref )
    else:
        cmp_exons = xcmp_exons

    cmp_len = cmp_exons[-1].mPeptideTo
    ref_len = ref_exons[-1].mPeptideTo    
    
    # default values
    result.mNumDifferenceExons = len(ref_exons) - len(cmp_exons)

    # number of matching boundaries (do not use first and last)
    boundaries_cmp = map( lambda x: x.mPeptideFrom, cmp_exons) + map( lambda x: x.mPeptideTo, cmp_exons)
    boundaries_ref = map( lambda x: x.mPeptideFrom, ref_exons) + map( lambda x: x.mPeptideTo, ref_exons)
    
    boundaries_cmp = RemoveRedundantEntries( boundaries_cmp[1:-1] )
    boundaries_ref = RemoveRedundantEntries( boundaries_ref[1:-1] )    
    
    result.mNumMissedCmpBoundaries = CountMissedBoundaries( boundaries_cmp, boundaries_ref,
                                                            threshold_slipping_exon_boundary,
                                                            threshold_terminal_exon, cmp_len - threshold_terminal_exon + 1 )
    result.mNumMissedRefBoundaries = CountMissedBoundaries( boundaries_ref, boundaries_cmp,
                                                            threshold_slipping_exon_boundary,
                                                            threshold_terminal_exon, ref_len - threshold_terminal_exon + 1 )                                                            
    result.mNumCmpBoundaries = len(boundaries_cmp)    
    result.mNumRefBoundaries = len(boundaries_ref)

    in_sync = False

    # assume that peptides start at 1.
    cmp_from = 1
    ref_from = 1
    cmp_to = 0
    for e in cmp_exons:
        cmp_to = max( cmp_to, e.mPeptideTo)
    cmp_to = cmp_to / 3 + 1
    ref_to = 0
    for e in ref_exons:
        ref_to = max( ref_to, e.mPeptideTo)
    ref_to = ref_to / 3 + 1

    e,r = 0,0

    while e < len(cmp_exons) and r < len(ref_exons):
        
        percent_identity = 0
        percent_similarity = 0
        is_good_exon = False

        exon1 = cmp_exons[e]
        exon2 = ref_exons[r]
        
        if e < len(cmp_exons)-1 :
            next_cmp_exon = cmp_exons[e+1]
        else:
            next_cmp_exon = None

        if r < len(ref_exons) - 1:
            next_ref_exon = ref_exons[r+1]
        else:
            next_ref_exon = None

        ## cmp can have 0,0 entries, these are exons that were impossible to align
        if exon1.mPeptideFrom == 0 and exon1.mPeptideTo == 0:
            if r == 0:
                result.mNumInsertedNExons += 1
            else:
                result.mNumInsertedExons += 1                
            e += 1
        elif exon1.mPeptideFrom < 0 or exon1.mPeptideTo < 0:
            result.mNumSkippedExons += 1
            e += 1
        elif exon2.mPeptideFrom < 0 or exon2.mPeptideTo < 0:
            result.mNumSkippedExons += 1
            r += 1
        elif exon2.mPeptideTo <= exon1.mPeptideFrom + threshold_slipping_exon_boundary:
            ## no overlap
            if e == 0:
                result.mNumDeletedNExons += 1
            else:
                result.mNumDeletedExons += 1
            r += 1
        elif exon1.mPeptideTo <= exon2.mPeptideFrom + threshold_slipping_exon_boundary:
            ## no overlap
            if r == 0:
                result.mNumInsertedNExons += 1
            else:
                result.mNumInsertedExons += 1
            e += 1
        else:
            ## overlap
            dfrom = int(abs(exon1.mPeptideFrom - exon2.mPeptideFrom))
            dto = int(abs(exon1.mPeptideTo - exon2.mPeptideTo))

            overlap  = min(exon1.mPeptideTo,  exon2.mPeptideTo) - max(exon1.mPeptideFrom, exon2.mPeptideFrom)
            coverage = 100 * overlap / (max(exon1.mPeptideTo,  exon2.mPeptideTo) - min(exon1.mPeptideFrom, exon2.mPeptideFrom))
            
            is_good_exon = True

            link = ComparisonResult.Link()
            link.mId1 = e
            link.mId2 = r
            link.mIsGoodExon = is_good_exon
            link.mIsSameFrame = exon1.frame = exon2.frame
            link.mCoverage = coverage
            link.mOverlap = overlap
            
            ## get percent identity
            if cmp_seq and ref_seq and map_ref2cmp:
                
                tmp_ali = alignlib_lite.py_makeAlignmentVector()

                xquery_from = max( exon2.mPeptideFrom / 3, exon1.mPeptideFrom / 3) 
                xquery_to = min(1 + exon2.mPeptideTo / 3, 1 + exon1.mPeptideTo / 3)

                alignlib_lite.py_copyAlignment( tmp_ali, map_ref2cmp, xquery_from, xquery_to )

                percent_identity = alignlib_lite.py_calculatePercentIdentity( tmp_ali,
                                                                      ref_seq,
                                                                      cmp_seq ) * 100
                percent_similarity = alignlib_lite.py_calculatePercentSimilarity( tmp_ali ) * 100

                if percent_identity <= threshold_min_pide:
                    is_good_exon = False

                link.mPercentIdentity = percent_identity
                link.mPercentSimilarity = percent_similarity

            result.mEquivalences.append( link )
            
            ## adjust regions for terminal exons
            if e == 0 and r == 0:
                if dfrom <= (cmp_from - 1) * 3 and dfrom > 0:
                    if is_good_exon:                        
                        result.mNumTruncatedNExons = dfrom
                    dfrom = 0
                elif dfrom <= (ref_from - 1) * 3 and dfrom > 0:                
                    if is_good_exon:                        
                        result.mNumExtendedNExons = dfrom
                    dfrom = 0
                
            if e == len(cmp_exons)-1 and r == len(ref_exons)-1:
                if dto <= (ref_to - cmp_to) * 3 and dto > 0:
                    if is_good_exon:
                        result.mNumTruncatedCExons = dto
                    dto = 0
                elif dto <= (cmp_to - ref_to) * 3 and dto > 0:
                    if is_good_exon:
                        result.mNumExtendedCExons = dto
                    dto = 0

            ## do not count deviations for terminal query exons
            if e == 0 and dfrom <= (cmp_from - 1) * 3 and dfrom > 0:
                dfrom = 0

            if e == len(cmp_exons)-1 and dto <= (ref_to - cmp_to) * 3 and dto > 0:
                dto = 0
                
            ## deal with different boundary conditions:
            if dfrom <= threshold_slipping_exon_boundary and dto <= threshold_slipping_exon_boundary:
                if is_good_exon: result.mNumIdenticalExons += 1
                e += 1
                r += 1
            ## next exon within this exon2_exon
            elif exon1.mPeptideTo < exon2.mPeptideTo and \
                     next_cmp_exon and \
                     next_cmp_exon.mPeptideTo <= exon2.mPeptideTo + threshold_slipping_exon_boundary:
                if is_good_exon: result.mNumInsertedIntrons += 1
                e += 1
                in_sync = True
                dto = 0
            ## next exon2_exon within this exon
            elif exon2.mPeptideTo < exon1.mPeptideTo and \
                     next_ref_exon and \
                     next_ref_exon.mPeptideTo <= exon1.mPeptideTo + threshold_slipping_exon_boundary:
                if is_good_exon: result.mNumDeletedIntrons += 1
                r += 1
                in_sync = True
                dto = 0
            else:
                e += 1
                r += 1
                if in_sync:
                    dfrom = 0

            if is_good_exon:
                result.mSumBoundaryDifferences += dfrom + dto
                result.mMaxBoundaryDifferences = max( dfrom, result.mMaxBoundaryDifferences )
                result.mMaxBoundaryDifferences = max( dto, result.mMaxBoundaryDifferences )                
            else:
                result.mNumDubiousExons += 1
                
    while e < len(cmp_exons):
        e += 1
        result.mNumInsertedCExons += 1

    while r < len(ref_exons):
        r += 1        
        result.mNumDeletedCExons += 1

    return result

###################################################################################
def MapExons( exons, map_a2b ):
    """map peptide coordinates of exons with map.

    returns a list of mapped exons.
    """

    long_map = alignlib_lite.py_makeAlignmentVector()

    for x in range(map_a2b.getRowFrom(), map_a2b.getRowTo() + 1):
        y = map_a2b.mapRowToCol(x)
        if y:
            alignlib_lite.py_addDiagonal2Alignment( long_map, 3 * (x - 1) + 1, 3 * x, 3 * (y - x) )

    def MyMapLeft( a, x):
        while x >= a.getRowFrom():
            c = a.mapRowToCol( x ) 
            if c: return c
            x -= 1
        else:
            return 0

    def MyMapRight( a, x):
        while x <= a.getRowTo():
            c = a.mapRowToCol( x ) 
            if c: return c
            x += 1
        else:
            return 0

    new_exons = []
    for e in exons:
        
        a = MyMapRight( long_map, e.mPeptideFrom + 1)
        b = MyMapLeft( long_map, e.mPeptideTo + 1)        
        
        ne = e.GetCopy()
        ne.mPeptideFrom, ne.mPeptideTo = a-1, b-1
        new_exons.append( ne )
            
    return new_exons

###################################################################################
def CountMissedBoundaries( cmp_boundaries, reference_boundaries,
                           max_slippage = 9,
                           min_from = 0,
                           max_to = 0):
    """count missed boundaries comparing cmp to ref.
    """

    nmissed = 0
    last_x = None
    ## check if all exon boundaries are ok
    for x in cmp_boundaries:
        if x <= min_from: continue
        if max_to and x >= max_to: continue
        is_ok = 0
        for c in reference_boundaries:
            if abs(x-c) < max_slippage:
                is_ok = 1
                break
        if not is_ok:
            nmissed += 1
                
    return nmissed

###################################################################################
def GetExonsRange( exons,
                   first, last,
                   full = True,
                   min_overlap = 0,
                   min_exon_size = 0 ):
    """get exons in range (first:last) (peptide coordinates).

    Set full to False, if you don't require full coverage.
    """

    new = []
    me = 3 * min_exon_size
    mo = 3 * min_overlap
    for e in exons:
        if e.mPeptideFrom > last or e.mPeptideTo < first:
            continue
        if full and ( e.mPeptideTo > last or e.mPeptideFrom < first):
            continue
        overlap = min(e.mPeptideTo, last) - max(e.mPeptideFrom, first)
        if overlap < mo: continue
        if e.mPeptideTo - e.mPeptideFrom < me: continue
        new.append( e.GetCopy() )
        
    return new

###################################################################################
def ClusterByExonIdentity( exons,
                           max_terminal_num_exons = 3,
                           min_terminal_exon_coverage = 0.0,
                           max_slippage = 0,
                           loglevel = 0):

    """build clusters of transcripts with identical exons.

    The boundaries in the first/last exon can vary.

    Returns two maps map_cluster2transcripts and
    map_transcript2cluster
    """

    ########################################################
    ## build maps
    num_exons = {}
    map_transcript2cluster = {}
    map_cluster2transcripts = {}
    list_of_exons = []
    
    for k, ee in exons.items():
        if k not in map_transcript2cluster:
            map_transcript2cluster[k] = k
            map_cluster2transcripts[k] = [k,]
        num_exons[k] = len(ee)
        
        for e in ee:
            list_of_exons.append( e )

    ########################################################
    ## sort exons
    list_of_exons.sort( lambda x, y: cmp( (x.mSbjctToken, x.mSbjctStrand, x.mGenomeFrom, x.mGenomeTo),
                                          (y.mSbjctToken, y.mSbjctStrand, y.mGenomeFrom, y.mGenomeTo)))

    if loglevel >= 1:
        print "# sorted %i exons" % len(list_of_exons)

    ########################################################
    ## cluster by overlap
    l = list_of_exons[0]
    last_id = l.mQueryToken
    
    for e in list_of_exons[1:]:
        
        identity = False
        coverage = 0
        overlap = 0

        ## check if they are already clustered:
        last_cluster = map_transcript2cluster[last_id]
        if e.mQueryToken in map_cluster2transcripts[last_cluster]:
            l = e
            last_id = e.mQueryToken
            continue
        
        if l.mSbjctStrand == e.mSbjctStrand and \
           l.mSbjctToken == e.mSbjctToken:
            
            ## overlap and coverage of larger exon
            overlap = min(e.mGenomeTo, l.mGenomeTo) - max( e.mGenomeFrom, l.mGenomeFrom)
            coverage = float(overlap) / max( (e.mGenomeTo-e.mGenomeFrom), (l.mGenomeTo-l.mGenomeFrom) )

            left_ok  = abs( l.mGenomeFrom - e.mGenomeFrom ) <= max_slippage
            right_ok = abs( l.mGenomeTo   - e.mGenomeTo ) <= max_slippage
            if overlap > 0:
                ## join if exon boundaries are identical
                if left_ok and right_ok:
                    identity = True
                ## join single exon genes if minimum coverage 
                elif (num_exons[e.mQueryToken] == 1 or num_exons[l.mQueryToken] == 1) and \
                         coverage >= min_terminal_exon_coverage:
                    identity = True
                ## join if one boundary is correct for terminal exons
                ## (for up to three exon genes)
                elif ( left_ok or right_ok ) and \
                      (num_exons[e.mQueryToken] < max_terminal_num_exons and \
                       num_exons[l.mQueryToken] < max_terminal_num_exons) and \
                      coverage >= min_terminal_exon_coverage :
                    identity = True
                ## join, if all exons are overlapping (minimum coverage)
                elif CheckContainedAinB( exons[l.mQueryToken], exons[e.mQueryToken],
                                         min_terminal_exon_coverage, loglevel=loglevel) or \
                     CheckContainedAinB( exons[e.mQueryToken], exons[l.mQueryToken],
                                         min_terminal_exon_coverage, loglevel=loglevel): 
                     identity = True
                    
        if loglevel >= 3:
            print "# ClusterByExonIdentity"
            print "#", identity, overlap, coverage, overlap, num_exons[e.mQueryToken], num_exons[l.mQueryToken]
            print "# l=", str(l)
            print "# e=", str(e)
            
        if identity:
            ## add current cluster to previous cluster.
            this_cluster = map_transcript2cluster[e.mQueryToken]
            map_cluster2transcripts[last_cluster] += map_cluster2transcripts[this_cluster]
            for x in map_cluster2transcripts[this_cluster]:
                map_transcript2cluster[x] = last_cluster
            map_cluster2transcripts[this_cluster] = []
            # print e.mQueryToken, last_cluster, map_cluster2transcripts

        l = e
        last_id = e.mQueryToken

    return map_cluster2transcripts, map_transcript2cluster

##----------------------------------------------------------------------------------------
def ClusterByExonOverlap( exons,
                          min_overlap = 0,
                          min_min_coverage = 0,
                          min_max_coverage = 0,
                          loglevel = 0 ):
    """build clusters of transcripts with overlapping exons.

    Exons need not be identical.

    Returns two maps map_cluster2transcripts and
    map_transcript2cluster
    """
    
    map_transcript2cluster = {}
    map_cluster2transcripts = {}
    list_of_exons = []

    for k, ee in exons.items():
        if k not in map_transcript2cluster:
            map_transcript2cluster[k] = k
            map_cluster2transcripts[k] = [k,]
        
        for e in ee:
            list_of_exons.append( e )
    
    list_of_exons.sort( lambda x, y: cmp( (x.mSbjctToken, x.mSbjctStrand, x.mGenomeFrom),
                                          (y.mSbjctToken, y.mSbjctStrand, y.mGenomeFrom)))

    if loglevel >= 1:
        print "# sorted %i exons" % len(list_of_exons)

    last_id = list_of_exons[0].mQueryToken
    last_from = list_of_exons[0].mGenomeFrom
    last_to = list_of_exons[0].mGenomeTo
    last_token = list_of_exons[0].mSbjctToken
    last_strand = list_of_exons[0].mSbjctStrand

    for e in list_of_exons[1:]:

        last_cluster = map_transcript2cluster[last_id]
        if e.mQueryToken in map_cluster2transcripts[last_cluster]:
            last_id = e.mQueryToken
            last_strand = e.mSbjctStrand
            last_token = e.mSbjctToken
            last_from = e.mGenomeFrom
            last_to = e.mGenomeTo
            continue

        o = min(e.mGenomeTo, last_to) - max(e.mGenomeFrom, last_from)
        o1 = float(o) / (e.mGenomeTo - e.mGenomeFrom)
        o2 = float(o) / (last_to - last_from)

        if last_strand == e.mSbjctStrand and \
               last_token == e.mSbjctToken and \
               o >= min_overlap and \
               min(o1,o2) >= min_min_coverage and \
               max(o1,o2) >= min_max_coverage:
            
            if loglevel >= 3:
                print "# ClusterByExonOverlap %i(%5.2f/%5.2f) between %s and %s: %i-%i with %i-%i" % \
                      ( o, o1, o2,
                        last_id, e.mQueryToken, last_from, last_to, e.mGenomeFrom, e.mGenomeTo)
                sys.stdout.flush()
            
            if e.mQueryToken not in map_cluster2transcripts[last_cluster]:
                this_cluster = map_transcript2cluster[e.mQueryToken]
                map_cluster2transcripts[last_cluster] += map_cluster2transcripts[this_cluster]
                for x in map_cluster2transcripts[this_cluster]:
                    map_transcript2cluster[x] = last_cluster
                map_cluster2transcripts[this_cluster] = []
                    
            last_from = min(e.mGenomeFrom, last_from)
            last_to = max(e.mGenomeTo, last_to)
        else:
            last_id = e.mQueryToken
            last_strand = e.mSbjctStrand
            last_token = e.mSbjctToken
            last_from = e.mGenomeFrom
            last_to = e.mGenomeTo
            
    return map_cluster2transcripts, map_transcript2cluster
    
##----------------------------------------------------------
def CheckOverlap( exons1,
                  exons2,
                  min_overlap = 1):
    """check if exons overlap.

    (does not check chromosome and strand.)
    """

    ## this could be quicker if sorted, but don't want to sort inplace.
    ## there are never that many exons anyway.
    for e1 in exons1:
        for e2 in exons2:
            if min(e1.mGenomeTo, e2.mGenomeTo) - max(e1.mGenomeFrom, e2.mGenomeFrom) >= min_overlap:
                return True
    return False

##----------------------------------------------------------
def CheckCoverage( exons1, exons2,
                   max_terminal_num_exons = 3,
                   min_terminal_exon_coverage = 0.0,
                   max_slippage = 0):
    """check if one set of exons covers the other.
    
    Note: does not check chromosome and strand, just genomic coordinates.
    """

    if len(exons1) > len(exons2):
        ee1, ee2 = exons1, exons2
    else:
        ee1, ee2 = exons2, exons1

    nexons1 = len(ee1)
    nexons2 = len(ee2)

    for e2 in ee2:
        identity = False
        
        for e1 in ee1:
            overlap = min(e1.mGenomeTo, e2.mGenomeTo) - max(e1.mGenomeFrom, e2.mGenomeFrom)

            if overlap <= 0: continue
            
            ## coverage of larger exon
            coverage = float(overlap) / max( (e1.mGenomeTo-e1.mGenomeFrom), (e2.mGenomeTo-e2.mGenomeFrom) )

            if overlap > 0:
                left_ok  = abs( e1.mGenomeFrom - e2.mGenomeFrom ) <= max_slippage
                right_ok = abs( e1.mGenomeTo   - e2.mGenomeTo ) <= max_slippage
                
                if e1.mRank == 0 and e2.mRank == 0:
                    ## for internal exons: the boundaries have to be within slippage distance
                    if left_ok and right_ok:
                        identity = True
                else:
                    ## for terminal exons and less than 4 exons, at least one
                    ## boundary has to be exact
                    if (left_ok or right_ok) and \
                        (nexons1 <= max_terminal_num_exons and nexons2 >= max_terminal_num_exons) and \
                        coverage >= min_terminal_exon_coverage :
                        identity = True

            if identity: break
            
        if not identity:
            return False

    return True

##----------------------------------------------------------
def CheckContainedAinB( exons1, exons2,
                        min_terminal_exon_coverage = 0.0,
                        loglevel = 0 ):
    """check if all exons in exons1 are contained in exons2.

    Note: does not check contig and strand.
    """
    
    for e1 in exons1:
        found = False
        for e2 in exons2:
            overlap = min(e1.mGenomeTo, e2.mGenomeTo) - max(e1.mGenomeFrom, e2.mGenomeFrom)
            ## coverage of smaller exon
            min_length = min((e1.mGenomeTo-e1.mGenomeFrom), (e2.mGenomeTo-e2.mGenomeFrom) )
            if min_length > 0:
                coverage = float(overlap) / min_length
            else:
                if overlap > 0:
                    coverage = 1
                else:
                    coverage = 0
                    
            if coverage >= min_terminal_exon_coverage:
                found = True
                break
        if not found: return False
    return True
    
##----------------------------------------------------------
def CheckCoverageAinB( exons1, exons2,
                       min_terminal_num_exons = 3,
                       min_terminal_exon_coverage = 0.0,
                       max_slippage = 0,
                       loglevel = 0 ):
    """check if exons1 are all in exons2
    
    Note: does not check contig and strand.
    """

    nexons1 = len(exons1)
    nexons2 = len(exons2)

    for e1 in exons1:
        identity = False
        
        for e2 in exons2:
            overlap = min(e1.mGenomeTo, e2.mGenomeTo) - max(e1.mGenomeFrom, e2.mGenomeFrom)
            ## coverage of larger exon
            coverage = float(overlap) / max( (e1.mGenomeTo-e1.mGenomeFrom), (e2.mGenomeTo-e2.mGenomeFrom) )

            left_ok  = abs( e1.mGenomeFrom - e2.mGenomeFrom ) <= max_slippage
            right_ok = abs( e1.mGenomeTo   - e2.mGenomeTo ) <= max_slippage

            if overlap > 0:
                if left_ok and right_ok:
                    identity = True
                    # terminal exons
                elif left_ok or right_ok:
                    if (e1.mRank != 0 or e2.mRank != 0) and \
                           coverage >= min_terminal_exon_coverage and \
                           nexons1 >= min_terminal_num_exons and \
                           nexons2 >= min_terminal_num_exons:
                        identity = True
                        
            if identity: break

        if loglevel >= 2:
            if identity:
                print "# exon found", str(e1), coverage
            else:
                print "# exon not found", str(e1), coverage
            
        if not identity:
            return False

    return True

##----------------------------------------------------------
def GetPeptideLengths( exons ):
    """for all exons get maximum length in coding nucleotides."""

    lengths = {}
    for k, ee in exons.items():
        lengths[k] = max( map( lambda x: x.mPeptideTo, ee) ) / 3 + 1
        
    return lengths
        
##----------------------------------------------------------
def GetGenomeLengths( exons ):
    """for all exons get maximum nucleotide."""

    lengths = {}
    for k, ee in exons.items():
        if ee[0].mGenomeFrom < ee[0].mGenomeTo:
            mi = min( map( lambda x: x.mGenomeFrom, ee) )
            ma = max( map( lambda x: x.mGenomeTo, ee) )        
        else:
            mi = min( map( lambda x: x.mGenomeTo, ee) )
            ma = max( map( lambda x: x.mGenomeFrom, ee) )        
            
        lengths[k] = ma - mi
        
    return lengths
        
##----------------------------------------------------------
def CalculateStats( exons ):
    """calculate some statistics for all exons.

    minimum/maximum intron/exon length, number of exons
    gene length
    """
    stats = {}
    for k, ee in exons.items():
        r = {}
        ee.sort( lambda x,y: cmp(x.mGenomeFrom, y.mGenomeFrom) )
        
        last_to = ee[0].mGenomeTo
        max_i = 0
        min_i = 0
        max_e = ee[0].mGenomeTo - ee[0].mGenomeFrom
        min_e = max_e
        for e in ee[1:]:
            max_e = max( max_e, e.mGenomeTo - e.mGenomeFrom )
            min_e = min( min_e, e.mGenomeTo - e.mGenomeFrom )
            max_i = max( max_i, e.mGenomeFrom - last_to )
            if min_i:
                min_i = min( min_i, e.mGenomeFrom - last_to )
            else:
                min_i = e.mGenomeFrom - last_to
                
        r['NumExons'] = len(ee)
        r['MaxIntronLength'] = max_i
        r['MinIntronLength'] = min_i
        r['MaxExonLength'] = max_e
        r['MinExonLength'] = min_e
        r['GeneLength'] = ee[-1].mGenomeTo - ee[0].mGenomeFrom
        
        stats[k] = r
        
    return stats
        
##----------------------------------------------------------        
def MatchExons( map_a2b, in_exons1, in_exons2, threshold_slipping_boundary = 9):
    """returns a list of overlapping exons (mapped via map_a2b)."""

    cmp_exons = MapExons( in_exons1, map_a2b )
    ref_exons = in_exons2
    
    e,r = 0,0
    matches = []
    
    while e < len(cmp_exons) and r < len(ref_exons):
        
        exon1 = cmp_exons[e]
        exon2 = ref_exons[r]

        ## no overlap
        ## cmp can have 0,0 entries, these are exons that were impossible to align
        if exon1.mPeptideFrom == 0 and exon1.mPeptideTo == 0:
            e += 1
        elif exon1.mPeptideFrom < 0 or exon1.mPeptideTo < 0:
            e += 1
        elif exon2.mPeptideFrom < 0 or exon2.mPeptideTo < 0:
            r += 1
        else:
            overlap  = min(exon1.mPeptideTo,  exon2.mPeptideTo) - max(exon1.mPeptideFrom, exon2.mPeptideFrom)
            if overlap > 0:
                matches.append( (e,r) )

            if exon1.mPeptideTo < exon2.mPeptideTo:
                e += 1
            elif exon2.mPeptideTo < exon1.mPeptideTo:
                r += 1
            else:
                e += 1
                r += 1

    return matches

###################################################################################
    
