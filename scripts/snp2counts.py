################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: snp2counts.py 2872 2010-03-03 10:21:13Z andreas $
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
"""
snp2counts.py - count number SNPs in geneset
============================================

:Author: Andreas Heger
:Release: $Id: snp2counts.py 2872 2010-03-03 10:21:13Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

read a list of genomic point locations (SNPs) and count the number
of SNPs falling in pre-defined windows.

The windows are given in gtf format. 

.. note::
   The script will be able to count snps in disjoint segments
   using the gene_id field in gtf format. It will not check 
   if these segments are non-overlapping. 

   In case of a gene set, make sure to first flatten the gene set by combining
   all transcript/exons per gene.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

.. note::
   The script currently uses ``variant`` in two meanings:
   
   1. a variable site (SNP/INDEL)

   2. a transcript variant (a transcript sequence that differs from the wild type)

   I have started calling the latter ``allele``, though it is not
   consistent across the whole script. However, the output is consistent and calls
   the former ``variant_site`` and the latter ``allele``.

""" 

import os
import sys
import re
import optparse
import collections

import numpy
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import pysam
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IndexedGenome as IndexedGenome
import CGAT.Genomics as Genomics
import CGAT.GTF as GTF
import alignlib

CdsResult = collections.namedtuple( 'CdsResult',
                                    '''strand, start, end,
                                    exon_id, exon_start, exon_end,
                                    cds_start, cds_end, cds_phase,
                                    intron_id, intron_start, intron_end,
                                    prev_exon_end, next_exon_start,
                                    cds_seq, cds_seq_start, cds_seq_end,
                                    nc_seq, nc_start, nc_end, exon_skipping
                                    ''')

Variant = collections.namedtuple( 'Variant',
                                  'code,sequence,start,end')

SpliceEffect = collections.namedtuple( 'SpliceEffect',
                                       '''exon_id, orig_seq, variant_seq''' )

SpliceChange = collections.namedtuple( 'SpliceChange',
                                       '''exon_id,
                                       is_frameshift,
                                       orig_name, orig_seq5, orig_seq3,
                                       variant_name, variant_seq5, variant_seq3''')

CdsEffect = collections.namedtuple( 'CdsEffect',
                                    '''exon_id, orig_seq, variant_seq,
                                    codon_orig_seq, codon_variant_seq''' )

CdsVariant = collections.namedtuple( 'CdsVariant',
                                     '''transcript_id,
                                     cds_start, cds_end,
                                     code,
                                     reference_seq,variant_seq,
                                     is_homozygous''' )

SpliceVariant = collections.namedtuple( 'SpliceVariant',
                                        '''transcript_id,
                                        intron_id,
                                        nc_start, nc_end,
                                        code,
                                        reference_seq,variant_seq,
                                        is_homozygous''' )

TranscriptVariant = collections.namedtuple( 'TransciptVariant',
                                            '''cds_variants,splice_variants''')

# counters for the cumulative effect of variants on the translated sequence
TranslationEffect = collections.namedtuple( 'TranslationEffect','''
                   ncodons,
                   ninserted_bases,
                   ninserted_codons,
                   ndeleted_bases,
                   ndeleted_codons,
                   nincomplete_codons,
                   noframe_codons,
                   nwrong_frames,
                   ncorrected_frames,
                   first_stop,
                   nstops,
                   nunaffected_codons,
                   nsynonymous_codons,
                   nnonsynonymous_codons,
                   nstop_codons''' )

SplicingEffect = collections.namedtuple( 'SplicingEffect',
                                         '''nintrons,
                                         ncanonical,
                                         nframeshifts,
                                         nnoncanonical,
                                         nunchanged_frames,
                                         ncorrected_frames,
                                         nuncorrected_frames,
                                         nunchanged,
                                         nsynonymous,
                                         nnonsynonymous,
                                         ndisrupted,
                                         nnovel,
                                         nnunknown,
                                         ninserted_codons,
                                         codes,
                                         last_exon''' )


def iterateOverFrames( variant_seq ):
    '''return tuples of segments within/without
    frame.

    Yields only coordinates in frame (multiples of 3)
    Everything that is out-of-frame is yielded together.
    '''

    frame_at_start = len("".join(variant_seq[:3]) ) % 3
    frame_at_end = 0
    start = 0
    for x in range( 3, len(variant_seq ), 3 ):

        var_codon = "".join( variant_seq[x:x+3] ).upper()
        lvar = len(var_codon)

        frame_at_end = (frame_at_start + lvar) % 3
        # print x, frame_at_start, frame_at_end, start

        # check for frame change
        if frame_at_end != frame_at_start:
            if frame_at_start == 0:
                # exclude current codon
                yield( (True, start, x) )
                start = x
            elif frame_at_start != 0 and frame_at_end == 0:
                # include current codon
                yield( (False, start, x+3) )
                start = x + 3
            else:
                # nothing to be done if frame changes
                # between out-of-frame frames
                pass
        frame_at_start = frame_at_end

    if start != len(variant_seq):
        yield( (frame_at_end == 0, start,len(variant_seq)) )

def countEffectsOnTranscript( var_seq, ref_seq, 
                              is_seleno = False ):

    '''count effects on transcript.

    *var_seq* is a list of characters according to a known cds and
    thus determines the reference frame.

    Insertions contain more than one character at a position, deletions
    are empty.

    The function returns a namedtuple of type TranscriptEffect. Counts are in terms of base/codons
    in the reference sequence.

    Counting will continue after a stop-codon is encountered.

    Note that codons inserted within a codon do not count as a frame shift. Instead
    these will be recorded as an inserted codon.

    ncodons
       number of codons in transcript
    ninserted_bases
       number of inserted bases
    ninserted_codons
       number of fully inserted codons
    ndeleted_bases
       number of deleted bases
    nincomplete_codons
       number of incomplete codons at the end
    ndeleted_codons
       number of fully deleted codons
    noframe_codons
       number of codons that are out-of-frame. This will include all codons where
       at least one base is out-of-frame. In case of an in+del, the codon will
       still be out-of-frame.
    nwrong_frames
       number of times the sequence gets out of frame
    ncorrected_frames
       number of times the frame is recovered
    nstops
       number of stop codons in translation
    nsynonymous_codons:
       number of reference codons that have synonymous mutations
    nnonsynonymous_codons:
       number of referenec codons that have non-synonymous mutations
    nstop_codons
       number of reference codons that now encode for a stop
    nunaffected_codons
       number of reference codons that are still the same
    first_stop
       codon position of first stop codon in variant sequence
    '''

    assert len(var_seq) == len(ref_seq)

    # values to fill
    ncodons = 0
    ninserted_bases, ndeleted_bases = 0, 0
    ninserted_codons, ndeleted_codons = 0, 0
    nincomplete_codons = 0
    noframe_bases, noframe_codons = 0, 0
    nwrong_frames, ncorrected_frames = 0, 0
    nsynonymous_codons, nnonsynonymous_codons, nunaffected_codons = 0, 0, 0
    nstop_codons = 0
    last_exon_start = 0

    # build sequences
    var_seq_na = "".join( var_seq ).upper()
    ref_seq_na = "".join( ref_seq ).upper()

    lrefseq = len(ref_seq)
    lvarseq = len(var_seq_na)

    ncodons = lrefseq // 3

    # truncate incomplete base at end of reference sequence
    if lvarseq % 3 != 0: var_seq_na = var_seq_na[:-(lvarseq%3)]
    var_seq_aa = Genomics.translate( var_seq_na,
                                     is_seleno = is_seleno )

    # check protein coding sequence for the first stop
    nstops = 0
    first_stop = len(var_seq_aa)
    ntruncated_codons_stop = 0
    
    for pos, c in enumerate( var_seq_aa ):
        if c == "X":
            nstops += 1
            first_stop = min( pos, first_stop )

    # start position for out-of-frame region
    var_pos = 0

    map_ref2var = alignlib.makeAlignmentVector()
    alignator = alignlib.makeAlignatorDPFull( alignlib.ALIGNMENT_GLOBAL,
                                              -10.0,
                                              -2.0 )

    was_in_frame = True
    for in_frame, start, end in iterateOverFrames( var_seq ):
        varseq = "".join(var_seq[start:end]).upper()
        refseq = "".join(ref_seq[start:end]).upper()

        # print in_frame, start, end
        if in_frame:
            for x in range(start, end,3):
                # ignore incomplete codons at the end:
                if x+3 > lrefseq: break

                var_codon = "".join(var_seq[x:x+3]).upper()

                assert len(var_codon) % 3 ==  0
                ref_codon = "".join(ref_seq[x:x+3]).upper()

                assert len(ref_codon) == 3
                    
                d = len(var_codon) - 3
                
                for y in var_seq[x:x+3]:
                    if y == "": ndeleted_bases += 1
                    if len(y) > 1: ninserted_bases += len(y) - 1

                if var_codon == "":
                    ndeleted_codons -= d // 3
                elif len(var_codon) > len(ref_codon):
                    # deal with in-frame inserted codons
                    ninserted_codons += d // 3
                    nunaffected_codons += 1
                elif var_codon == ref_codon:
                    nunaffected_codons += 1
                else:
                    var_aa = Genomics.translate( var_codon )
                    ref_aa = Genomics.translate( ref_codon )
                    if var_aa == "X":
                        nstop_codons += 1
                    elif ref_aa == var_aa:
                        nsynonymous_codons += 1
                    else:
                        nnonsynonymous_codons += 1
                var_pos += len(var_codon)
                
        else:
            # count deletions/insertions in the variant
            for x in range(start, end, 3):
                var_codon = "".join(var_seq[x:x+3]).upper()
                # count insertion and deletion separately to avoid them
                # compensating
                for y in var_seq[x:x+3]:
                    if y == "": ndeleted_bases += 1
                    if len(y) > 1:
                        ninserted_bases += len(y) - 1
                        ninserted_codons += (len(y) - 1) // 3
                    
                # completely deleted codons
                if var_codon == "":
                    ndeleted_codons += 1
                else:
                    noframe_codons += 1

            # count effect on protein coding sequence
            var_frag_aa = Genomics.translate( varseq )
            ref_frag_aa = Genomics.translate( refseq )

            # count effect on protein coding sequence
            var_s = alignlib.makeSequence(var_frag_aa)
            ref_s = alignlib.makeSequence(ref_frag_aa)

            diff_length = abs(len(ref_frag_aa) - len(var_frag_aa))
            # very heuristic - might lead strange effects
            alignment_band = max(10, diff_length * 2)
            iterator = alignlib.makeIterator2DBanded( -alignment_band, +alignment_band )
            alignlib.setDefaultIterator2D( iterator )

            E.debug( "alignment: reference(%i) with variant(%i) (diff=%i) within diagonals %i and %i" % \
                     (len(ref_frag_aa), len(var_frag_aa), diff_length, -alignment_band, alignment_band ))
            
            alignator.align(  map_ref2var, ref_s, var_s )
            # print alignlib.AlignmentFormatExplicit( map_ref2var, ref_s, var_s )
            for x, ref_aa in enumerate( ref_frag_aa ):
                p = map_ref2var.mapRowToCol(x)
                if p < 0: continue
                var_aa = var_frag_aa[p]

                if var_aa == ref_aa:
                    nsynonymous_codons += 1
                else:
                    nnonsynonymous_codons += 1
                
            nwrong_frames += 1
            ncorrected_frames += 1
            
        was_in_frame = in_frame
        # if var_pos > first_stop * 3: break

    if lvarseq % 3 != 0:
        nincomplete_codons += 1
        
    # reduce corrected frames by one if we do not end on frame        
    if not was_in_frame and lvarseq % 3 != 0: ncorrected_frames -= 1

    return TranslationEffect._make( (ncodons,
                                     ninserted_bases,
                                     ninserted_codons,
                                     ndeleted_bases,
                                     ndeleted_codons,
                                     nincomplete_codons, 
                                     noframe_codons,
                                     nwrong_frames,
                                     ncorrected_frames,
                                     first_stop,
                                     nstops,
                                     nunaffected_codons,
                                     nsynonymous_codons,
                                     nnonsynonymous_codons,
                                     nstop_codons) )


def getCDSPosition( exons, start, end, fasta = None, lcontig = None ):
    '''return cds information for a (positive strand) genomic position.

    exons is a list of exons in GTF format.
    start, end: are the coordinates of the variant in forward strand coordinates.

    if the first exon is not in frame, cds_start, cds_end will be
    not multiples of 3, but that is correct, as cds_start and cds_end
    start counting from 0.

    If a region spans a whole intron, the region will be treated as
    a single coding sequence variant. Such deletions usually concern
    short frame-shifting introns.
    
    '''
    strand = exons[0].strand
    contig = exons[0].contig
    is_positive_strand = Genomics.IsPositiveStrand(strand)
    
    if is_positive_strand:
        coordinates = [ (e.start, e.end, int(e.frame)) for e in exons ]
    else:
        if not fasta and not lcontig:
            raise ValueError( "no fasta or lcontig option given for a negative strand transcript" )

        if fasta:
            # convert all to negative strand coordinates
            lcontig = fasta.getLength( contig )
            
        start, end = lcontig - end, lcontig - start
        coordinates = [ (lcontig- e.end, lcontig - e.start, int(e.frame)) for e in exons ]
    coordinates.sort()

    # phase is the complement to frame (i.e. position in codon, not base to next codon)
    cds_start, cds_end, cds_phase = None, None, None
    # coordinates for previous/next exons for snps spanning splice sites
    prev_exon_end, next_exon_start = None,None
    # intron positions
    intron_start, intron_end = None, None
    # start, end of feature within coding sequence
    # sequence is a sequence of all codons that cover the feature
    cds_seq_start, cds_seq_end, cds_seq = None, None, None

    # flag to denote exon skipping
    exon_skipping = False
    
    # start, end of feature within non-coding sequence
    nc_seq_start, nc_seq_end, nc_seq = None, None, None

    exon_id = None
    
    # empty result in case of no overlap
    if end <= coordinates[0][0] or start >= coordinates[-1][1]:
        return None

    intron_id, exon_id = None, 0
    nexons = len(coordinates)
    last_exon = nexons-1
    start_phase = (3 - coordinates[0][2]) % 3
    # start within frame
    cds_start = start_phase

    # find exon overlapping the region or exon immediately after it
    while exon_id < len(coordinates):

        exon_start, exon_end, exon_frame = coordinates[exon_id]
        if start < exon_end: break
        cds_start += exon_end - exon_start
        exon_id += 1

    if end <= exon_start:
        # overlap with intron only
        cds_start, cds_end = None, None
        if exon_id > 0:
            intron_start, intron_end = coordinates[exon_id-1][1], exon_start
            nc_seq_start, nc_seq_end = start, end
            intron_id = exon_id - 1
    else:
        # overlap with exon
        #
        # position of first complete codon in this feature:
        first_full_codon_start = exon_start + exon_frame
        
        # correction of frame at end of exon to due frame-shifting codons
        frame_correction = False
        
        # special treatment if region spans the complete intron
        if exon_id < last_exon and end > coordinates[exon_id+1][0]:
            if end > coordinates[exon_id+1][1]:
                raise ValueError( "can not deal with variants spanning multiple exons")
            # simply increase the current exon
            exon_end = coordinates[exon_id+1][1]
            # in order to adjust the frame, add the intron towards the exon
            frame_correction = True

        cds_x, cds_y = max(start, exon_start), min(end, exon_end)
        cds_start += cds_x - exon_start
        cds_end = cds_start + cds_y - cds_x
        cds_phase = (3 - exon_frame + cds_x - exon_start) % 3

        # print "exon_id=", exon_id, "start=", start, "end=", end, \
        #     "codon_start", cds_x, "codon_end", cds_y, \
        #     "cds_phase=", cds_phase, "cds_start", cds_start, "cds_end", cds_end, \
        #           "exon_start=",exon_start, "exon_end=", exon_end, "exon_frame=", exon_frame

        # incomplete bases in a codon at 3' end of this feature:
        # = frame of next feature
        # (3 - ( exon_end - exon_start - exon_frame) % 3 ) % 3
        if exon_id < last_exon:
            last_full_codon_end = exon_end - (3 - coordinates[exon_id+1][2]) % 3
            if not frame_correction:
                assert (3 - ( exon_end - exon_start - exon_frame + frame_correction ) % 3 ) % 3 == coordinates[exon_id+1][2], \
                       "frame mismatch between exons %i and %i" % (exon_id, exon_id+1)
        else:
            last_full_codon_end = 0

        # link to previous next exons in the case of split exons
        is_split_start = exon_start <= cds_x < first_full_codon_start and exon_id > 0
        is_split_end = last_full_codon_end < cds_y <= exon_end and exon_id < last_exon

        if is_split_start:
            prev_exon_end = coordinates[exon_id-1][1]
        if is_split_end:
            next_exon_start = coordinates[exon_id+1][0]
            next_frame = coordinates[exon_id+1][2]
            
        # sort out the sequence
        cds_seq = []

        # position of variant in cds_sequence
        cds_seq_start, cds_seq_end = cds_phase, cds_y - cds_x + cds_phase

        if fasta:

            # link to previous/next exons in the case of split codons
            if is_split_start:
                codon_start = prev_exon_end - (3 - exon_frame) % 3
                codon_end = prev_exon_end
                cds_seq.insert( 0, fasta.getSequence( contig,
                                                      strand,
                                                      codon_start,
                                                      codon_end))

            codon_start = cds_x - cds_phase
            codon_end = cds_y + (3 - (cds_end % 3) ) % 3

            # deal with incomplete codon at start
            if codon_start < exon_start and exon_id == 0:
                assert exon_frame != 0
                cds_seq.extend( list("X" * (exon_start - codon_start) ))
                
            # print "exon_id=", exon_id, "start=", start, "end=", end, "codon_start", codon_start, "codon_end", codon_end, "cdsx", cds_x, "cdsy", cds_y, cds_phase, "cds_start", cds_start, "cds_end", cds_end, \
            #      "exon_start=",exon_start, "exon_end=", exon_end, "start_phase=", start_phase, "first_start", first_full_codon_start, "last_end", last_full_codon_end, \
            #      "split_start", is_split_start, "split_end", is_split_end
            
            cds_seq.extend( list( fasta.getSequence( contig,
                                                     strand,
                                                     max(exon_start, codon_start),
                                                     min(exon_end, codon_end) ) ) )

            # fill up, if last codon is incomplete
            if codon_end > exon_end and exon_id == last_exon:
                cds_seq.extend( list("X" * (codon_end - exon_end ) ) )

            if is_split_end:
                cds_seq.append( fasta.getSequence( contig,
                                                   strand,
                                                   next_exon_start,
                                                   next_exon_start + next_frame 
                                                   ))
                    
        cds_seq = "".join( cds_seq )
        lnoncoding = (end - start) - (cds_y - cds_x)

        if start <= exon_start and end >= exon_end:
            # special treatment if region spans the complete exon
            if exon_id < nexons - 1:
                nc_seq_start, nc_seq_end = start, end
                intron_start = prev_exon_end
                intron_end = next_exon_start
                intron_id  = exon_id - 1
            else:
                # unless it is last exon - truncate, but only
                # it it extends into intron
                if start < exon_start:
                    intron_start, inron_end = prev_exon_end, exon_start
                    nc_seq_start, nc_seq_end = start, exon_start
                    intron_id = exon_id - 1
                    
            exon_skipping = True

        elif start < exon_start and exon_id > 0:
            # disrupted intronic sequence
            intron_start, intron_end = coordinates[exon_id-1][1], exon_start
            nc_seq_start, nc_seq_end = exon_start - lnoncoding , exon_start
            intron_id = exon_id - 1

        elif end > exon_end and exon_id < nexons-1:
            # disrupted intronic sequence
            intron_start, intron_end = exon_end, coordinates[exon_id+1][0]
            nc_seq_start, nc_seq_end = exon_end, exon_end + lnoncoding
            intron_id = exon_id

    if fasta and nc_seq_start != None:
        nc_seq = fasta.getSequence( contig, strand, nc_seq_start, nc_seq_end )

    # subtract starting frame
    if cds_start != None:
        cds_start -= start_phase
        cds_end -= start_phase

    return CdsResult._make( (strand, start, end,
                             exon_id, exon_start, exon_end,
                             cds_start, cds_end, cds_phase,
                             intron_id, intron_start, intron_end,
                             prev_exon_end, next_exon_start,
                             cds_seq, cds_seq_start, cds_seq_end,
                             nc_seq, nc_seq_start, nc_seq_end,
                             exon_skipping ))

class Counter(object):
    '''annotator for single bases in the genome.'''

    mHeader = ()

    def __init__(self, fasta = None, pattern = "%s", *args, **kwargs):
        self.mFasta = fasta
        self.mFilenamePattern = pattern
        
    def __str__( self ):
        return ""

    def getHeader( self ):
        '''return header'''
        return "\t".join(self.mHeader)

class CounterGenes( Counter ):
    '''count SNPs per gene that it overlaps with.'''

    mHeader = [ "exons_%s" % x for x in ( "ntranscripts", "nused", "pos" )]

    def __init__(self, filename_exons, *args, **kwargs):
        Counter.__init__(self, *args,**kwargs)

        exons = IndexedGenome.IndexedGenome()
        nexons = 0

        inf = IOTools.openFile( filename_exons, "r") 
        for g in GTF.iterator( inf ):
            exons.add( g.contig, g.start, g.end, g )
            nexons += 1
        inf.close()

        self.mExons = exons
        E.info( "indexed %i exons on %i contigs" % (nexons, len(exons) ) )

        # create counter
        self.mCounts = collections.defaultdict( int )

    def update( self, snp ):
        '''update with snp.'''

        exons = list(self.mExons.get( snp.chromosome, snp.pos, snp.pos+1))

        if exons:
            for start, end, gtf in exons:
                self.mCounts[gtf.gene_id] += 1

    def writeTable( self, outfile ):

        outfile.write( "gene_id\tnsnps\n" )
        for key in sorted(self.mCounts.keys()):
            outfile.write( "\t".join( (key, str(self.mCounts[key] ))) + "\n" )

class CounterTranscripts( Counter ):
    '''count SNPs per transcripts that it overlaps with.

        Variants are not phased, so is not always clear which of the two allelles of a transcript
        is affected. Thus, the following heuristic is adopted:
        
        1 Only homozygous variants:   locus flagged as homozygous. Both alleles are assumed to be the same and
                                      different from the wild type.
                                      
        2 Only heterozygous variants: locus flagged as heterozygous. One allele is assumed to be the wild type,
                                      the other one is a variant.
                                      
        3 Mixture of homo- and heterozygous variants: locus flagged as mixture. A mixed allele is constructed with
                                      all variants.

    Columns

    transcript_id
       the transcript_id
    cds_len
       length of cds in bases
    ncodons
       number of codons in wild type
    last_exon_start
       start (cds coordinates) of last exon (useful for detecting nonsense-mediated decay)
    max_variants
       maximum number of variants per site
    nvariant_sites
       number of variable sites within the ntranscript
    genotype
       the genotype
    nalleles
       number of variants (1 = either homozygote variant or heterozygote variant/wild type)
    stop_min
       number of codons truncated either due to disrupted splice signal and/or stop codon.
       This is the minimum between two transcripts. If the wildtype is still present,
       this value will be 0.
    stop_max
       number of codons truncated either due to disrupted splice signal or stop codon.
       This is the maximum between two transcripts.

    Columns are prefixed with ``cds_`` and ``splice_`` for cds and splice variants,
    respectively. Without prefix, it refers to the effecs of cds and splice variants
    combined.
    '''

    # outfile.write( 
    mHeader = [ "transcript_id",
                "cds_len", 
                "ncodons",
                "last_exon_start",
                "cds_max_variants","cds_nvariant_sites","cds_genotype","cds_nalleles",
                "cds_stop_min", "cds_stop_max",
                "splice_max_variants","splice_nvariant_sites","splice_genotype","splice_nalleles",
                "splice_stop_min", "splice_stop_max",
                "max_vars","nvariant_sites","genotype","nalleles",
                "stop_min", "stop_max" ]
        
    # add this area to check for overlap with splice signals
    # This should be larger than the longest deletion.
    mSize = 500

    # introns smaller than this size are considered to be frameshifts
    mMinIntronSize = 5

    def __init__(self, filename_exons, seleno, *args, **kwargs):
        Counter.__init__(self, *args,**kwargs)

        transcripts = IndexedGenome.IndexedGenome()
        self.mExons = {}
        nexons = 0
        ntranscripts = 0
        inf = IOTools.openFile( filename_exons, "r")
        for gtfs in GTF.transcript_iterator( GTF.iterator( inf )):
            start, end = min( [x.start for x in gtfs]), max([x.end for x in gtfs ] )
            transcripts.add( gtfs[0].contig, start, end, gtfs )
            nexons += len(gtfs)
            ntranscripts += 1
            self.mExons[gtfs[0].transcript_id] = gtfs
        inf.close()

        self.mTranscripts = transcripts
        self.mSeleno = seleno
        
        E.info( "indexed %i transcripts and %i exons on %i contigs" % (ntranscripts, nexons, len(transcripts) ) )
        E.info( "received %i selenoprotein transcripts" % (len(self.mSeleno )) )

        # create counter
        self.mCounts = collections.defaultdict( int )

        self.mOutfileIntron = IOTools.openFile( self.mFilenamePattern % "intron", "w" )
        self.mOutfileIntron.write("transcript_id\tcontig\tsnp_position\tvariant_type\tvariant_code\tvariant_seq\texon_id\tnexon\tcode\torig_name\torig_seq5\torig_seq3\tvariant_name\tvariant_seq5\tvariant_seq3\tintron_start\tintron_end\tstrand\tnc_start\tnc_end\n")

        self.mOutfileCds = IOTools.openFile( self.mFilenamePattern % "cds", "w" )
        self.mOutfileCds.write( "\t".join( (
                    "transcript_id",
                    "contig",
                    "snp_position", 
                    "reference", 
                    "variant_type", 
                    "variant_code", 
                    "variant_bases", 
                    "exon_id", 
                    "nexons", 
                    "code", 
                    "orig_seq", 
                    "orig_na", 
                    "orig_codons", 
                    "variant_seq", 
                    "variant_na", 
                    "variant_codons", 
                    "cds_phase", 
                    "cds_start", 
                    "cds_end", 
                    "cds_len") ) + "\n" )

        self.mOutfileTranscripts = IOTools.openFile( self.mFilenamePattern % "translation", "w" )
        self.mOutfileTranscripts.write( "transcript_id\tvariant_id\tlast_exon_start\t%s\tseq_na\tseq_aa\n" % "\t".join(TranslationEffect._fields) )

        self.mOutfileSplicing = IOTools.openFile( self.mFilenamePattern % "splicing", "w" )
        self.mOutfileSplicing.write("transcript_id\tvariant_id\t%s\n" % "\t".join(SplicingEffect._fields) )

        self.mTranscriptVariants = {}
        
    def getVariantRange( self, snp ):
        '''return effective range of a variant.
        
        The effective range is a single base in case of a SNP.
        It is two bases if it is an insertion in case it is a
        coding SNP.
        Deletions are as large as the deletion.
        '''

        contig = snp.chromosome
        lcontig = self.mFasta.getLength( contig )

        reference_base = snp.reference_base

        start, end = lcontig, 0
        
        # process according to variant type
        # indels need to be treated differently from SNPs as
        # they have larger effects
        if reference_base == "*":
            variants = snp.genotype.split("/")
            for variant in variants:
                if variant[0] == "*":
                    continue
                elif variant[0] == "+":
                    start = min( start, snp.pos )
                    end = max( end, snp.pos + 2 )
                elif variant[0] == "-":
                    # deletions are after the base denoted by snp.pos
                    start = min( start, snp.pos + 1 )
                    # pos + 1 + len(var) - 1 = pos + len(var)
                    end = max( end, snp.pos + len(variant) )
                else:
                    raise ValueError( "unknown variant sign '%s'" % variant[0] )
        else:
            # a single base SNP
            start = min( start, snp.pos )
            end = max( end, snp.pos + 1 )

        start, end = max(0,start),min(end,lcontig)
        if start == end:
            return None, None
        else:
            return start, end

    def getSequence( self, snp, r, variant):
        '''return sequence of snp taking into account strandedness of transcript.'''
        
        contig = snp.chromosome
        
        # collect sequences (resolving strandedness)
        reference_base = snp.reference_base

        if reference_base != "*": 
            variant_bases = Genomics.resolveAmbiguousNA( variant.sequence )
            assert len(variant_bases) == 1
        else:
            variant_bases = []

        variant_seq = variant.sequence

        if not Genomics.IsPositiveStrand( r.strand) :
            variant_seq = Genomics.complement( variant_seq )
            variant_bases = [ Genomics.complement( base ) for base in variant_bases ]
            reference_base = Genomics.complement( reference_base )

        return reference_base, variant_seq, variant_bases

    def collectSplicingEffects( self, snp, r, variant, reference_base, variant_seq, variant_bases ):
        '''compute effects of a variant on a transcript.

        The effects are independent of any other variants.

        return a list of splicing effects.
        '''
        intron_effects, intron_changes = [], []
        # collect splicing effects only
        if r.nc_start == None: return intron_effects, intron_changes

        contig = snp.chromosome

        lvariant = len(variant_seq)

        intron_seq = self.mFasta.getSequence( contig, r.strand, r.intron_start, r.intron_end ).upper()
        is_frameshift = len(intron_seq) < self.mMinIntronSize
            
        intron_name, intron_seq5, intron_seq3 = Genomics.GetIntronType( intron_seq )
        variant_introns = []

        if (r.nc_start - r.intron_start) >= len(intron_seq5) and (r.intron_end - r.nc_end) >= len(intron_seq3):
            # intronic variant - ignore if not actually overlapping with splice site
           pass
        else:
            E.debug( "cds=%s, variant=%s" % (str(r), str(snp)))
            variant_intron_seq = list(intron_seq)
            x,y = r.nc_start-r.intron_start, r.nc_end-r.intron_start

            if variant.code == "=":
                # add SNP
                assert y-x == 1, "expect only single base substitutions"
                if intron_seq[x:y] != reference_base:
                    raise ValueError( "expected=%s, got=%s:%s:%s, snp=%s, cds=%s" %\
                                      (reference_base,
                                       intron_seq[x-3:x],
                                       intron_seq[x:y],
                                       intron_seq[y:y+3],
                                       str(snp), str(r)) )

                # record multiple substitutions
                for base in variant_bases:
                    if base != reference_base:
                        variant_intron_seq[x:y] = base
                        variant_introns.append( "".join( variant_intron_seq) )

                        intron_effects.append( SpliceEffect._make( (r.intron_id,
                                                                    reference_base,
                                                                    base ) ) )

            elif variant.code == "+":
                # add insertion
                # If the insertion is at an intron/exon boundary
                # y -x = 1. In this case attribute this to a
                # coding sequence change and ignore
                if y - x == 2:
                    # python inserts before the index
                    variant_intron_seq[y:y] = list(variant_seq)
                    variant_introns.append( "".join( variant_intron_seq) )
                    intron_effects.append( SpliceEffect._make( (r.intron_id,
                                                                "",
                                                                variant_seq ) ) )
                else:
                    if y - x != 1:
                        raise ValueError( "expected an insert of length 1 or 2, got %i for %s" % (y-x, str(snp)) )
                                            
            elif variant.code == "-":
                # add deletion
                if x == 0 and y == r.intron_end - r.intron_start:
                    # deletion covers full length of intron
                    if r.intron_id < r.exon_id:
                        # truncate from start if intron preceding exon
                        xx, yy = 0, y-x
                    else:
                        # truncate from end if intron succceding exon
                        xx, yy = lvariant - (y-x), lvariant
                elif x == 0: 
                    # deletion at 3' end of intron: truncate from the end
                    xx, yy = lvariant - (y-x), lvariant
                else: xx, yy = 0, y-x

                if intron_seq[x:y] != variant_seq[xx:yy]:
                    raise ValueError( "expected=%s, got=%s:%s:%s, %i:%i, %i:%i, snp=%s, cds=%s" %\
                                      (variant_seq[xx:yy],
                                       intron_seq[x-3:x],
                                       intron_seq[x:y],
                                       intron_seq[y:y+3],
                                       x,y,
                                       xx, yy,
                                       str(snp), str(r)) )

                intron_effects.append( SpliceEffect._make( (r.intron_id,
                                                            variant_intron_seq[x:y],
                                                            "" ) ) )

                del variant_intron_seq[x:y]
                variant_introns.append( "".join( variant_intron_seq) )

        for variant_intron_seq in variant_introns:
            variant_intron_name, variant_intron_seq5, variant_intron_seq3 = Genomics.GetIntronType( variant_intron_seq )

            # if intron is a frameshift, the full intron seq is returned
            #if is_frameshift: reference_seq, variant_seq = intron_seq, variant_inseq
            intron_changes.append( SpliceChange._make( (r.exon_id-1,
                                                        is_frameshift,
                                                        intron_name, intron_seq5, intron_seq3,
                                                        variant_intron_name, variant_intron_seq5, variant_intron_seq3) ) )

        return intron_effects, intron_changes

    def collectCodingEffects( self, snp, r, variant, reference_base, variant_seq, variant_bases ):
        '''compute effects of a variant on a transcript.

        The effects are independent of any other variants.

        return a list of cds effects
        '''

        coding_effects = []
        ## process coding effects, return empty if none
        if r.cds_start == None: return coding_effects
        
        contig = snp.chromosome
        lvariant = len(variant_seq)
                
        cds_seq = r.cds_seq.upper()
        variant_cds_seq = list(cds_seq)
        x,y = r.cds_seq_start, r.cds_seq_end

        if len(cds_seq) % 3 != 0:
            raise ValueError( "expected codon sequence, got=%s (%i), %s:%s:%s, %i:%i, snp=%s, cds=%s" %\
                              (cds_seq,
                               len(cds_seq),
                               cds_seq[:x],
                               cds_seq[x:y],
                               cds_seq[y:],
                               x,y,
                               str(snp), str(r)) )

        if variant.code == "=":
            # process substitution
            assert y-x == 1, "expect only single base substitutions"
            if cds_seq[x:y] != reference_base:
                raise ValueError( "expected=%s, got=%s:%s:%s, %i:%i, snp=%s, cds=%s" %\
                                  (reference_base,
                                   cds_seq[:x],
                                   cds_seq[x:y],
                                   cds_seq[y:],
                                   x,y,
                                   str(snp), str(r)) )

            # record multiple substitutions
            for base in variant_bases:
                if base != reference_base:
                    variant_cds_seq[x] = base
                    coding_effects.append( CdsEffect._make( (r.exon_id,
                                                             reference_base,
                                                             base,
                                                             cds_seq,
                                                             "".join( variant_cds_seq ),
                                                             )))

        elif variant.code == "+":
            # add insertion - python inserts before index
            variant_cds_seq[y:y] = variant_seq
            coding_effects.append( CdsEffect._make( (r.exon_id,
                                                     "",
                                                     variant_seq,
                                                     cds_seq,
                                                     "".join( variant_cds_seq) )))

        elif variant.code == "-":
            # add deletion
            if r.exon_skipping:
                xx, yy = r.exon_start - r.nc_start, lvariant - (r.nc_end - r.exon_end)
            elif r.nc_start != None:
                # deletion at exon boundary
                if r.intron_id < r.exon_id:
                    # deletion at 5' end of exon, take only 3' bases of variant
                    xx, yy = lvariant - (y-x), lvariant
                else:
                    # deletion at 3' end of exon, take only 5' bases of variant
                    xx, yy = 0, y - x                  
                # removed the following condition: "and r.nc_start != r.intron_start:"
                # deletion at 3' end of intron boundary - delete last bases
                # xx, yy = lvariant - (y-x), lvariant
            elif r.cds_start == 0:
                # deletion at first codon - take only 3' bases of variant
                xx, yy = lvariant - (y-x), lvariant
            else:
                # deletion after - delete last bases
                xx, yy = 0, y-x
                
            if cds_seq[x:y] != variant_seq[xx:yy]:
                raise ValueError( "expected=%s, got=%s:%s:%s, %i:%i, %i:%i, snp=%s, cds=%s" %\
                                  (variant_seq[xx:yy],
                                   cds_seq[:x],
                                   cds_seq[x:y],
                                   cds_seq[y:],
                                   x,y,
                                   xx,yy,
                                   str(snp), str(r)) )

            del variant_cds_seq[x:y]
            coding_effects.append( CdsEffect._make( (r.exon_id,
                                                     cds_seq[x:y],
                                                     "",
                                                     cds_seq,
                                                     "".join( variant_cds_seq )) ))

        return coding_effects

                    
    def update( self, snp ):
        '''update with snp.'''

        # get effective range of snp
        snp_start, snp_end = self.getVariantRange( snp )

        # ignore snps that are out-of-range
        if snp_start == None: return
        
        contig = snp.chromosome
        
        transcripts = list(self.mTranscripts.get( snp.chromosome, snp_start, snp_end))
        
        if not transcripts: return

        reference_base = snp.reference_base

        ## collect all variants at this position
        ## indels and deletions might effect more than this
        ## position
        variants_to_test = []

        variant_types = []
        
        is_homozygous = True
        
        if reference_base == "*":
            variants = snp.genotype.split("/")
            codes = [ x[0] for x in variants ]

            # variant is hetorozygous if wildtype is present, codes/sequences of
            # variants are not identical.
            if ("*" in codes) or (variants[0] != variants[1]):
                is_homozygous = False

            # note that I found an inconsistency between the genotype field and the second-allele field
            # genotype='-GGG/-GGG', first_allelle='-GGG', second_allele='-GGGG'
            # In other cases it is correct, even with longer deletions.
            for variant in set(variants):
                if variant[0] == "*":
                    variant_types.append( "W")
                elif variant[0] == "+":
                    variant_types.append( "I")
                    # insertions affect the base before and after the insertion
                    variants_to_test.append( Variant._make( (variant[0], variant[1:], snp.pos, snp.pos + 1)) )
                elif variant[0] == "-":
                    variant_types.append( "D") 
                    # deletions are after the base denoted by snp.pos
                    start = snp.pos + 1 
                    # pos + 1 + len(var) - 1 = pos + len(var)
                    end = snp.pos + len(variant)
                    variants_to_test.append( Variant._make( (variant[0], variant[1:], start, end) ) )
        else:
            if snp.genotype in 'ACGTacgt':
                # homozygous substitution
                variant_types.append( "O" )
            else:
                # heterozygous substitution
                variant_types.append( "E" )
                is_homozygous = False

            for base in Genomics.resolveAmbiguousNA( snp.genotype ).upper():
                if base == snp.reference_base: continue
                variants_to_test.append( Variant._make( ("=", base, snp.pos, snp.pos+1 ) ) )

        self.mVariantTypes = variant_types
        E.debug( "snp: %s:%i variants_to_test=%i, transcripts=%i, is_homozygous=%s" %\
                 (snp.chromosome, snp.pos,
                  len(variants_to_test), len(transcripts), str(is_homozygous) ))

        counts = E.Counter()

        ## intersect all transcripts in the gene with the possible substitutions
        for transcript_start, transcript_end, exons in transcripts:

            transcript_id = exons[0].transcript_id

            all_splice_changes, all_splice_effects, all_cds_effects = [], [], []
            for variant in variants_to_test:

                E.debug( "snp: %s:%i variant=%i:%i:%s:%s, transcript=%s" % (snp.chromosome, snp.pos,
                                                                            variant.start,
                                                                            variant.end,
                                                                            variant.code,
                                                                            variant.sequence,
                                                                            transcript_id))


                r = getCDSPosition( exons,
                                    variant.start, variant.end,
                                    self.mFasta )

                if not r: continue

                reference_base, variant_seq, variant_bases = self.getSequence( snp, r, variant )
                
                # assert variant_seq.lower() in r.cds_seq.lower(), \
                #     "variant sequence %s not in cds seq %s: %s" % (variant_seq, r.cds_seq, str(r))

                cds_effects = self.collectCodingEffects( snp, r, variant,
                                                         reference_base, variant_seq, variant_bases )

                splice_effects, splice_changes = self.collectSplicingEffects( snp, r, variant,
                                                                              reference_base, variant_seq, variant_bases )
                
                if len(splice_effects) + len(cds_effects) == 0:
                    counts.no_effect += 1
                    continue

                all_splice_effects.extend( splice_effects )
                all_cds_effects.extend( cds_effects )
                all_splice_changes.extend( splice_changes )
                
            if all_splice_changes:
                self.outputSpliceEffects( snp, exons, variant, all_splice_changes, r )
            if all_cds_effects:
                self.outputCDSEffects( snp, exons, variant, all_cds_effects, r )

            if len(all_splice_effects) + len(all_cds_effects) == 0:
                continue
            
            self.updateVariantTranscripts( transcript_id, snp,
                                           exons, variant,
                                           all_splice_effects,
                                           all_cds_effects,
                                           r, is_homozygous )

    def updateVariantTranscripts( self, transcript_id, snp, exons, variant, splice_effects, cds_effects, r, is_homozygous):
        '''collect variation for each transcript.
        '''

        if transcript_id not in self.mTranscriptVariants:
            self.mTranscriptVariants[transcript_id] = TranscriptVariant._make( ([], []) )
        v = self.mTranscriptVariants[transcript_id] 
            
        for e in cds_effects:
            # splice variants cause all residues after a modified splice site to be deleted             
            v.cds_variants.append( 
                CdsVariant._make( (transcript_id,
                                   r.cds_start, r.cds_end,
                                   variant.code,
                                   e.orig_seq, e.variant_seq,
                                   is_homozygous) ))
            
        for e in splice_effects:
            # for splice effects save the full snps to sort out the intron sequence later.
            # due to deletions, etc, the resolving might be difficult.
            v.splice_variants.append( 
                SpliceVariant._make( (transcript_id,
                                      e.exon_id,
                                      r.nc_start - r.intron_start, r.nc_end - r.intron_start,
                                      variant.code,
                                      e.orig_seq, e.variant_seq,
                                      is_homozygous ) ) )
            
                                                

        
            
    def getSpliceCode( self, splice_name, new_splice_name ):
        '''assign one-letter code to a splice-signal change.'''
        
        if splice_name == "unknown" and new_splice_name == "unknown":
            # unknown splice signal
            code = "U"
        elif new_splice_name == "unknown":
            # disrupted splice site
            code = "D"
        elif splice_name == "unknown":
            # newly created splice site
            code = "C"
        elif splice_name == new_splice_name:
            # synonymous change
            code = "S"
        elif splice_name != new_splice_name:
            # non-synonymous change
            code = "N"
        return code
    
    def outputSpliceEffects( self, snp, exons, variant, splice_effects, r ):
        '''output effects of variants affecting splice sites.'''
        
        for e in splice_effects:
            self.mOutfileIntron.write( "%s\n" % "\t".join( \
                ( exons[0].transcript_id,
                  snp.chromosome,
                  "%i" % snp.pos,
                  ",".join( self.mVariantTypes ),
                  variant.code,
                  variant.sequence,
                  "%i" % e.exon_id,
                  "%i" % len(exons),
                  self.getSpliceCode( e.orig_name, e.variant_name ),
                  str(e.orig_name),
                  e.orig_seq5,
                  e.orig_seq3,
                  str(e.variant_name),
                  e.variant_seq5,
                  e.variant_seq3,
                  "%i" % r.intron_start,
                  "%i" % r.intron_end,
                  r.strand,
                  "%i" % r.nc_start,
                  "%i" % r.nc_end,
                  ) ) )

    def getSubstitutionCode( self, original_codons, variant_codons ):
        '''assign one-letter code codon change.

        '''

        if variant_codons == "!":
            # variant creates a frameshift
            code = "F"
        elif original_codons == variant_codons:
            # a synonymous substitution
            code = "S"
        elif "X" in variant_codons:
            # variant creates a stop codon
            code = "X"
        elif "U" in variant_codons:
            # variant creates a stop codon - that might a selenocysteine
            code = "U"
        elif original_codons in variant_codons:
            # a synonymous insertion
            code = "I"
        elif len(variant_codons) == 0 or variant_codons in original_codons:
            # a synonymous deletion
            code = "D"
        else:
            # a non-synonymous variant (substition or indel)
            # removing the original codon and replacing it with others
            code = "N"
        return code
    
    def outputCDSEffects( self, snp, exons, variant, cds_effects, r ):

        cds_len = sum( [x.end - x.start for x in exons ] )

        is_seleno = exons[0].transcript_id in self.mSeleno         
        for e in cds_effects:

            assert len(e.codon_orig_seq) % 3 == 0
            assert e.codon_orig_seq != e.codon_variant_seq
            
            orig_codons = Genomics.translate( e.codon_orig_seq,
                                              is_seleno = is_seleno )
            
            if len(e.codon_variant_seq) % 3 == 0:
                variant_codons = Genomics.translate( e.codon_variant_seq,
                                                     is_seleno = is_seleno )
            else:
                variant_codons = "!"
            
            self.mOutfileCds.write( "%s\n" % "\t".join( \
                ( exons[0].transcript_id,
                  snp.chromosome,
                  "%i" % snp.pos,
                  snp.reference_base,
                  ",".join( self.mVariantTypes ),
                  variant.code,
                  variant.sequence,
                  "%i" % e.exon_id,
                  "%i" % len(exons),
                  self.getSubstitutionCode( orig_codons, variant_codons ),
                  str(e.orig_seq),
                  str(e.codon_orig_seq),
                  orig_codons,
                  str(e.variant_seq),
                  str(e.codon_variant_seq),
                  variant_codons,
                  "%i" % r.cds_phase,
                  "%i" % r.cds_start,
                  "%i" % r.cds_end,
                  "%i" % cds_len)))

    def buildCDSVariantsPerPosition( self, transcript_id, cds_variants, cds_len ):
        '''count the number of variants.
        
        '''
        
        variants_per_position = numpy.zeros( cds_len )

        ncds_variants = len(cds_variants)

        for v in cds_variants:
            assert v.cds_end <= cds_len
            variants_per_position[v.cds_start:v.cds_end] += 1

        return variants_per_position

    def buildIntronsVariantsPerPosition( self, transcript_id, variants, intron_seqs ):
        '''count the number of within introns

        (variants have already been filtered to only include those that
        affect splicing).
        '''
        s = self.mSize
        lengths = [ len(x) for x in intron_seqs ]
        # only count 2 * s positions within intron
        variants_per_position = numpy.zeros( 2 * s * len(lengths) )

        nvar = len(variants)
        
        for v in variants:
            offset = v.intron_id * 2 * s
            l = lengths[v.intron_id]
            start, end = v.nc_start, v.nc_end
            if start < s:
                assert end < s, "variant (%i) larger than mSize (%i)" % (end, s)
            elif l - end < s:
                assert l - start < s, "variant (%i) larger than mSize (%i)" % (l-start, s)
                offset += s
                start, end = l-end, l-start
            else:
                raise ValueError("count out of range")
            
            variants_per_position[offset + start:offset + end] += 1
            
        return variants_per_position

    def getGenotype( self, variants, variants_per_position, counts ):
        '''compute the genotype and number of variants.

        *variants_per_position* is a vector of variants affecting a position.

        returns a genotype and the number of variants.
        '''

        max_variants_per_position = max(variants_per_position)

        nvar = len(variants)
        homo = [ x.is_homozygous for x in variants ]
        nhomo = len([x for x in homo if x])
        nhetero = len(homo) - nhomo 
        if nhomo == nvar and max_variants_per_position == 1:
            # all are homozygous, one variant only
            genotype = "O"
            counts.is_homozygous += 1
            counts.is_resolvable += 1
            nvariants = 1
        elif nhomo == 0 and nvar == 1 and max_variants_per_position == 1:
            # one heterozygous position, rest is wild type
            genotype = "W"
            counts.is_heterozygous += 1
            counts.is_resolvable += 1
            nvariants = 1
        elif nhomo == nvar-1 and max_variants_per_position == 1:
            # one heterozygous allowed if the rest are homozygous
            genotype = "E"
            counts.is_heterozygous += 1
            counts.is_resolvable += 1
            nvariants = 2
        elif nvar == 1 and max_variants_per_position == 2:
            # if there is only one heterozygous variant, which does not include the wild type
            genotype = "E"
            counts.is_heterozygous += 1
            counts.is_resolvable += 1
            nvariants = 2
        elif nhetero == nvar and max_variants_per_position == 1:
            # if all are heterozygous and one allele is always the wild type
            # resolve towards one allele, though it might be ambiguous
            genotype = "V"
            counts.is_heterozygous += 1
            counts.is_ambiguous += 1
            nvariants = 1
        elif max_variants_per_position == 1:
            # if there is only one variant at each position but more than two 
            # heterozygous variants in total
            # resolve towards two alleles
            genotype = "v"
            counts.is_heterozygous += 1
            counts.is_ambiguous += 1
            nvariants = 2
        else:
            genotype = "M"
            counts.is_mixture += 1
            counts.is_unresolvable += 1
            nvariants = 2
            
        return genotype, nvariants, max_variants_per_position

    def buildCDSVariants( self, 
                          transcript_id,
                          cds_variants,
                          reference_seq_na,
                          offset,
                          nvariants ):
        '''build variants for the coding sequence.

        offset: offset to correct for starting frame != 0
        '''

        variant_cds_seqs = []

        # the following code works with two variants at most
        assert 0 < nvariants <= 2, "expected 1 or 2 variants, got %i" % nvariants

        for x in range(nvariants):
            variant_cds_seqs.append( list(reference_seq_na) )
        n = 0
        for v in cds_variants:

            # ignore variants at incomplete codons
            if v.cds_start+offset < 0:
                E.warn("skipping variant in %s in first out-frame codon: %s." % (transcript_id, str(v)))
                continue
            
            if v.is_homozygous: toupdate = range(nvariants)
            else: toupdate = (0,)

            if v.code == "=":
                assert len(v.variant_seq) == 1
                assert reference_seq_na[v.cds_start+offset] == v.reference_seq.lower(), "transcript %s: base mismatch: %s != %s at %i, %s" %\
                       (transcript_id, reference_seq_na[v.cds_start+offset], v.reference_seq.lower(), v.cds_start, str(v))
                for x in toupdate:
                    variant_cds_seqs[x][v.cds_start+offset] = v.variant_seq
            elif v.code == "+":
                # indels are done without disrupting the frame
                # prepend.
                for x in toupdate:
                    variant_cds_seqs[x][v.cds_start+offset] = v.variant_seq + variant_cds_seqs[x][v.cds_start+offset]
            elif v.code == "-":
                # indels are done without disrupting the frame
                for x in toupdate:
                    for y in range(v.cds_start, v.cds_end):
                        variant_cds_seqs[x][y+offset] = ""
            n += 1

        if E.global_options.loglevel >= 10:
            for x in range(nvariants):
                Genomics.printPrettyAlignment( reference_seq_na, 
                                               variant_cds_seqs[x] )

        return variant_cds_seqs

    def buildIntronVariants( self, transcript_id, splice_variants,
                             reference_seqs_na, nvariants ):
        '''build all intron variants.

        Returns a list of variants. Each variant is a list of introns. Introns that are unchanged
        are None.

        The first entry in the list is the wildtype.

        returns a list of variants.
        '''
        variant_intron_seqs = []

        # the following code works with one or two variants
        assert 0 < nvariants <= 2, "expected 1 or 2 variants, got %i" % nvariants
        nintrons = len(reference_seqs_na)

        for x in range(nvariants):
            variant_intron_seqs.append( [ None for y in reference_seqs_na] )

        n = 0
        for v in splice_variants:

            E.debug( "transcript_id=%s: splice=%s" % (transcript_id, str(v) ))

            if v.is_homozygous: toupdate = range(nvariants)
            else: toupdate = (0,)

            intron_id = v.intron_id
            assert 0 <= intron_id < len(reference_seqs_na), "intron id `%i` out of range" % intron_id

            # instantiate intron sequence
            for x in toupdate:
                if variant_intron_seqs[x][intron_id] == None:
                    variant_intron_seqs[x][intron_id] = list(reference_seqs_na[intron_id])
                    
            if v.code == "=":
                assert len(v.variant_seq) == 1
                assert reference_seqs_na[intron_id][v.nc_start] == v.reference_seq.lower(), \
                       "transcript %s: base mismatch: %s != %s at %i:%i" %\
                       (transcript_id, reference_seqs_na[v.nc_start], v.reference_seq.lower(),
                        v.intron_id, v.nc_start)
                for x in toupdate:
                    variant_intron_seqs[x][intron_id][v.nc_start] = v.variant_seq
            elif v.code == "+":
                # indels are done without disrupting the frame
                # prepend to second residue
                assert (v.nc_end - v.nc_start) == 2
                for x in toupdate:
                    variant_intron_seqs[x][intron_id][v.nc_end] = v.variant_seq + \
                        variant_intron_seqs[x][intron_id][v.nc_end]
            elif v.code == "-":
                # indels are done without disrupting the frame
                for x in toupdate:
                    for y in range(v.nc_start, v.nc_end):
                        variant_intron_seqs[x][intron_id][y] = ""
            n += 1

        return variant_intron_seqs

    def countEffectsOnSplicing( self, variant_intron_seqs, reference_intron_seqs, min_intron_size = 5 ):
        '''collect all effects per intron

        return a count for each intron.
        '''

        nintrons = len(reference_intron_seqs)
        ncorrected_frames, nsynonymous, nnonsynonymous, ncanonical = 0, 0, 0, 0
        ndisrupted, nunknown, nunchanged, nnovel = 0, 0, 0, 0
        ncorrected_frames, nuncorrected_frames = 0, 0
        nframeshifts, ninserted_codons, nunchanged_frames = 0, 0, 0
        nnoncanonical = 0
        codes = []
        last_exon = nintrons + 1

        for intron_id, reference_seq in enumerate(reference_intron_seqs):
            
            reference_name, reference_seq5, reference_seq3 = Genomics.GetIntronType( reference_seq )
            e = 0

            variant_seq = variant_intron_seqs[intron_id]
            # process frameshift introns
            if len(reference_seq) < min_intron_size: 

                nframeshifts += 1

                if variant_seq == None:
                    variant_name, variant_seq5, variant_seq3 = reference_name, reference_seq5, reference_seq3
                    nunchanged_frames += 1
                    codes.append( "." )
                    continue

                variant_seq = "".join( variant_seq )

                # there might be both sequences of mod 3 and not
                fullseq = "".join( variant_seq )
                
                if len(variant_seq) % 3 == 0:
                    # a fixed frame shift
                    ncorrected_frames += 1
                    # note that the inserted codon sequence might contian stops
                    # needs to be tested with the other exons as it might not be
                    # in frame.
                    code = "F"
                    ninserted_codons += len(variant_seq) // 3
                else:
                    code = "P"
                    nuncorrected_frames += 1

            # process real introns
            else:
                if reference_name != "unknown":
                    ncanonical += 1
                else:
                    nnoncanonical += 1
                    
                if variant_seq == None:
                    variant_name, variant_seq5, variant_seq3 = reference_name, reference_seq5, reference_seq3
                    nunchanged += 1
                    codes.append( "." )
                    continue

                variant_seq = "".join( variant_seq )

                variant_name, variant_seq5, variant_seq3 = Genomics.GetIntronType( variant_seq )

                code = self.getSpliceCode( reference_name, variant_name )
                if code == "D":
                    last_exon = min( last_exon, intron_id )
                    ndisrupted += 1
                elif code == "C": nnovel += 1
                elif code == "N": nnonsynonymous += 1
                elif code == "S": nsynonymous += 1
                elif code == "U": nunknown += 1

            codes.append( code )

        return SplicingEffect._make( (nintrons,
                                      ncanonical,
                                      nframeshifts,
                                      nnoncanonical,
                                      nunchanged_frames,
                                      ncorrected_frames,
                                      nuncorrected_frames,
                                      nunchanged,
                                      nsynonymous,
                                      nnonsynonymous,
                                      ndisrupted,
                                      nnovel,
                                      nunknown,
                                      ninserted_codons,
                                      "".join( codes ),
                                      last_exon) )

    def getTruncatedCodons( self, is_homozygous, stops, ncodons ):
        '''return codons that are truncated due to stop codons.

        Note that if two variants are present and there is
        a homozygous variant causing a stop codon, both variants
        will have the same stop codon registered automatically.

        return for each variant.
        '''
        if len(stops) == 0: return 0, 0
        
        # one stop - one variant
        if len(stops) == 1:
            # if homozygous: both allelles have the stop
            if is_homozygous:
                stop_min = stop_max = ncodons - stops[0]
            else: # wildtype still present
                stop_min, stop_max = 0, ncodons - stops[0]
        else:
            stop_min, stop_max = ncodons - max(stops), ncodons - min(stops)

        return max(0,stop_min), max(0,stop_max)

    def fixCDSTermini( self, variant_cds_seqs, contig, strand, start, end ):
        '''if the first codon in a sequence has been deleted, add
           sequence from the UTR.

           Not implemented yet - needs to take into account indels in 
           the UTR as well.
           '''

        return variant_cds_seqs
    
        for vairant_cds_seq in variant_cds_seqs:
            x = 0

            # find first base that is not deleted
            while x < len(variant_cds_seq) and variant_cds_seq[x] == "": x += 1

            # note: to be correct, this should take into account indels as well.
            extra_seq = self.mFasta.getSequence( contig, strand, start - x, start )

            for xx in range(0,xx):
                pass
    

    def writeTable( self, outfile ):
        '''output summary for each transcript.

        Output three tables;

        1. mOutfileTranscripts: translation information
        2. mOutfileSplicing: splicing 
        3. mOutfile: counts
        '''

        cds_counts = E.Counter()
        splice_counts = E.Counter()
        all_counts = E.Counter()

        # TODO: the current code is not consistent when it comes
        # to counting premature stop codons as it also includes
        # the wild type as variant 0.

        for transcript_id, exons in self.mExons.iteritems():

            ###################################################
            ###################################################
            ###################################################            
            ## sort out exons and get some chromosomal coordinates
            exons = self.mExons[transcript_id]
            exons.sort( key = lambda x: x.start )
            cds_len = sum( [x.end - x.start for x in exons ] )
            ncodons = cds_len // 3
            contig = exons[0].contig
            lcontig = self.mFasta.getLength( contig )
            strand = exons[0].strand
            is_positive_strand = Genomics.IsPositiveStrand( strand )
                
            # obtain cds sequences            
            reference_seq_na = GTF.toSequence( exons, self.mFasta ).lower()

            # obtain intron sequences
            intron_intervals = GTF.toIntronIntervals( exons )
                
            if not is_positive_strand: 
                intron_intervals = [ (lcontig - end, lcontig - start) for start,end in intron_intervals ]
                intron_intervals.reverse()

            intron_sequences = [ self.mFasta.getSequence( contig, strand, x[0], x[1]).lower() for x in intron_intervals ]
            nintrons = len(intron_intervals)
            
            is_seleno = transcript_id in self.mSeleno         

            # result variables - set to wildtype
            all_genotype, all_nalleles, all_max_variants = "", 0, 0
            cds_genotype, cds_nalleles, cds_max_variants = "", 0, 0
            splice_genotype, splice_nalleles, splice_max_variants = "", 0, 0
            cds_nvariant_positions, splice_nvariant_positions = 0, 0
            variant_intron_seqs, splice_variants_per_position = [], []
            variant_cds_seqs, cds_variants_per_position = [], []

            exon2cds = []
            if is_positive_strand:
                frame = int(exons[0].frame)
                cds_pos = frame
                for x in exons:
                    exon2cds.append( cds_pos )
                    cds_pos += x.end - x.start
            else:
                frame = int(exons[-1].frame)
                cds_pos = frame
                for x in exons[::-1]:
                    exon2cds.append( cds_pos )
                    cds_pos += x.end - x.start
            last_exon_start = exon2cds[-1]
            exon2cds.append( cds_len )

            if transcript_id in self.mTranscriptVariants:
                
                variants = self.mTranscriptVariants[transcript_id]

                E.debug( "processing %s with %i cds effects and %i splice effects started" % \
                         (transcript_id, len(variants.cds_variants), len(variants.splice_variants) ) )

                # we should have some variants
                assert len(variants.cds_variants) + len(variants.splice_variants) > 0

                # correct for frame at start - truncate the reference_seq_na
                if frame != 0:
                    E.debug( "transcript_id %s - correcting frame %i" % (transcript_id, frame) )
                    reference_seq_na = reference_seq_na[frame:]
                    # all coordinates need to modified by this amount
                    offset = -frame
                else:
                    offset = 0

                reference_seq_aa = Genomics.translate( reference_seq_na, 
                                                       is_seleno = is_seleno )

                cds_nvariant_positions = len(variants.cds_variants)
                splice_nvariant_positions = len(variants.splice_variants)
                
                ###################################################
                ###################################################
                ###################################################
                ## build coding sequence variants
                if len(variants.cds_variants) > 0:

                    ###################################################
                    ###################################################
                    ###################################################            
                    # decide what variants to build
                    # 1. homozygous: 1 variant per position, all also flagged as homozygous
                    # 2. heterozygous + wildtype: only 1 variant per position, all flagged as heterozygous
                    # 3. heterozygous: 2 variants per position, but only if there is only one position modified
                    # 4: mixture: rest
                    cds_variants_per_position = self.buildCDSVariantsPerPosition(
                        transcript_id,
                        variants.cds_variants,
                        cds_len )

                    cds_genotype, cds_nalleles, cds_max_variants = self.getGenotype(
                        variants.cds_variants,
                        cds_variants_per_position,
                        cds_counts ) 

                    variant_cds_seqs = self.buildCDSVariants( transcript_id,
                                                              variants.cds_variants,
                                                              reference_seq_na,
                                                              offset,
                                                              cds_nalleles )

                ###################################################
                ###################################################
                ###################################################            
                ## build intron variants

                ###################################################
                ###################################################
                ###################################################            
                ## collect all intron sequences
                if len(variants.splice_variants) > 0:
                    ###################################################
                    ###################################################
                    ###################################################            
                    ## collect genotype and variants to build
                    splice_variants_per_position = self.buildIntronsVariantsPerPosition(
                        transcript_id,
                        variants.splice_variants,
                        intron_sequences )

                    splice_genotype, splice_nalleles, splice_max_variants = self.getGenotype( 
                        variants.splice_variants,
                        splice_variants_per_position,
                        splice_counts )

                    variant_intron_seqs = self.buildIntronVariants(
                        transcript_id,
                        variants.splice_variants,
                        intron_sequences,
                        splice_nalleles )

                ###################################################
                ###################################################
                ###################################################            
                ## collect overall genotype
                all_genotype, all_nalleles, all_max_variants = self.getGenotype( 
                    variants.cds_variants + variants.splice_variants,
                    numpy.concatenate( (cds_variants_per_position, splice_variants_per_position) ),
                    all_counts )

            ###################################################
            ###################################################
            ###################################################            
            # add the wild type at top of both cds and intron variants
            #
            # This is necessary so that stop codons originally present
            # in the sequence will be taken into account.
            #
            # Note that this invalidates the cds_stop_min below.
            #
            # A better way would be to merge variants and only 
            # add the wild type if there is only variant allele.
            #
            # Then, treat the wildtype separately to get numbers for
            # for the wildtype.
            if len(variant_cds_seqs) == 0:
                variant_cds_seqs = [ list(reference_seq_na),
                                     list(reference_seq_na) ]

            elif len(variant_cds_seqs) == 1:
                if cds_genotype == "O":
                    # is homozygous - duplicate allele
                    variant_cds_seqs.append( variant_cds_seqs[0] )
                else:
                    # add wildtype
                    variant_cds_seqs[0:0] = [ list(reference_seq_na),]

            if len(variant_intron_seqs) == 0:
                variant_intron_seqs = [ [ None for x in range(nintrons) ],
                                        [ None for x in range(nintrons) ] ]

            elif len(variant_intron_seqs) == 1:
                if splice_genotype == "O":
                    # is homozygous - duplicate allele
                    variant_intron_seqs.append( variant_intron_seqs[0] )
                else:
                    # add wildtype
                    variant_intron_seqs[0:0] = [ [ None for x in range(nintrons) ],]

            assert len(variant_cds_seqs) == 2
            assert len(variant_intron_seqs) == 2

            ###################################################
            ###################################################
            ###################################################            
            ## output information on splice/cds variants per transcript
            ## output also the wild type (variant_id = 0)
            ###################################################                            
            cds_stops, splice_stops = [], []

            for variant_id, variant_seq in enumerate(variant_intron_seqs):

                variant_result = self.countEffectsOnSplicing( variant_seq,
                                                              intron_sequences )

                self.mOutfileSplicing.write( "%s\t%i\t%s\n" % \
                                  (transcript_id,
                                   variant_id,
                                   "\t".join( map(str, variant_result) ) ) )

                splice_stops.append( exon2cds[variant_result.last_exon] // 3 )

            # estimate effect on protein coding sequence for each variant and output
            for variant_id, variant_seq in enumerate(variant_cds_seqs):
                variant_result = countEffectsOnTranscript( variant_seq,
                                                           reference_seq_na,
                                                           is_seleno = is_seleno )
                s = "".join( variant_seq )
                self.mOutfileTranscripts.write( \
                    "%s\t%i\t%i\t%s\t%s\t%s\n" %\
                    (transcript_id,
                     variant_id,
                     last_exon_start,
                     "\t".join( map(str, variant_result) ),
                     "".join( s ),
                     Genomics.translate( s, is_seleno = is_seleno ),
                     ) )

                cds_stops.append( variant_result.first_stop )

            ###################################################
            ###################################################
            ###################################################            
            ## compute the shortest transcript variants
            ## due to splicing and cds changes separately and
            ## combined.
            ###################################################                                
            if splice_nalleles > 0:
                splice_stop_min, splice_stop_max = \
                                 self.getTruncatedCodons( splice_genotype == "O", splice_stops, ncodons )
            else:
                splice_stop_min, splice_stop_max = 0, 0

            if cds_nalleles > 0:
                cds_stop_min, cds_stop_max = \
                    self.getTruncatedCodons( cds_genotype == "O", cds_stops, ncodons )
            else:
                cds_stop_min, cds_stop_max = 0, 0

            # combine stops between cds and slice variants
            # the two variants will have the overall maxima
            all_stop_min, all_stop_max = ( max(splice_stop_min, cds_stop_min),
                                           max(splice_stop_max, cds_stop_max) )

            ###################################################
            ###################################################
            ###################################################            
            # output stats per transcript
            ###################################################            
            outfile.write( "%s\n" % "\t".join( ( \
                transcript_id,
                "%i" % cds_len,
                "%i" % ncodons,
                "%i" % last_exon_start,
                "%i" % cds_max_variants,
                "%i" % cds_nvariant_positions,
                "%s" % cds_genotype,
                "%i" % cds_nalleles,
                "%i" % cds_stop_min,
                "%i" % cds_stop_max,
                "%i" % splice_max_variants,
                "%i" % splice_nvariant_positions,
                "%s" % splice_genotype,
                "%i" % splice_nalleles,
                "%i" % splice_stop_min,
                "%i" % splice_stop_max,
                "%i" % all_max_variants,
                "%i" % (cds_nvariant_positions + splice_nvariant_positions),
                "%s" % all_genotype,
                "%i" % all_nalleles,
                "%i" % all_stop_min,
                "%i" % all_stop_max,
                )))

            E.debug( "processing %s with %i cds effects and %i splice effects finished" % \
                     (transcript_id, cds_nvariant_positions, splice_nvariant_positions ))
                
        E.info( "cds counts: %s" % (str(cds_counts)) )
        E.info( "splice counts: %s" % (str(splice_counts)) )
        E.info( "combined counts: %s" % (str(all_counts)) )

class CounterContigs( Counter ):
    '''count variants across the genome per chromosome.'''

    mHeader = [ "genome_%s" % x for x in ( "ntranscripts", "nused", "pos" )]

    def __init__(self, *args, **kwargs):
        Counter.__init__(self, *args,**kwargs)

        # create counter
        self.mCountsSNPs = collections.defaultdict( int )
        self.mCountsIndels = collections.defaultdict( int )

    def update( self, snp ):
        '''update with snp.'''

        if snp.reference_base == "*":
            self.mCountsIndels[snp.chromosome] += 1
        else:
            self.mCountsSNPs[snp.chromosome] += 1

    def writeTable( self, outfile ):

        outfile.write( "contig\tsize\tnindels\tnsnps\n" )
        total_snps, total_indels, total_length = 0, 0, 0
        for key in sorted(self.mCountsSNPs.keys()):
            total_snps += self.mCountsSNPs[key]
            total_indels += self.mCountsIndels[key]
            total_length += self.mFasta.getLength(key)
            outfile.write( "\t".join( (key,
                                       "%i" % self.mFasta.getLength(key),
                                       "%i" % self.mCountsIndels[key],
                                       "%i" % self.mCountsSNPs[key] ) ) + "\n" )

        outfile.write( "\t".join( ("total",
                                   "%i" % total_length,
                                   "%i" % total_indels,
                                   "%i" % total_snps ) ) + "\n" )



def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: snp2counts.py 2872 2010-03-03 10:21:13Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default]."  )
    parser.add_option("-f", "--filename-exons", dest="filename_exons", type="string",
                      help="filename with exon information (gtf formatted file)  [default=%default]."  )
    parser.add_option("-s", "--filename-seleno", dest="filename_seleno", type="string",
                      help="filename of a list of transcript ids that are selenoproteins [default=%default]."  )
    parser.add_option("-c", "--filename-vcf", dest="filename_vcf", type="string",
                      help="vcf file to parse [default=%default]."  )
    parser.add_option("-m", "--module", dest="modules", type="choice", action="append",
                      choices=("gene-counts", "transcript-effects", "contig-counts"),
                      help="modules to apply [default=%default]."  )
    parser.add_option("-i", "--input-format", dest="input_format", type="choice",
                      choices = ("pileup", "vcf" ),
                      help="input format [default=%default]."  )
    parser.add_option( "--vcf-sample", dest="vcf_sample", type="string",
                      help="sample id in vcf file to analyse [default=%default]."  )

    parser.set_defaults(
        genome_file = None,
        filename_exons = None,
        filename_seleno = None,
        filename_vcf = None,
        modules = [],
        input_format = "pileup",
        vcf_sample = None,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

    ninput, nskipped, noutput = 0, 0, 0

    ################################
    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )
    else:
        fasta = None

    if options.filename_seleno:
        seleno = set( IOTools.readList( IOTools.openFile(options.filename_seleno, "r")) )
    else:
        seleno = {}

    # setup iterator
    if options.input_format == "pileup":
        iterator = pysam.Pileup.iterate( options.stdin)
    elif options.input_format == "vcf":
        if not options.vcf_sample:
            raise ValueError( "vcf format requires sample id (--vcf-sample) to be set" )
        if not options.filename_vcf:
            raise ValueError( "reading from vcf requires vcf filename (--filename-vcf) to be set)" )

        iterator = pysam.Pileup.iterate_from_vcf(options.filename_vcf, options.vcf_sample )

    ################################
    modules = []
    for module in options.modules:
        if module == "gene-counts":
            if not options.filename_exons:
                raise ValueError("please supply exon information (--filename-exons)")
            modules.append( CounterGenes( options.filename_exons, fasta=fasta ) )

        elif module == "transcript-effects":
            if not options.filename_exons:
                raise ValueError("please supply exon information (--filename-exons)")
            modules.append( CounterTranscripts( options.filename_exons, fasta=fasta,
                                                pattern = options.output_filename_pattern,
                                                seleno = seleno ) )
            
        elif module == "contig-counts":
            modules.append( CounterContigs( fasta=fasta ) )
            
    options.stdout.write( "\t".join( [x.getHeader() for x in modules]) + "\n" )

    for snp in iterator:
        ninput += 1

        # translate chromosome according to fasta
        if fasta: 
            snp = snp._replace( chromosome=fasta.getToken( snp.chromosome ))

        for module in modules:
            module.update( snp )

        # if ninput > 1000: break
        
    for module in modules:
        module.writeTable( options.stdout )
        
    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput,nskipped) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
