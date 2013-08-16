################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: gtf2alleles.py 2886 2010-04-07 08:47:46Z andreas $
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
gtf2alleles.py - predict effects of variants on gene set
========================================================

:Author: Andreas Heger
:Release: $Id: gtf2alleles.py 2886 2010-04-07 08:47:46Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

build from a gene set in gtf format an allelic gene set in which
each transcript is split into separate alleles. 

Alleles are built from variants stored in an sqlite database (see the
options ``--database`` and ``--tablename``).

Alleles are built from transcripts in a gene set (*--filename-exons*).

Caveats
-------

* no phasing - heterozygous variants are assigned to the first or second
  allele simply based on the input code.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

""" 

import os
import sys
import re
import optparse
import collections
import sqlite3

import numpy
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Genomics as Genomics
import CGAT.GTF as GTF
import CGAT.Blat as Blat
import CGAT.Variants as Variants
import alignlib
import pysam

Allele = collections.namedtuple( 'Allele',
                                 '''cds, 
                                    peptide, 
                                    nexons,
                                    cds_starts,
                                    exon_starts,
                                    frames,
                                    is_nmd_knockout,
                                    is_splice_truncated,
                                    is_stop_truncated,
                                    nframeshifts,
                                    ncorrected_frameshifts,
                                    nuncorrected_frameshits,
                                    peptide_first_stop,
                                    peptide_len,
                                    cds_first_stop,
                                    cds_len,
                                    reference_first_stop_start,
                                    reference_first_stop_end,
                                    cds_original_len,
                                    nsplice_noncanonical,
                                 ''')


class VariantGetter(object):
    '''base class for objects returning variants.'''
    pass

class VariantGetterSqlite(VariantGetter):
    '''retrieve variants from an sqlite table in pileup format.'''

    def __init__(self, dbname, tablename ):
        self.dbname = dbname
        self.tablename = tablename
        
        self.dbhandle = sqlite3.connect( dbname )

        self.statement = '''SELECT 
                   pos, reference, genotype 
                   FROM %(tablename)s
                   WHERE contig = '%(contig)s' AND 
                   pos BETWEEN %(start)s and %(end)s
                '''
        
    def __call__( self, contig, start, end ):

        cc = self.dbhandle.cursor()
        tablename = self.tablename
        cc.execute( self.statement % locals() )
        variants = map( Variants.Variant._make, cc.fetchall())
        cc.close()
        return variants

class VariantGetterPileup(VariantGetter):
    '''retrieve variants from file in pileup format.'''

    def __init__(self, filename ):
        self.tabix = pysam.Tabixfile( filename )
        self.contigs = set(self.tabix.contigs)

    def __call__(self, contig, start, end ):

        variants = []
        if contig not in self.contigs: return []

        for line in self.tabix.fetch( contig, start, end ):
            data = line[:-1].split()
            contig, pos, reference, genotype = data[:4]
            # fix 1-ness
            pos = int(pos) - 1
            variants.append( Variants.Variant._make( (pos, reference, genotype) ) )
        return variants
    
class VariantGetterVCF( VariantGetter ):
    '''retrieve variants from tabix indexed vcf file.'''
    
    def __init__(self, filename, sample ):
        self.sample = sample
        self.vcf = pysam.VCF()
        self.vcf.connect( filename )

        if sample not in self.vcf.getsamples():
            raise KeyError( "sample %s not vcf file" % sample )
        
    def __call__(self, contig, start, end ):
        
        variants = []
        s = self.sample

        try:
            iter = self.vcf.fetch( contig, start, end )
        except ValueError:
            # contigs not in variants, ignore
            return variants

        for row in iter:
            result = pysam.Pileup.vcf2pileup( row, s )
            if not result: continue
            variants.append( Variants.Variant._make( (result.pos, result.reference_base, result.genotype) ))

        return variants

def collectExonIntronSequences( transcripts, fasta ):
    '''collect all the wild type sequences for exons and introns
    
    exons and introns are indexed by their respective positions.

    The function changes coordinates in ``gtfs`` to reverse coordinates.
    '''

    contig = transcripts[0][0].contig 
    strand = transcripts[0][0].strand
    lcontig = fasta.getLength( contig )

    all_exons, all_introns = {}, {}
    for exons in transcripts:
        for exon in exons:
            exon.invert( lcontig )
            start, end = exon.start, exon.end
            key = start, end
            if key not in all_exons:
                all_exons[key] = fasta.getSequence( contig, strand, start, end ).lower()
                    
        intron_intervals = GTF.toIntronIntervals(exons)
        for start, end in intron_intervals:
            key = start, end
            if key not in all_introns:
                all_introns[key] = fasta.getSequence( contig, strand, start, end ).lower()
    return all_exons, all_introns


def buildCompactVariantSequences( variants, sequences ):
    '''build variant sequences by inserting ``variants`` into ``sequences``.

    The original frame of the sequence is maintained by
    converting the input sequence to a list. Each entry
    in the list corresponds to a position in a wild type.

    The wild type (WT) sequence is lower case

    SNP:                    variant (ambiguity codes for variants)
    homozygous insertion:   upper-case bases after lower-case (WT) base
    heterozygous insertion: lower-case bases after lower-case (WT) base
    homozygous deletion:    empty fields
    heterozygous deletion:  "-" after lower-case (WT) base

    returns a dictionary of lists.
    '''

    result = {}
    for key, sequence in sequences.iteritems():
        variant_seq = list( sequence.lower() )
        start, end = key
        
        # get all variants that overlap with sequences
        for var_start, var_end, values in variants.find( start, end ):
            reference, action, has_wildtype, variantseqs = values

            is_homozygous = len(variantseqs) == 1 and not has_wildtype
            
            rel_start, rel_end = var_start - start, var_end - start
            startoffset = max(0, start-var_start)
            endoffset = max(0,var_end-end)
            
            if action == "=":
                assert rel_start >= 0
                assert sequence[rel_start].upper() == reference, \
                    'reference base mismatch: expected %s, got %s at %i-%i' % \
                    (sequence[rel_start].upper(), reference, 
                     var_start, var_end)

                if is_homozygous:
                    variant_seq[rel_start] = variantseqs[0]
                else:
                    variant_seq[rel_start] = Genomics.resolveReverseAmbiguousNA( "".join(variantseqs) )
                    
            elif action == "-":

                xstart, xend = max(0,rel_start), min(len(sequence),rel_end)

                for variant in variantseqs:
                    # truncated for variants of unequal lengths (-AA/-AAA)
                    refseq = sequence[xstart:xend].upper()[:len(variant)]

                    assert refseq == variant[startoffset:len(variant) - endoffset], \
                    'reference base mismatch at deletion: expected %s %s %s, got %s[%i:%i] at %i-%i (%i-%i), action=%s' % \
                    (sequence[xstart-10:xstart],
                     refseq, 
                     sequence[xend:xend+10],
                     variant, startoffset, len(variant)-endoffset,
                     var_start, var_end, start, end, 
                     action)

                    l = len(variant) - startoffset - endoffset

                    if is_homozygous:
                        variant_seq[xstart:xend] = [""] * l
                    else:
                        for x in range(xstart,xend):
                            if variant_seq[x].endswith("-"):
                                assert not has_wildtype
                                variant_seq[x] = ""
                            else:
                                variant_seq[x] += "-"

            elif action == "+":

                if is_homozygous:
                    variant_seq[rel_start] += variantseqs[0].upper()
                else:
                    if has_wildtype:
                        variant_seq[rel_start] += variantseqs[0].upper()
                    else:
                        # merge indels like +AAA/+AA
                        a,b = variantseqs
                        if a.startswith(b):
                            variant_seq[rel_start] += b.upper() + a[len(b):].lower()
                        elif b.startswith(a):
                            variant_seq[rel_start] += a.upper() + b[len(a):].lower()
                        else:
                            raise ValueError("don't know how to encode variant: %s" % variantseqs)

        result[(start,end)] = variant_seq
    return result

def buildVariantSequences( indexed_variants, sequences ):
    '''build variant sequences by inserting ``variants`` into ``sequences``.

    For each sequence, two alleles are returned. Both alleles are initialized
    as wildtype sequences. In the absence of any phasing information, variants 
    are preferably added to the second allele, such that the wild-type status
    of the first allele is preserved as much as possible

    returns a dictionary of lists.
    '''

    result = {}
    for key, sequence in sequences.iteritems():

        feature_start, feature_end = key

        variants = [ (x,y,) + z for (x,y,z) in indexed_variants.find( feature_start, feature_end )]
        allele1, allele2 = Variants.buildAlleles( sequence, 
                                                  variants, 
                                                  reference_start = feature_start )

        result[(feature_start,feature_end)] = (allele1,allele2)

    return result

def buildCDSSequence( transcript, exons ):
    '''build transcript sequences from exons.'''

    cds = []
    for exon in transcript:
        cds.append( exons[exon.start,exon.end] )
    return "".join(cds)

def buildAlleles( transcript, 
                  variant_exons, 
                  variant_introns,
                  reference_exons,
                  reference_introns,
                  offsets,
                  is_seleno = False,
                  frameshiftsize = 5,
                  reference_coordinates = False):
    '''reconstitute variant transcript sequences from
    variant exons and introns. This method returns alleles
    as they are estimated to look like taking into account
    splice effects and nonsense mediated decay.

    Transcripts are build in the following way:
    
    1. The exonic sequences are examinated for stop codons.

    If they contain a stop-codon in the last exon, the sequence 
    will be truncated at the stop. If a stop-codon exists
    in exons before the last codon, the allele is thought to
    have been knocked out due to NMD and it is set to 0.

    2. If a splice-site is disrupted, the transcript is 
    truncated (it is not extended to the first stop, as
    splicing might be recovered).

    If ``variant_exons`` and ``variant_introns`` are set to
    None, the two wildtype alleles will be returned.

    Scenarios that this method ignores:
       1. insertion after 3' splice site are counted as belonging
       to the intron and not the exon, though they could be either.
       2. insertion before the 5' splice are counted as belonging
       to the exon and not the intron.

    If ``reference_coordinates`` is set to true, exon coordinates are
    taken from the reference. This ensures, that the base that is
    derived from the same reference sequence gets the same coordinate
    in the two alleles. The alternative is allelic coordinates, that
    take into account intron size changes. For example, the reference
    has two exons, (100,200) and (300,400). Due to variants, the
    first allele has now an intron size of 110. Thus, in allelic coordinates,
    the coordinates are ((100,200), (310,410)) and ((100,200),(300,400)),
    while in reference coordinates, both alleles whould have the same
    coordinates. Note that due to insertions/deletion, the coordinates
    might change within an exon, too.

    returns two alleles
    '''

    result = []

    def _buildAllele( allele_id, 
                      transcript, exons, 
                      introns, offsets, 
                      virtual_coordinates = False,
                      reference_exons = None ):

        def _getOffset( pos, offsets):
            x = 0
            while x < len(offsets) and offsets[x][0] <= pos: x += 1
            x -= 1
            if x >= 0: return offsets[x][1]
            else: return 0

        def _sumIndels( ss ):
            '''sum indels within ss'''
            c = 0
            for s in ss:
                c += len(s) - 1 
            return c

        def _getEndOffsets( ss ):
            '''get the offset at exons due to deletions at
            start/end of exon.'''
            l = len(ss)
            x = 0
            while x < l and ss[x] == "": x += 1
            start_offset = x

            x = l-1
            while x >= 0 and ss[x] == "": x -= 1
            if x >= 0: 
                return start_offset, (l-1) - x
            else:
                return start_offset, 0

        def _addCds2Reference( map_cds2reference,
                               cds_start,
                               cds_seq, 
                               reference_start ):
            '''add cds to reference'''
            c, r = cds_start, reference_start
            for x in cds_seq:
                l = len(x)
                if l == 0: 
                    r += 1
                else:
                    map_cds2reference.addPair( c, r )
                    c += l
                    r += 1
        # counts
        is_splice_truncated = False
        is_nmd_knockout = False
        is_stop_truncated = False
        nuncorrected_frameshifts = 0
        ncorrected_frameshifts = 0
        nframeshifts = 0
        nsplice_noncanonical = 0
        reference_first_stop_start = -1
        reference_first_stop_end = -1

        # map between the new cds sequence and the reference
        # sequence
        map_cds2reference = alignlib.makeAlignmentBlocks()

        ###################################################
        # process first exon
        exon = transcript[0]
        transcript_id = exon.transcript_id

        # collect offset for exon.start
        genome_start = exon.start
        genome_start += _getOffset( genome_start, offsets )
        lcds, cds = 0, []
        cds_starts = [0]

        # still need to deal with deletions of first base:
        exon_starts = [genome_start]
        exon_key = (exon.start, exon.end)
        exon_sequence = exons[exon_key] 
        exon_seq = "".join(exon_sequence) 
        
        cds.append(exon_seq)
        _addCds2Reference( map_cds2reference, 
                           lcds,
                           exon_sequence,
                           exon.start )
        lcds = len(exon_seq)

        if len(exon_seq) != exon.end - exon.start:
            nframeshifts += 1

        # add first exon to genome position
        genome_pos = genome_start + len(exon_seq)
        last_end = exon.end 

        # correct for deletions at start/end of exon
        start_offset, end_offset = _getEndOffsets( exon_sequence )

        # length of original transcript
        loriginal = sum( [x.end - x.start for x in transcript] )

        if E.global_options.loglevel >= 8:
            print "%i: exon_indels (%i-%i):" % (allele_id, exon.start, exon.end)
            for x,c in enumerate(exons[exon_key]):
                if len(c) != 1: print x + exon.start, ":%s:" % c 
            print
            print exons[exon_key]
            print "genome_pos=", genome_pos, \
                ",exon=%i-%i" % (genome_pos, genome_pos + len(exon_seq)), \
                ", len(exon_seq)=", len(exon_seq), \
                ", len(exon)=", exon.end - exon.start, \
                ", offsets=%i,%i," % (start_offset, end_offset), \
                ", offset at start=",_getOffset( exon.start, offsets), \
                ", offset at end=",_getOffset( exon.end, offsets)
            

        for exon in transcript[1:]:

            last_exon_sequence = exon_sequence
            last_start_offset, last_end_offset = start_offset, end_offset

            # get the next intron/exon parameters
            exon_key = (exon.start, exon.end)
            exon_sequence = exons[exon_key]
            start_offset, end_offset = _getEndOffsets( exon_sequence )
            intron_key = (last_end, exon.start)
            
            if last_end == exon.start:
                # catch empty introns
                intron_sequence = []
                intron_key = None
            else:
                intron_sequence = introns[intron_key] 

            intron_seq = "".join( intron_sequence )

            ###################################################
            ###################################################
            ###################################################
            # add preceding intron
            new_exon = True

            if len(intron_seq) > frameshiftsize:
                
                intron_name, intron_seq5, intron_seq3 = Genomics.GetIntronType( intron_seq )
                if intron_name == "unknown":
                    if intron_seq[:2].islower() and intron_seq[-2:].islower():
                        E.debug( "%s: transcript has unknown splice signal - kept because not a variant: %s: %s:%s" % \
                                     (transcript_id, intron_name, intron_seq5, intron_seq3 ) )
                        nsplice_noncanonical += 1
                    else:
                        is_splice_truncated = True
                        E.debug( "%s: transcript has splice truncated allele: %s: %s:%s" % \
                                     (transcript_id, intron_name, intron_seq5, intron_seq3 ) )
                        break
                # start a new exon
                cds_starts.append( lcds )
                
            else:
                # treat as frameshifting intron
                #
                # frame-shifting introns are checked if they are
                # fixed by indels either in the intron itself or
                # the terminal exon sequence. To this end, the effective
                # size of the intron is computed:
                # effective size of intron = 
                # indels at terminal x bases at previous exon
                # + size of intron
                # + indels at terminal x bases at next exon 
                effective_intron_size = len(intron_seq)
                previous_indels = _sumIndels( last_exon_sequence[max(0,-frameshiftsize):] )
                next_indels = _sumIndels( exon_sequence[:frameshiftsize] )
                effective_intron_size += previous_indels + next_indels

                if previous_indels + next_indels == 0 and len(intron_seq) % 3 == 0:
                    has_stop = "X" in Genomics.translate(intron_seq.upper(), 
                                                         is_seleno=is_seleno )
                else:
                    has_stop = False

                if effective_intron_size % 3 == 0 and not has_stop:
                    E.debug( "%s: fixed frame-shifting intron %i-%i of size %i (size:%i, indels:%i,%i)" % \
                                 (transcript_id, last_end, exon.start,
                                  effective_intron_size, 
                                  len(intron_seq),
                                  previous_indels, next_indels,))
                                  
                    # add to previous exon
                    cds.append( intron_seq )
                    lcds += len(intron_seq)
                    ncorrected_frameshifts += 1
                    new_exon = False
                else:
                    E.debug( "%s: could not fix frame-shifting intron %i-%i of size %i (size:%i, indels:%i,%i, has_stop=%i)" % \
                                 (transcript_id, last_end, exon.start,
                                  effective_intron_size, 
                                  len(intron_seq),
                                  previous_indels, next_indels,
                                  has_stop))

                    nuncorrected_frameshifts += 1
                    # start a new exon
                    cds_starts.append( lcds )

            if E.global_options.loglevel >= 8:
                print "%i: intron_indels (%i-%i):" % (allele_id, last_end, exon.start)
                if intron_key:
                    for x,c in enumerate(introns[intron_key]):
                        if len(c) != 1: print x + last_end, ":%s:" % c 
                    print
                    print introns[intron_key]
                    print "genome_pos=", genome_pos, \
                        ",intron=%i-%i" % (genome_pos, genome_pos + len(intron_seq)), \
                        ", len(intron_seq)=", len(intron_seq), \
                        ", len(intron)=", exon.start-last_end, \
                        ", offset at start=",_getOffset( last_end, offsets), \
                        ", offset at end=",_getOffset( exon.start, offsets)
                else:
                    print "empty intron"

            genome_pos += len(intron_seq)
            
            # assertion - check if genomic coordinate of intron is consistent with offset
            test_offset =_getOffset( exon.start, offsets )
            is_offset = genome_pos - exon.start
            assert is_offset == test_offset, "intron offset difference: %i != %i" % (is_offset, test_offset)

            ###################################################
            ###################################################
            ###################################################
            # add the exon
            exon_seq = "".join(exon_sequence)
            cds.append( exon_seq )

            if len(exon_seq) != exon.end - exon.start:
                nframeshifts += 1

            if new_exon:
                if reference_coordinates:
                    exon_starts.append( exon.start + start_offset )
                else:
                    exon_starts.append( genome_pos )

            _addCds2Reference( map_cds2reference, 
                               lcds, 
                               exon_sequence, 
                               exon.start )

            lcds += len(exon_seq)
            last_end = exon.end

            if E.global_options.loglevel >= 8:
                print "%i: exon_indels (%i-%i):" % (allele_id, exon.start, exon.end)
                for x,c in enumerate(exons[exon_key]):
                    if len(c) != 1: print x + exon.start, ":%s:" % c 
                print
                print exons[exon_key]
                print "genome_pos=", genome_pos, \
                    ",exon=%i-%i" % (genome_pos, genome_pos + len(exon_seq)), \
                    ", len(exon_seq)=", len(exon_seq), \
                    ", len(exon)=", exon.end - exon.start, \
                    ", offsets=%i,%i," % (start_offset, end_offset), \
                    ", offset at start=",_getOffset( exon.start, offsets), \
                    ", offset at end=",_getOffset( exon.end, offsets)

            genome_pos += len(exon_seq)

            test_offset =_getOffset( exon.end, offsets )
            is_offset = genome_pos - exon.end
            assert is_offset == test_offset, "exon offset difference: %i != %i" % (is_offset, test_offset)

        cds = "".join( cds )
        assert lcds == len(cds)

        # fix incomplete codons at the end of the sequence
        if lcds % 3 != 0:
            offset = lcds % 3
            cds=cds[:-offset]
        
        # add frame correction for transcripts that do not start at frame=0
        start_frame = (3 - (int(transcript[0].frame) % 3)) % 3

        # n are ignored (? in sequence to deal with genes like Muc2)
        peptide = Genomics.translate( "n" * start_frame + cds, 
                                      is_seleno = is_seleno, 
                                      prefer_lowercase = False,
                                      ignore_n = True )

        # find the first stop codon
        if start_frame != 0:
            # ignore first, potentially incomplete base
            pep_first_stop = peptide.upper().find( "X", 1 )
        else:
            pep_first_stop = peptide.upper().find( "X" )

        E.debug("%s: translated peptide = %s, first stop at %i" % (transcript_id, peptide, pep_first_stop))

        peptide = peptide.replace( "?", "x" )

        if E.global_options.loglevel >= 8:
            E.debug( "peptide=%s" % peptide )
            E.debug( "cds=%s" % cds )

        E.debug( "%s: start_frame=%i, first stop at %i/%i" % (transcript_id, 
                                                              start_frame,
                                                              pep_first_stop, 
                                                              len(peptide)))


        lpeptide, lcds = len(peptide), len(cds)

        # check for non-sense mediated decay
        if pep_first_stop != -1:
            cds_first_stop = pep_first_stop * 3 - start_frame
            if cds_first_stop < cds_starts[-1]:
                if ncorrected_frameshifts or nuncorrected_frameshifts:
                    E.warn( "nmd knockout transcript %s has frameshifts: %i corrected, %i uncorrected" %\
                                ( transcript_id,
                                  ncorrected_frameshifts, 
                                  nuncorrected_frameshifts) )
                is_nmd_knockout = True
                cds = peptide = ""
                lpeptide, lcds = 0, 0
                reference_first_stop_start, reference_first_stop_end = \
                    (map_cds2reference.mapRowToCol( cds_first_stop ),
                     map_cds2reference.mapRowToCol( cds_first_stop + 3 ) )
            elif pep_first_stop < len(peptide) - 1:
                is_stop_truncated = True
                cds = cds[:cds_first_stop]
                peptide[:pep_first_stop]
                lpeptide, lcds = len(peptide), len(cds)
                reference_first_stop_start, reference_first_stop_end = \
                    (map_cds2reference.mapRowToCol( cds_first_stop ),
                     map_cds2reference.mapRowToCol( cds_first_stop + 3 ) )
            else:
                E.warn( "first stop at %i(cds=%i) ignored: last exon start at %i" % \
                            (pep_first_stop, 
                             cds_first_stop,
                             cds_starts[-1]) )

        else:
            # -1 for no stop codon found
            pep_first_stop = -1
            cds_first_stop = -1
            lpeptide, lcds = len(peptide), len(cds)

        if peptide == None and nframeshifts == 0:
            E.warn( "transcript %s is knockout, though there are no indels - must be nonsense mutation" % (transcript_id))

        # build frames
        frames = [ start_frame ]
        start = start_frame
        l = 0 
        for end in cds_starts[1:]:
            l += end - start
            frames.append( (3 - l % 3) % 3 )
            start = end

        return Allele._make( (cds, 
                              peptide, 
                              len(cds_starts),
                              cds_starts,
                              exon_starts,
                              frames,
                              is_nmd_knockout,
                              is_splice_truncated,
                              is_stop_truncated,
                              nframeshifts,
                              ncorrected_frameshifts,
                              nuncorrected_frameshifts,
                              pep_first_stop,
                              lpeptide,
                              cds_first_stop,
                              lcds,
                              reference_first_stop_start,
                              reference_first_stop_end,
                              loriginal,
                              nsplice_noncanonical,
                              )), map_cds2reference

    if variant_exons or variant_introns:
        for allele in range(0,2):
            exons = dict( [(x,y[allele]) for x,y in variant_exons.iteritems() ] )
            introns = dict( [(x,y[allele]) for x,y in variant_introns.iteritems() ] )
            result.append( _buildAllele( allele, transcript, exons, introns, offsets[allele] ) )
    else:
        a = _buildAllele( 0, transcript, reference_exons, reference_introns, [] )
        result.append( a )
        result.append( a )

    return result



def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: gtf2alleles.py 2886 2010-04-07 08:47:46Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default]."  )
    parser.add_option("-t", "--tablename", dest="tablename", type="string",
                      help="tablename to get variants from (in samtools pileup format) [default=%default]."  )
    parser.add_option("-d", "--database", dest="database", type="string",
                      help="sqlite3 database [default=%default]."  )
    parser.add_option("-f", "--filename-exons", dest="filename_exons", type="string",
                      help="filename with transcript model information (gtf formatted file)  [default=%default]."  )
    parser.add_option("-r", "--filename-reference", dest="filename_reference", type="string",
                      help="filename with transcript models of a reference gene set. Stop codons that do not"
                      " overlap any of the exons in this file are ignore (gtf-formatted file)  [default=%default]."  )
    parser.add_option( "--filename-vcf", dest="filename_vcf", type="string",
                      help="filename with variants in VCF format. Should be indexed by tabix  [default=%default]."  )
    parser.add_option( "--filename-pileup", dest="filename_pileup", type="string",
                      help="filename with variants in samtools pileup format. Should be indexed by tabix  [default=%default]."  )
    parser.add_option( "--vcf-sample", dest="vcf_sample", type="string",
                      help="sample id for species of interest in vcf formatted file [default=%default]."  )
    parser.add_option("-s", "--filename-seleno", dest="filename_seleno", type="string",
                      help="filename of a list of transcript ids that are selenoproteins [default=%default]."  )
    parser.add_option("-m", "--module", dest="modules", type="choice", action="append",
                      choices=("gene-counts", "transcript-effects"),
                      help="modules to apply [default=%default]."  )
    parser.add_option("-o", "--output", dest="output", type="choice", action="append",
                      choices=("all", "peptide", "cds", "table", "gtf", "map"),
                      help="sections to output [default=%default]."  )
    parser.add_option("-k", "--with-knockouts", dest="with_knockouts", action="store_true",
                      help="add alleles that are knocked out to fasta and gtf files [default=%default]." )

    parser.set_defaults(
            genome_file = None,
            filename_exons = None,
            filename_referenec = None,
            filename_seleno = None,
            modules = [],
            border = 200,
            separator = "|",
            tablename = None,
            database = "csvdb",
            output = [],
            with_knockouts = False,
            filename_vcf = None,
            vcf_sample = None,
            )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

    ninput, nskipped, noutput = 0, 0, 0

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )
    else:
        fasta = None

    if options.filename_seleno:
        seleno = set( IOTools.readList( open(options.filename_seleno, "r")) )
    else:
        seleno = {}

    infile_gtf = GTF.gene_iterator( GTF.iterator( options.stdin) )

    # acquire variants from SQLlite database
    if options.tablename:
        if not options.database:
            raise ValueError("please supply both database and tablename")
        variant_getter = VariantGetterSqlite( options.database, options.tablename )
    elif options.filename_pileup:
        variant_getter = VariantGetterPileup( options.filename_pileup )
    elif options.filename_vcf:
        variant_getter = VariantGetterVCF( options.filename_vcf, options.vcf_sample )
    else:
        raise ValueError("please specify a source of variants." )

    if len(options.output) == 0 or "all" in options.output:
        output_all = True
    else:
        output_all = False

    if "cds" in options.output or output_all:
        outfile_cds = E.openOutputFile( "cds.fasta" )
    else: outfile_cds = None

    if "map" in options.output or output_all:
        outfile_map = E.openOutputFile( "map.psl" )
    else: outfile_map = None

    if "peptide" in options.output or output_all:
        outfile_peptides = E.openOutputFile( "peptides.fasta" )
    else: outfile_peptides = None

    if "table" in options.output or output_all:
        outfile_alleles = E.openOutputFile( "table" )
        outfile_alleles.write( "\t".join( 
            ("gene_id", 
             "transcript_id", "allele_id", "contig", "strand", 
             "is_wildtype", 
             ("\t".join( Allele._fields)))) + "\n" )
    else: outfile_alleles = None

    if "gtf" in options.output or output_all:
        outfile_gtf = E.openOutputFile( "gtf" )
    else: outfile_gtf = None

    ## id separatar
    separator = options.separator

    for transcripts in infile_gtf:

        gene_id = transcripts[0][0].gene_id

        overall_start = min( [min( [x.start for x in y]) for y in transcripts] )
        overall_end = max( [max( [x.end for x in y]) for y in transcripts] )
        contig = transcripts[0][0].contig
        strand = transcripts[0][0].strand
        is_positive_strand = Genomics.IsPositiveStrand( strand )
        lcontig = fasta.getLength( contig )
        E.info( "%s: started processing on %s:%i..%i (%s)" % \
                    (gene_id, contig, overall_start, overall_end, strand ))

        ninput += 1
        extended_start = max( 0, overall_start - options.border )
        extended_end = min( lcontig, overall_end + options.border )

        # if contig.startswith("chr"): contig = contig[3:]

        variants = variant_getter( contig, extended_start, extended_end )

        E.debug( "%s: found %i variants in %s:%i..%i" % \
                     (gene_id, len(variants), contig, extended_start, extended_end ))

        if E.global_options.loglevel >= 10:
            print "# collected variants:", variants

        # collect intron/exon sequences
        # coordinates are forward/reverse
        # also updates the coordinates in transcripts
        all_exons, all_introns = collectExonIntronSequences( transcripts, fasta )

        # update variants such that they use the same coordinates
        # as the transcript
        variants = Variants.updateVariants( variants, lcontig, strand )

        # deal with overlapping but consistent variants
        variants = Variants.mergeVariants( variants )

        E.debug( "%s: found %i variants after merging in %s:%i..%i" % \
                     (gene_id, len(variants), contig, extended_start, extended_end ))

        if E.global_options.loglevel >= 10:
            print "# merged variants:", variants

        # collect coordinate offsets and remove conflicting variants
        variants, removed_variants, offsets = Variants.buildOffsets( variants, contig = contig ) 

        if len(removed_variants) > 0:
            E.warn("removed %i conflicting variants" % len(removed_variants) )
            for v in removed_variants:
                E.info("removed variant: %s" % str(v))

        E.info( "%i variants after filtering" % len(variants) )

        if len(variants) > 0:
            # build variants
            indexed_variants = Variants.indexVariants( variants )

            # update exon sequences according to variants
            variant_exons = buildVariantSequences( indexed_variants, all_exons )

            # update intron sequences according to variants
            variant_introns = buildVariantSequences( indexed_variants, all_introns )

            if E.global_options.loglevel >= 10:
                for key in variant_exons:
                    print "exon", key
                    Genomics.printPrettyAlignment( 
                        all_exons[key],
                        variant_exons[key][0],
                        variant_exons[key][1],
                        )
                for key in variant_introns:
                    print "intron", key
                    Genomics.printPrettyAlignment( 
                        all_introns[key][:30] + all_introns[key][-30:],
                        variant_introns[key][0][:30] + variant_introns[key][0][-30:],
                        variant_introns[key][1][:30] + variant_introns[key][1][-30:] )
        
        else:
            variant_exons, variant_introns = None, None


        for transcript in transcripts:

            transcript.sort( key = lambda x: x.start)

            transcript_id = transcript[0].transcript_id
            alleles = buildAlleles( transcript, 
                                    variant_exons, 
                                    variant_introns,
                                    all_exons,
                                    all_introns,
                                    offsets,
                                    is_seleno = transcript_id in seleno,
                                    reference_coordinates = False,
                                    )

            ##############################################################
            ##############################################################
            ##############################################################
            # output
            for aid, al in enumerate(alleles):

                allele, map_cds2reference = al

                reference_cds_sequence = buildCDSSequence( transcript, all_exons )
                is_wildtype = reference_cds_sequence == allele.cds

                allele_id = str(aid)
                assert len(allele.exon_starts) == allele.nexons
                assert len(allele.cds_starts) == allele.nexons
                assert len(allele.frames) == allele.nexons

                # the output id
                outid = separator.join( (gene_id, transcript_id, allele_id) )

                # output map between cds and reference
                if outfile_map and map_cds2reference:
                    match = Blat.Match()
                    match.mQueryId = allele_id
                    match.mQueryLength = allele.cds_len
                    match.mSbjctId = contig
                    match.mSbjctLength = lcontig
                    match.strand = strand
                    match.fromMap( map_cds2reference, use_strand = True )
                    outfile_map.write( "%s\n" % str(match) )
                    
                # only output sequences for genes that have not been knocked out, unless required
                if not allele.is_nmd_knockout or options.with_knockouts:

                    if outfile_gtf:
                        gtf = GTF.Entry()
                        gtf.gene_id = gene_id
                        gtf.transcript_id = transcript_id
                        gtf.addAttribute( "allele_id", allele_id)
                        gtf.contig = contig
                        gtf.strand = strand
                        gtf.feature = "CDS"
                        gtf.source = "gtfxnsps"
                        l = 0
                        last_cds_start = allele.cds_starts[0]
                        gtf.start = allele.exon_starts[0]
                        gtf.frame = allele.frames[0]

                        for exon_start, cds_start, frame in zip( allele.exon_starts[1:], 
                                                                 allele.cds_starts[1:], 
                                                                 allele.frames[1:]):
                            cds_length = cds_start - last_cds_start
                            gtf.end = gtf.start + cds_length
                            if not is_positive_strand:
                                gtf.start, gtf.end = lcontig - gtf.end, lcontig - gtf.start
                            outfile_gtf.write( str(gtf) + "\n" )

                            gtf.start = exon_start
                            gtf.frame = frame

                            l += cds_length
                            last_cds_start = cds_start

                        cds_length = len(allele.cds) - last_cds_start
                        gtf.end = gtf.start + cds_length
                        if not is_positive_strand:
                            gtf.start, gtf.end = lcontig - gtf.end, lcontig - gtf.start
                        outfile_gtf.write( str(gtf) + "\n" )

                    if outfile_cds:
                        outfile_cds.write(">%s\n%s\n" % (outid, allele.cds ))
                    if outfile_peptides:
                        outfile_peptides.write(">%s\n%s\n" % (outid, allele.peptide ))

                # reformat for tabular output
                allele = allele._replace( 
                    cds_starts = ",".join(map(str,allele.cds_starts)),
                    exon_starts = ",".join(map(str,allele.exon_starts)),
                    frames = ",".join(map(str,allele.frames)) )
                
                # convert reference coordinates to positive strand coordinates
                if allele.reference_first_stop_start >= 0 and not is_positive_strand:
                    allele = allele._replace(
                        reference_first_stop_start = lcontig - allele.reference_first_stop_end,
                        reference_first_stop_end = lcontig - allele.reference_first_stop_start, )
                    
                if outfile_alleles:
                    outfile_alleles.write( "%s\t%s\n" % (\
                            "\t".join( (gene_id, 
                                        transcript_id, 
                                        allele_id,
                                        contig,
                                        strand,
                                        "%i" % is_wildtype) ),
                            "\t".join(map(str,allele)) ) )
                            
                noutput += 1
                # only output first allele (debugging)
                # break

    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput,nskipped) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )










