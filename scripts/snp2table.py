################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: snp2table.py 2861 2010-02-23 17:36:32Z andreas $
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
snp2table.py - annotate variants
================================

:Author: Andreas Heger
:Release: $Id: snp2table.py 2861 2010-02-23 17:36:32Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

This script reads a file in samtools pileup -c format and annotates the
SNPs in that file.

.. note::

   This script is still under construction. There are issues with the
   phasing of variants, using the wild-type sequence and merging
   splice and coding variants.

Usage
-----

Example::

   gunzip < snps.pileup.gz 
   | python snp2table.py  
        --genome-file=genome 
        --filename-exons=exons
        --filename-annotations=bases
        --filename-junctions=junctions
        --log=log 
   > result.out

Type::

   python <script_name>.py --help

for command line help.

Documenation
------------

The following annotations are available:

exons
   Given a gtf-file with the output from ``gtf2gff.py --method=exons``, each
   SNP that overlaps an exon is labeled according to its position within a transcript. 
   The columns are:

   ntranscripts
      number of transrcipts in the gene that the exon belongs to.

   nused
      number of transcripts that contain this exon 
   
   positions
      positions of exons within transcipts

annotation
   Given a fasta file with base-level annotations (output from ``gtf2fasta.py``),
   each SNP is labeled by its effect. The columns are:

   code
      one letter annotation of base. See :file:`gtf2fasta.py`.
   reference_codon
      the reference codon (in case sequence is coding)
   reference_aa
      the reference amina acid (in case sequence is coding)
   variant_type
      the variant types ('I' for insertion, 'D' for deletion', 'S' for substitution).
   variant_codon
      the variant codons (in case sequence is coding)
   variant_aa
      the variant amino acids (in case sequence is coding)

   Variants will be separated by ``,`` if there are several variants per position.

The script deconvolutes substitution and indel snps. The following columns are always output:

chromosome
   chromosome
position
   reference base position  (0-based)
reference_base
   reference base. For indels, this is ``*``.
genotype
   consensus base. For indels, this is the genotype
consensus_quality
   consensus quality
snp_quality
   snp quality
rms_mapping_qualit
   rms quality
coverage
   coverage
read_bases
   read bases. For indels, this is the first allele (``*`` for reference).
base_qualities
   base qualities. For indels, this is the second allele (``*`` for reference).

See the samtools documentation for the meaning of these.

Code
----

""" 

import os
import sys
import re
import optparse
import collections

import CGAT.Experiment as E

import pysam
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IndexedGenome as IndexedGenome
import CGAT.Genomics as Genomics

def readJunctions( filename_junctions ):
    '''read junctions from a tab-separated file.

    returns a dictionary of dictionaries for splice junctions per
    contig and position.
    '''
    
    junctions = collections.defaultdict( dict )
    infile = open( filename_junctions, "r" )
    njunctions = 0

    for line in infile:
        if line.startswith( "#" ): continue
        data = line[:-1].split("\t")
        if data[0] == "contig": continue
        contig, strand, end, start, frame = data[:5]
        start,end,frame = int(start), int(end), int(frame)
        j = junctions[contig]
        if end in j:
            j[end].add( (strand,start,frame) )
        else:
            j[end] = set( ((strand,start,frame),) )
        njunctions += 1

    infile.close()
    E.info( "read %i junctions for %i contigs" % (njunctions, len(junctions) ) )
    return junctions

class BaseAnnotator(object):
    '''annotator for single bases in the genome.'''

    mHeader = ()

    def __init__(self, fasta = None, *args, **kwargs):
        self.mFasta = fasta

    def __str__( self ):
        return ""

    def getHeader( self ):
        '''return header'''
        return "\t".join(self.mHeader)

class BaseAnnotatorSNP( BaseAnnotator ):
    '''output SNPs. 
    
    This annotator deconvolutes substitution and indel SNPs.
    '''

    mHeader = (
        "chromosome", 
        "position", 
        "reference_base", 
        "genotype",
        "consensus_quality",
        "snp_quality",
        "rms_mapping_quality",
        "coverage" )

    mAdditionalHeader = (
        "read_bases",
        "base_qualities" )

    def __init__(self, *args, **kwargs):
        BaseAnnotator.__init__(self, *args,**kwargs)
    
    def __str__(self ):
        # truncate the last two columns to make the snp output even length
        return "\t".join( map(str, self.mSNP)[:len(self.mHeader)] )

    def update( self, snp ):
        '''update with snp.'''
        self.mSNP = snp

class BaseAnnotatorExons( BaseAnnotator ):
    '''annotate SNP by exons that it overlaps with.'''

    mHeader = [ "exons_%s" % x for x in ( "ntranscripts", "nused", "pos" )]

    def __init__(self, filename_exons, *args, **kwargs):
        BaseAnnotator.__init__(self, *args,**kwargs)

        exons = IndexedGenome.IndexedGenome()
        nexons = 0
        for g in GTF.iterator( open( filename_exons, "r") ):
            exons.add( g.contig, g.start, g.end, g )
            nexons += 1

        self.mExons = exons

        E.info( "indexed %i exons on %i contigs" % (nexons, len(exons) ) )

    def update( self, snp ):
        '''update with snp.'''

        exons = list(self.mExons.get( snp.chromosome, 
                                      snp.pos, 
                                      snp.pos+1))

        if exons:
            for start, end, gff in exons:
                self.mOTranscripts = "%i" % gff["ntranscripts"]
                self.mOUsed = "%i" % gff["nused"]
                self.mOPos = "%s" % gff["pos"]
        else:
            self.mOTranscripts = "na"
            self.mOUsed = "na"
            self.mOPos = "na"

    def __str__(self):
        return "\t".join( (self.mOTranscripts, self.mOUsed, self.mOPos) )

class BaseAnnotatorSpliceSites( BaseAnnotator ):
    '''annotate SNP by splice sites that it overlaps with.

    Splice sites are read from filename_junctions.

    The module is aware of the following splice signals
    
       * GT/AG: A
       * GC/AG: B
       * AT/AC: C
       
    A SNP is labeled by splice site if it is directly overlapping.
    Furthermore, it will be labeled as synonymous (S), non-synomous (N,
    if it changes from one splice class to another) or as a deleterious (D).

    In particular, the codes are:
    
+---------------+----------------------------------------------------------------------+                                                                                                                                                                   
|code           | description                                                          |
+---------------+----------------------------------------------------------------------+                                                                                                                                                                   
|S              | synonymous SNP - the splice signal is unchanged                      |
+---------------+----------------------------------------------------------------------+
|N              | non-synonymous SNP - one splice site changes to another              |
+---------------+----------------------------------------------------------------------+
|D              | deleterious snp, a canonical  splice site is abrogated               |
+---------------+----------------------------------------------------------------------+
|U              | change from unknown to unknown splice                                |
+---------------+----------------------------------------------------------------------+
|C              | creation of an new canonical splice site                             |
+---------------+----------------------------------------------------------------------+


    All SNPS will be annotatoted if they are within 10 bases of a splice
    site. Furthermore, the overlap flag will be set for those
    SNPs overlapping one of the known splice signals.

    Note that a disrupted splice signal might still lead to a valid
    protein product.
    '''

    mHeader = [ "splice_%s" % x \
                for x in ( "transcript_ids", "gene_ids", "has_overlap", "is_homozygous", "old", "new", "distance", "code" )]
    mIntronTypes = ( ("U2-GT/AG", "GT", "AG" ),
                     ("U2-nc-GC/AG", "GC", "AG" ),
                     ("U12-AT/AC", "AT", "AC") )
    # size of splice signal
    mSize = 10

    def __init__(self, filename_junctions, *args, **kwargs):
        BaseAnnotator.__init__(self, *args,**kwargs)

        junctions = IndexedGenome.IndexedGenome()
        
        infile = open( filename_junctions, "r" )
        njunctions = 0

        for line in infile:
            if line.startswith( "#" ): continue
            data = line[:-1].split("\t")
            if data[0] == "contig": continue
            # end, start are the positions of the last base of the codon
            # 5' of the intron and first base of codon 3' of the intron.
            contig, strand, end, start, frame, gene_id, transcript_id = data
            start,end,frame = int(start), int(end), int(frame)
            # convert to intron coordinates
            intron_start, intron_end = end + 1, start
            # convert to positive strand coordinates
            if strand == "-":
                lcontig = self.mFasta.getLength( contig )
                intron_start,intron_end = lcontig-intron_end, lcontig-intron_start
            junctions.add( contig, intron_start, intron_start+self.mSize, (strand, intron_start, intron_end, gene_id, transcript_id) )
            junctions.add( contig, intron_end-self.mSize, intron_end, (strand, intron_start, intron_end, gene_id, transcript_id) )
            njunctions += 1
        infile.close()
        
        self.mJunctions = junctions
        
        E.info( "read and indexed %i junctions for %i contigs" % (njunctions, len(junctions) ) )

    def hasSpliceMotif( self, seq5, seq3 ):
        """find a splice motif in seq.

        returns the name of the splice motif or None
        """
        for name, prime5, prime3 in self.mIntronTypes:
            if seq5.startswith( prime5 ) and seq3.endswith( prime3 ):
                return name, prime5, prime3
        return None, None, None

    def buildSequenceVariants( self, seq, strand, pos, snp ):
        '''build new sequence by modifying a sequence fragment in seq at
        pos with snp.

        It is assumed that seq is already oriented according to strand.
        The strand is used to revert the snp if necessary.

        Note that only sequences different from seq will be returned.

        returns is_homozygous, seqs
        '''
        is_negative_strand = Genomics.IsNegativeStrand( strand )
        reference_base = snp.reference_base

        if reference_base != "*" and is_negative_strand:
            reference_base = Genomics.complement( reference_base )

        new_sequences = []
        is_homozygous = True
        if reference_base != "*":
            if seq[pos].upper() != reference_base.upper():
                raise ValueError( "base mismatch at snp %i, expected %s, got %s in %s at position %i; snp=%s" % \
                                  (snp.pos, reference_base, seq[pos], seq, pos,
                                   ";".join(map(str, snp))) )

            # single base changes
            variant_bases = Genomics.resolveAmbiguousNA( snp.genotype )
            if len(variant_bases) == 1:
                is_homozygous = True
            else:
                is_homozygous = False

            for variant_base in variant_bases:
                if is_negative_strand:
                    variant_base = Genomics.complement( variant_base )

                s = list(seq)
                s[pos] = variant_base
                s = "".join(s)
                if s != seq: new_sequences.append( s )
        else:
            variants = snp.genotype.split("/")
            is_homozygous = False
            for variant in variants:
                
                s = list(seq)
                # samtools denotes insert/deletion after position
                # while python is before/at position, hence the pos+1
                if variant[0] == "+":
                    toinsert = variant[1:].upper()
                    if is_negative_strand:
                        toinsert = Genomics.complement( toinsert )
                        s.insert( pos, toinsert )
                    else:
                        s.insert( pos+1, toinsert )

                elif variant[0] == "-":
                    # pos+1+len(x)-1 = pos+len(x)
                    todelete = variant[1:].upper()
                    l = len(todelete)
                    if is_negative_strand:
                        # delete left of pos
                        xstart = max(0, pos-l)
                        xend = pos
                        todelete = todelete[:min(l,pos)]
                    else:
                        # delete right of pos
                        xstart = pos+1
                        xend = min(self.mSize, pos + 1 + l)
                        todelete = todelete[:self.mSize-(pos+1)]
                            
                    deleted = "".join(s[xstart:xend])
                    
                    if is_negative_strand:
                        deleted = Genomics.complement( deleted )

                    if deleted != todelete:
                        raise ValueError( "base mismatch at indel %i, expected %s, got %s in %s at position %i(%i:%i); is_negative_strand=%s, snp=%s" % \
                                (snp.pos, todelete, deleted, seq, pos, xstart, xend,
                                 is_negative_strand,
                                 ";".join(map(str, snp))) )
                    del s[xstart:xend]
                    
                elif variant[0] == "*":
                    is_homozygous = True
                else:
                    raise ValueError( "unknown variant sign '%s'" % variant[0] )

                s = "".join(s)
                if s != seq: new_sequences.append( s )
            
        return is_homozygous, new_sequences

    def update( self, snp ):
        '''update with snp.'''

        junctions = list(self.mJunctions.get( snp.chromosome, 
                                              snp.pos, 
                                              snp.pos+1))

        self.mResults = []

        reference_base = snp.reference_base
        
        if junctions:
            for start, end, data in junctions:
                strand, intron_start, intron_end, gene_id, transcript_id = data
                snp_pos = snp.pos
                
                if strand == "-":
                    lcontig = self.mFasta.getLength( snp.chromosome )
                    intron_start, intron_end = lcontig - intron_end, lcontig - intron_start
                    snp_pos = lcontig - snp_pos - 1

                seq5 = self.mFasta.getSequence( snp.chromosome, strand, intron_start, intron_start + self.mSize ).upper()
                seq3 = self.mFasta.getSequence( snp.chromosome, strand, intron_end - self.mSize, intron_end ).upper()

                splice_name, prime5, prime3 = self.hasSpliceMotif( seq5, seq3 )
                if splice_name:
                    if intron_start <= snp_pos < intron_start + len(prime5):
                        has_overlap = "5"
                    elif intron_end - len(prime3) <= snp_pos < intron_end:
                        has_overlap = "3"
                    else:
                        has_overlap = ""
                else:
                    has_overlap = ""
                    
                E.debug( "building a splice site for %s:%s:%i, splice=%i:%i:%s (snp=%i), splice_name=%s, has_overlap=%s" %\
                         (snp.chromosome, strand, snp.pos,
                          intron_start, intron_end, strand, snp_pos,
                          str(splice_name),
                          str(has_overlap) ))

                if intron_start <= snp_pos < intron_start + self.mSize:
                    is_homozygous, new_seq5 = self.buildSequenceVariants( seq5, strand, snp_pos - intron_start, snp )
                    new_seqs = [ (x, seq3) for x in new_seq5 ]
                    distance = snp_pos - intron_start
                elif intron_end - self.mSize <= snp_pos < intron_end:
                    is_homozygous, new_seq3 = self.buildSequenceVariants( seq3, strand, self.mSize - (intron_end - snp_pos), snp )
                    new_seqs = [ (seq5, x) for x in new_seq3 ]
                    distance = intron_end - snp_pos - 1
                else:
                    raise ValueError( "no overlap?")

                E.debug( "building a splice site for %s:%s:%i, splice=%i:%i:%s, splice_name=%s, distance=%i, has_overlap=%s, is_homozygous=%s, new_seqs=%s" %\
                         (snp.chromosome, strand, snp.pos,
                          intron_start, intron_end, strand,
                          str(splice_name),
                          distance,
                          str(has_overlap),
                          str(is_homozygous),
                          str(new_seqs) ) )

                # indels at the end of size are out-of-range - ignore
                if new_seqs == []: continue

                for new_seq5, new_seq3 in new_seqs:
                    new_splice_name, new_prime5, new_prime3 = self.hasSpliceMotif( new_seq5, new_seq3 )

                    if splice_name == None and new_splice_name == None:
                        code = "U"
                    elif new_splice_name == None:
                        code = "D"
                    elif splice_name == None:
                        code = "C"
                    elif splice_name == new_splice_name:
                        code = "S"
                    elif splice_name != new_splice_name:
                        code = "N"

                    if is_homozygous: is_homozygous = "1"
                    else: is_homozygous = "0"
                    if not splice_name: splice_name = ""
                    if not new_splice_name: new_splice_name = ""
                    self.mResults.append( (transcript_id, gene_id, has_overlap, is_homozygous,
                                           splice_name, new_splice_name, str(distance), code ) )

    def __str__(self):
        if self.mResults:
            return "\t".join( [ ";".join( x ) for x in zip( *self.mResults ) ] )
        else:
            return "\t" * 8
        
class BaseAnnotatorCodon( BaseAnnotator ):
    '''annotate snp by coding sequence that it overlaps with

    Variant types are annotated with a one-letter code.

+---------------+----------------------------------------------------------------------+                                                                                                                                                                   
|code           | description                                                          |
+---------------+----------------------------------------------------------------------+                                                                                                                                                                   
|O              | homozygous substitution                                              |                                                                                                                                                                   
+---------------+----------------------------------------------------------------------+                                                                                                                                                                   
|E              | heterozygous substitution                                            |
+---------------+----------------------------------------------------------------------+
|I              | insertion with respect to reference genome                           |
+---------------+----------------------------------------------------------------------+
|D              | deletion with respect to reference genome                            |
+---------------+----------------------------------------------------------------------+
|W              | wild type (reference genome)                                         |
+---------------+----------------------------------------------------------------------+

    This function will not annotate indels accurately.
    '''

    mHeader = [ "code", "reference_codon", "reference_aa", "variant_type", "variant_codon", "variant_aa" ]

    def __init__(self, annotations_file, junctions, *args, **kwargs):

        BaseAnnotator.__init__(self, *args,**kwargs)

        self.mAnnotations = IndexedFasta.IndexedFasta( annotations_file )
        self.mJunctions = junctions

    def updateSNPs( self, snp, is_negative_strand, pos ):
        '''update SNPs.'''

        contig = snp.chromosome
        lcontig = self.mFasta.getLength( contig )
        reference_base = snp.reference_base

        if snp.genotype in 'ACGTacgt':
            # homozygous substitution
            self.mVariantType.append( "O" )
        else:
            # heterozygous substitution
            self.mVariantType.append( "E" )

        # switch reference strand codon to correct strand
        if reference_base != "*" and is_negative_strand:
            reference_base = Genomics.complement( reference_base )

        # collect all possible variants of reference codons
        for reference_codon in self.mReferenceCodons:

            self.mReferenceAAs.append( Genomics.translate( reference_codon ) )

            # process single base changes
            variant_bases = Genomics.resolveAmbiguousNA( snp.genotype )

            if reference_codon[ pos ] != reference_base:
                raise ValueError( "base mismatch at %i (codon=%s,%i): codon:%s != genome:%s; `%s`" % \
                                      (snp.pos, reference_codon, pos, reference_codon[pos], reference_base, ";".join(map(str, snp))) )

            for variant_base in variant_bases:
                if is_negative_strand:
                    variant_base = Genomics.complement( variant_base )

        self.mVariantAAs.extend( [ Genomics.translate( x ) for x in self.mVariantCodons ] )

    def updateIndels( self, snp, is_negative_strand ):

        contig = snp.chromosome
        lcontig = self.mFasta.getLength( contig )

        # get location of insertion/deletion. The location
        # is after position, hence get position and position + 1
        code = self.mAnnotations.getSequence( contig, "+", snp.pos, snp.pos+2 )
        self.mCode = code

        variants = snp.genotype.split("/")
        for variant in variants:

            if variant[0] == "*":
                self.mVariantType.append( "W" )

            elif variant[0] == "+":
                toinsert = variant[1:] 
                self.mVariantType.append( "I" )
                
            elif variant[0] == "-":
                todelete = variant[1:]
                # deletions need to be looked at in a wider range
                self.mVariantType.append( "D" )
                
            else:
                raise ValueError( "unknown variant sign '%s'" % variant[0] )


        # ignore non-coding Indels
        if code[0] and code[1] not in 'abcABC': return

        if is_negative_strand:
            variants = [ Genomics.complement(x) for x in variants ]

        for reference_codon in self.mReferenceCodons:

            variants = snp.genotype.split("/")
            variants = [ x[1:] for x in variants ]

            for variant in variants:
                if len(variant) % 3 != 0:
                    self.mVariantCodons.append( "!" )
                else:
                    self.mVariantCodons.append( variant )

            self.mVariantAAs.extend( [ Genomics.translate( x ) for x in self.mVariantCodons ] )

    def update( self, snp ):
        '''update with snp.'''

        contig = snp.chromosome
        lcontig = self.mFasta.getLength( contig )

        # check for codon:
        self.mReferenceCodons = []
        self.mVariantCodons = []
        self.mReferenceAAs = []
        self.mVariantAAs = []
        self.mVariantType = []

        if snp.pos >= lcontig:
            E.warn("snp is out of range: %s:%i => %s" % (contig, lcontig, ";".join( map(str, snp)) ))
            self.mCode = "X"
            return

        reference_base = snp.reference_base

        ######
        code = self.mAnnotations.getSequence( contig, "+", snp.pos, snp.pos+1 )
        self.mCode = code

        is_negative_strand = False

        pos = None
        # ignore non-coding SNPs
        if code in 'abcABC': 

            pos = "abcABC".find(code)
            is_negative_strand = code in 'ABC'

            if is_negative_strand:
                pos -= 3
                position = lcontig - snp.pos - 1
                strand = "-"
            else:
                position = snp.pos
                strand = "+"

            # start, end are forward/negative strand coordinates
            offset = pos % 3 
            start, end = position - offset, position - offset + 3

            # If a codon is split by an intron, several 
            # reference codons are possible. 
            self.mReferenceCodons = set()
            if self.mJunctions:
                if contig in self.mJunctions:
                    j = self.mJunctions[contig]
                    for x in xrange(start, end):
                        if x in j:
                            # found a junction within codon
                            for junction_strand, next_start, frame in j[x]:

                                # wrong strand
                                if (junction_strand == "-") != is_negative_strand:
                                    continue

                                # codon is complete, junction within frame
                                if frame == 0:
                                    self.mReferenceCodons.add( self.mFasta.getSequence( contig, strand, start, end ).upper() )                                    
                                    continue

                                codon = self.mFasta.getSequence( contig, strand, start, x + 1) + \
                                    self.mFasta.getSequence( contig, strand, next_start, next_start + end - (x+1))

                                # position within codon
                                y = x + 1 - start
                                E.debug( "building a split codon for %s:%s:%i(%i) start=%i, split=%i, cont=%i, end=%i: codon=%s:%s" % \
                                         (contig, strand, snp.pos, position, start, x+1, next_start, next_start + end - (x+1),
                                         codon[:y], codon[y:]))

                                assert len(codon) == 3
                                self.mReferenceCodons.add( codon.upper() )

            # ordinary codon                
            if not self.mReferenceCodons:
                self.mReferenceCodons.add( self.mFasta.getSequence( contig, strand, start, end ).upper() )

            # fixate the order of reference codons
            self.mReferenceCodons = tuple( self.mReferenceCodons )

        # process according to variant type
        # indels need to be treated differently from SNPs as
        # they have larger effects
        if reference_base == "*":
            self.updateIndels( snp, is_negative_strand )
        else:
            self.updateSNPs( snp, is_negative_strand, pos )

    def __str__(self):
        return "\t".join( (self.mCode, 
                           ",".join(self.mReferenceCodons), 
                           ",".join(self.mReferenceAAs),
                           ",".join(self.mVariantType), 
                           ",".join(self.mVariantCodons), 
                           ",".join(self.mVariantAAs), ) )

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: snp2table.py 2861 2010-02-23 17:36:32Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default]."  )
    parser.add_option("-a", "--filename-annotations", dest="filename_annotations", type="string",
                      help="filename with base annotations (output from gtf2fasta.py) [default=%default]."  )
    parser.add_option("-f", "--filename-exons", dest="filename_exons", type="string",
                      help="filename with exon information (gff formatted file)  [default=%default]."  )
    parser.add_option("-j", "--filename-junctions", dest="filename_junctions", type="string",
                      help="filename with junction information (filename with exon junctions)  [default=%default]."  )
    parser.add_option("-c", "--filename-vcf", dest="filename_vcf", type="string",
                      help="vcf file to parse [default=%default]."  )
    parser.add_option("-i", "--input-format", dest="input_format", type="choice",
                      choices = ("pileup", "vcf" ),
                      help="input format [default=%default]."  )
    parser.add_option( "--vcf-sample", dest="vcf_sample", type="string",
                      help="sample id in vcf file to analyse [default=%default]."  )

    parser.set_defaults(
        genome_file = None,
        filename_annotations = None,
        filename_exons = None,
        filename_junctions = None,
        input_format = "pileup",
        vcf_sample = None,
        filename_vcf = None,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    ninput, nskipped, noutput = 0, 0, 0

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )
    else:
        fasta = None

    if options.filename_junctions:
        junctions = readJunctions( options.filename_junctions )
    else:
        junctions = None

    # setup iterator
    if options.input_format == "pileup":
        iterator = pysam.Pileup.iterate(sys.stdin)
    elif options.input_format == "vcf":
        if not options.vcf_sample:
            raise ValueError( "vcf format requires sample id (--vcf-sample) to be set" )
        if not options.filename_vcf:
            raise ValueError( "reading from vcf requires vcf filename (--filename-vcf) to be set)" )

        iterator = pysam.Pileup.iterate_from_vcf(options.filename_vcf, options.vcf_sample )
    

    modules = []
    modules.append( BaseAnnotatorSNP() )

    if options.filename_exons:
        modules.append( BaseAnnotatorExons( options.filename_exons, fasta=fasta ) )
    if options.filename_annotations:
        modules.append( BaseAnnotatorCodon( options.filename_annotations, fasta=fasta, junctions=junctions ) )
    if options.filename_junctions:
        modules.append( BaseAnnotatorSpliceSites( options.filename_junctions, fasta=fasta ) )
                        
    options.stdout.write( "\t".join( [x.getHeader() for x in modules]) + "\n" )

    for snp in iterator:
        ninput += 1

        # translate chromosome according to fasta
        if fasta: 
            try:
                snp = snp._replace( chromosome=fasta.getToken( snp.chromosome ))
            except KeyError:
                E.warn( "unknown contig `%s` for snp `%s`" % (snp.chromosome, str(snp) ) )
                continue

        for module in modules:
            module.update( snp )

        options.stdout.write( "\t".join( map(str, modules) ) + "\n" )

        noutput += 1

    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput,nskipped) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

