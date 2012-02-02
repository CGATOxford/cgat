################################################################################
#
#   Gene prediction pipeline 
#
#   $Id: gtf2gff.py 2861 2010-02-23 17:36:32Z andreas $
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
gtf2gff.py - flatten a gtf file
===============================

:Author: Andreas Heger
:Release: $Id: gtf2gff.py 2861 2010-02-23 17:36:32Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

This scripts converts a :term:`gtf` formatted file into a :term:`gff` formatted file. 
In other words, a gene set (gtf), which constitutes a hierarchal set of annotations, 
will be converted in into a flat list of genomic segments.

Various methods can be used to do the conversion (see command line argument *method*):

genome
   annotate genome with gene set. Genomic segments are labeled 
   ``intronic``, ``intergenic``, ...

territories
   build territories around gene set.

exons
   annotate exons. Exonic segments are classified according 
   to the transcript structure

promoters
   declare promoter regions. These segments might be overlapping. A promotor
   is the region x kb upstream of a transcription start site. The option
   ``--promotor`` sets the region width.

regulons
   declare regulatory regions. Regulatory regions contain the region x kb of
   upstream and downstream of a transciption start site. The option ``--promotor``
   sets the region width.

Genome
++++++

If ``--method=genome``, the gene set is used to annotate the complete genome.

.. note::
   The gtf file has to be sorted first by contig and then by position.

A segment in the genome will either be covered by:

cds
   a coding exon (also: CDS, start_codon).

utr
   a UTR (also: stop_codon) 

5flank, 3flank, flank
   an upstream/downstream segment of defined size. If 
   the intergenic region is too small to accomodate a flank, the regions is
   just 'flank'.
                
intergenic
   intergenic region.

5telomeric, 3telomeric
   telomeric region (before/after first/last gene).

intronic
   intronic region.

frameshift
   frameshift. Introns of less than 10 residues length

ambiguous
   in case of overlapping genes, regions are designated ambiguous

All segments are annotated by their closest gene. Intergenic regions are
annotated with their two neighbouring genes. The upstream gene is listed
in the attribute gene_id, the downstream one is listed in the attribute
downstream_gene_id.

Territories
+++++++++++

If ``--method=territories``, the gene set is used to define gene territories. 
Territories are segments around genes.

.. note::
   The gtf file has to be sorted first by contig and then by position.

.. note::
   Genes should already have been merged (gtf2gtf --merge-transcripts)

Exons
+++++

If ``--method=exons``, exons are annotated by their dispensibility.

.. note::
   The gtf file should be sorted by genes

For each exon, the following additional fields are added to the gtf file:

ntranscripts
   number of transcripts
nused
   number of transcripts using this exon
positions
   positions of exon within transcripts. This is a ``,`` separated list of
   tuples ``pos:total``. For example, ``1:10,5:8`` indicates an exon that appears
   in first position in a ten exon transcript and fifth position in an eight
   exon transcript. The position is according to the direction of transcription.

.. note::
   overlapping but non-identical exons, for example due to internal splice sites, are listed 
   as separate exons. Thus the output is not fully flat as some segments could be overlapping
   (see output variable ``noverlapping`` in the log file).

The following example uses an ENSEMBL gene set::

   gunzip < Mus_musculus.NCBIM37.55.gtf.gz | awk '$3 == "CDS"' | python gtf2gff.py --method=exons --restrict-source=protein_coding

Promoters
+++++++++

If ``--method=promotors``, putative promotor regions are output. A promoter is
a x kb segment upstream of the transcription start site. As the actual start site
is usually not known, the start of the first exon within a transcript is used as a proxy.
A gene can have several promotors associated with it, but overlapping promotor
regions of the same gene will be merged. A promoter can extend into an adjacent
upstream gene.

Only protein coding genes are used to annotate promotors.

The size of the promotor region can be specified by the command line
argument ``--promotor``

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----
"""

import os, sys, string, re, optparse, types, collections

import GFF, GTF
import Experiment as E
import IndexedFasta
import Genomics
import Intervals

def addSegment( feature, start, end, template, options ):
    """add a generic segment of type *feature*.
    """
    if start >= end: return 0

    entry = GTF.Entry()

    if type(template) == types.TupleType:
        entry.copy( template[0] )        
        entry.clearAttributes()
        entry.addAttribute( "downstream_gene_id", template[1].gene_id )
    else:
        entry.copy( template )
        entry.clearAttributes()

    entry.start, entry.end = start, end
    entry.feature = feature
    if feature not in ("exon", "CDS", "UTR", "UTR3", "UTR5"):
        entry.score = "."
    options.stdout.write( str(entry) + "\n" )

    return 1

def addFlank( start, end, template, options ):
    """add a flank.
    """    
    is_positive = Genomics.IsPositiveStrand( template.strand )
    is_before = end <= template.start
    if (is_before and is_positive) or (not is_before and not is_positive):
        name = "5flank"
    else:
        name = "3flank"
            
    return addSegment( name, start, end, template, options )

def addIntergenicSegment( last, this, fasta, options ):
    """add an intergenic segment between last and this.

    At telomeres, either can be None.
    """
    if not this and not last: return 0
        
    nadded = 0
    if not this:
        # last telomere
        try:
            lcontig = fasta.getLength( last.contig )
        except KeyError, msg:
            if options.ignore_missing:
                return nadded
            else:
                raise KeyError, msg
        flank = min( last.end + options.flank, lcontig )
        nadded += addFlank( last.end, flank, last, options )
        nadded += addSegment( "telomeric", flank, lcontig, last, options )
    elif not last:
        # first telomere
        flank = max(0, this.start - options.flank )
        nadded += addSegment( "telomeric", 0, flank, this, options )
        nadded += addFlank( flank, this.start, this, options )
    else:
        # intergenic region
        d = this.start - last.end
        flank = options.flank
        if d > flank * 2:
            nadded += addFlank( last.end, last.end + flank, last, options )
            nadded += addSegment( "intergenic", last.end + flank, this.start - flank, (last,this), options )
            nadded += addFlank( this.start - flank, this.start, this, options )
        else:
            # add short flank between two genes. If they can not agree
            # on the directionality, "flank" is used.
            is_positive1 = Genomics.IsPositiveStrand( last.strand ) 
            is_positive2 = Genomics.IsPositiveStrand( this.strand )
            if is_positive1 and not is_positive2:
                key = "3flank"
            elif not is_positive1 and is_positive2:
                key = "5flank"
            else:
                key = "flank"
            nadded += addSegment( key, last.end, this.start, (last,this), options )

    return nadded

##-----------------------------------------------------------------------------
def buildTerritories( iterator, fasta, options ):
    """build gene territories. 

    Exons in a gene are merged and the resulting 
    segments enlarged by --radius. Territories
    overlapping are divided in the midpoint between
    the two genes.
    """

    ninput, noutput, nambiguous = 0, 0, 0

    dr = 2 * options.radius

    prev_pos = 0
    last_contig = None
    gff = None

    def _iterator( iterator ):
        """yield gene plus the locations of the end of the previous gene and start of next gene"""

        last_end, prev_end = 0, 0
        last_contig = None
        last = None
        for matches in GFF.iterator_overlaps( iterator ):

            this_start = min( [ x.start for x in matches ] )
            this_end = max( [ x.end for x in matches ] )
            this_contig = matches[0].contig

            if last_contig != this_contig: 
                if last:
                    yield prev_end, last, fasta.getLength( last_contig )
                last_end, prev_end = 0, 0
            else:
                yield prev_end, last, this_start

            prev_end = last_end
            last_end = this_end
            last = matches
            last_contig = this_contig
            
        if last:
            yield prev_end, last, fasta.getLength( last_contig )

    for last_end, matches, next_start in _iterator( iterator ):

        gff = GTF.Entry().copy( matches[0] )

        start = min( [ x.start for x in matches ] )
        end = max( [ x.end for x in matches ] )

        d = start - last_end
        if d < dr:
            start -= d // 2
        else:
            start -= options.radius

        d = next_start - end
        if d < dr:
            end += d // 2
        else:
            end += options.radius

        gff.gene_id = ":".join( set( [x.gene_id for x in matches ] ) )
        gff.transcript_id = gff.gene_id
        gff.start, gff.end = start, end

        nsegments = len(matches)
        if nsegments > 1: 
            gff.addAttribute("ambiguous", nsegments )
            nambiguous += 1

        assert gff.start < gff.end, "invalid segment: %s" % str(gff)
        options.stdout.write( str(gff) + "\n" )
        noutput += 1
        
    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nambiguous=%i\n" % (ninput, noutput, nambiguous ))

##-----------------------------------------------------------------------------
def fullSegmentation( iterator, fasta, options ):
    """perform a full segmentation of the genome (UTR, exon, intron ...)
    """

    ninput, noutput, nadded, nambiguous, nframeshifts, nunknown = 0, 0, 0, 0, 0, 0
    last = None
    is_ambiguous = False

    for this in iterator:
        ninput += 1

        E.debug( "last=%s" % str(last))
        E.debug( "this=%s" % str(this))
        E.debug( "is_ambiguous=%s" % str(is_ambiguous))

        if last and last.contig == this.contig:
            # check if file is sorted correctly
            assert last.start <= this.start, "input file needs to be sorted by contig, start" 
            if last.end <= this.start:
                if not is_ambiguous:
                    if last.gene_id != this.gene_id:
                        nadded += addIntergenicSegment( last, this, fasta, options )
                    else:
                        d = this.start - last.end  
                        if d >= options.min_intron_length:
                            nadded += addSegment( "intronic", last.end, this.start, last, options )
                        elif d <= options.max_frameshift_length:
                            nframeshifts += addSegment( "frameshift", last.end, this.start, last, options )
                        else:
                            nunknown += addSegment( "unknown", last.end, this.start, last, options )
                else:
                    if last.feature == this.feature and last.gene_id == this.gene_id:
                        nambiguous += addSegment( last.feature, last.end, this.start, last, options )
                    else:
                        nambiguous += addSegment( "ambiguous", last.end, this.start, last, options )
                    is_ambiguous = False
                last = this
            elif last.end > this.start:
                if last.gene_id != this.gene_id:
                    # flag next region as ambiguous
                    is_ambiguous = True
                last.end = this.end
        else:
            nadded += addIntergenicSegment( last, None, fasta, options )
            nadded += addIntergenicSegment( None, this, fasta, options )
            last = this

        options.stdout.write( "%s\n" % str(this) )
        noutput += 1
            
    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nadded=%i, nambiguous=%i, nframeshifts=%i, nunknown=%i\n" % (ninput, noutput, nadded, nambiguous, nframeshifts, nunknown) )

def annotateExons( iterator, fasta, options ):
    """annotate exons within iterator."""

    gene_iterator = GTF.gene_iterator( iterator )

    ninput, noutput, noverlapping = 0, 0, 0

    for this in gene_iterator:
        ninput += 1
        intervals = collections.defaultdict( list )
        ntranscripts = len(this)

        is_negative_strand = Genomics.IsNegativeStrand( this[0][0].strand )

        for exons in this:
            # make sure these are sorted correctly
            exons.sort( key=lambda x: x.start )
            if is_negative_strand: exons.reverse()

            nexons = len(exons)
            for i, e in enumerate(exons):
                intervals[ (e.start, e.end) ].append( (i+1,nexons) )

        gtf = GTF.Entry()
        gtf.fromGFF( this[0][0], this[0][0].gene_id, this[0][0].gene_id )
        gtf.addAttribute( "ntranscripts", ntranscripts )        

        gtfs = []
        for r, pos in intervals.iteritems():

            g = GTF.Entry().copy(gtf)
            g.start, g.end = r
            g.addAttribute( "nused", len(pos) )
            g.addAttribute( "pos", ",".join( ["%i:%i" % x for x in pos ] ) )
            gtfs.append( g )
        
        gtfs.sort( key=lambda x: x.start )
        
        for g in gtfs:
            options.stdout.write( "%s\n" % str(g) )

        # check for exon overlap    
        intervals = [ (g.start, g.end) for g in gtfs ]
        nbefore = len(intervals)
        nafter = len( Intervals.combine( intervals ) )
        if nafter != nbefore: noverlapping += 1

        noutput += 1
            
    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, noverlapping=%i\n" % (ninput, noutput, noverlapping) )

def annotatePromoters( iterator, fasta, options ):
    """annotate promoters within iterator.

    Only protein_coding genes are annotated.
    """

    gene_iterator = GTF.gene_iterator( iterator )

    ngenes, ntranscripts, npromotors = 0, 0, 0

    for gene in gene_iterator:
        ngenes += 1
        is_negative_strand = Genomics.IsNegativeStrand( gene[0][0].strand )
        lcontig = fasta.getLength( gene[0][0].contig )
        promotors = []
        transcript_ids = []
        for transcript in gene:

            ntranscripts += 1
            mi, ma = min( [x.start for x in transcript ] ), max( [x.end for x in transcript ] )
            transcript_ids.append( transcript[0].transcript_id )
            # if tss is directly at start/end of contig, the tss will be within an exon.
            # otherwise, it is outside an exon.
            if is_negative_strand:
                promotors.append( (min( lcontig-options.promotor, ma), min(lcontig, ma + options.promotor)) )
            else:
                promotors.append( (max(0,mi - options.promotor), max(options.promotor,mi)) )

        if options.merge_promotors:
            # merge the promotors (and rename - as sort order might have changed)
            promotors = Intervals.combine( promotors )
            transcript_ids = ["%i" % (x+1) for x in range(len(promotors) )]
            
        gtf = GTF.Entry()
        gtf.fromGFF( gene[0][0], gene[0][0].gene_id, gene[0][0].gene_id )
        gtf.source = "promotor"

        x = 0
        for start, end in promotors:
            gtf.start, gtf.end = start, end
            gtf.transcript_id = transcript_ids[x]
            options.stdout.write( "%s\n" % str(gtf) )
            npromotors += 1
            x += 1

    E.info( "ngenes=%i, ntranscripts=%i, npromotors=%i" % (ngenes, ntranscripts, npromotors) )

def annotateRegulons( iterator, fasta, options ):
    """annotate regulons within iterator.

    Only protein_coding genes are annotated.
    """

    gene_iterator = GTF.gene_iterator( iterator )

    ngenes, ntranscripts, nregulons = 0, 0, 0

    for gene in gene_iterator:
        ngenes += 1
        is_negative_strand = Genomics.IsNegativeStrand( gene[0][0].strand )
        lcontig = fasta.getLength( gene[0][0].contig )
        regulons = []
        transcript_ids = []
        for transcript in gene:

            ntranscripts += 1
            mi, ma = min( [x.start for x in transcript ] ), max( [x.end for x in transcript ] )
            transcript_ids.append( transcript[0].transcript_id )
            # add range to both sides of tss
            regulons.append( (max( 0, mi - options.promotor), min( lcontig, ma + options.promotor) ) )

        if options.merge_promotors:
            # merge the regulons (and rename - as sort order might have changed)
            regulons = Intervals.combine( regulons )
            transcript_ids = ["%i" % (x+1) for x in range(len(regulons) )]
            
        gtf = GTF.Entry()
        gtf.fromGFF( gene[0][0], gene[0][0].gene_id, gene[0][0].gene_id )
        gtf.source = "promotor"

        x = 0
        for start, end in regulons:
            gtf.start, gtf.end = start, end
            gtf.transcript_id = transcript_ids[x]
            options.stdout.write( "%s\n" % str(gtf) )
            nregulons += 1
            x += 1

    E.info( "ngenes=%i, ntranscripts=%i, nregulons=%i" % (ngenes, ntranscripts, nregulons) )


def annotateTTS( iterator, fasta, options ):
    """annotate termination sites within iterator.

    Only protein_coding genes are annotated.
    """

    gene_iterator = GTF.gene_iterator( iterator )

    ngenes, ntranscripts, npromotors = 0, 0, 0

    for gene in gene_iterator:
        ngenes += 1
        is_negative_strand = Genomics.IsNegativeStrand( gene[0][0].strand )
        lcontig = fasta.getLength( gene[0][0].contig )
        tts = []
        transcript_ids = []
        for transcript in gene:

            ntranscripts += 1
            mi, ma = min( [x.start for x in transcript ] ), max( [x.end for x in transcript ] )
            transcript_ids.append( transcript[0].transcript_id )
            # if tts is directly at start/end of contig, the tss will be within an exon.
            # otherwise, it is outside an exon.
            if is_negative_strand:
                #tss.append( (min( lcontig-options.promotor, ma), min(lcontig, ma + options.promotor)) )
                tts.append( (max(0, mi-options.promotor), max(options.promotor,mi)) )
            else:
                #tss.append( (max(0,mi - options.promotor), max(options.promotor,mi)) )
                tts.append( (min(ma,lcontig-options.promotor), min(lcontig,ma + options.promotor)) )

        if options.merge_promotors:
            # merge the promotors (and rename - as sort order might have changed)
            tts = Intervals.combine( tts )
            transcript_ids = ["%i" % (x+1) for x in range(len(tts) )]
            
        gtf = GTF.Entry()
        gtf.fromGFF( gene[0][0], gene[0][0].gene_id, gene[0][0].gene_id )
        gtf.source = "tts"

        x = 0
        for start, end in tts:
            gtf.start, gtf.end = start, end
            gtf.transcript_id = transcript_ids[x]
            options.stdout.write( "%s\n" % str(gtf) )
            npromotors += 1
            x += 1

    if options.loglevel >= 1:
        options.stdlog.write( "# ngenes=%i, ntranscripts=%i, ntss=%i\n" % (ngenes, ntranscripts, npromotors) )

##------------------------------------------------------------
if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id: gtf2gff.py 2861 2010-02-23 17:36:32Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default]."  )

    parser.add_option("-i", "--ignore-missing", dest="ignore_missing", action="store_true",
                      help="Ignore transcripts on contigs that are not in the genome-file [default=%default]."  )

    parser.add_option("-s", "--restrict-source", dest="restrict_source", type="choice",
                      choices=("protein_coding", "pseudogene" ),
                      help="restrict input by source [default=%default]."  )

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=("full", "genome", "territories", "exons", "promotors", "tts", "regulons"),
                      help="method for defining segments [default=%default]."  )

    parser.add_option("-r", "--radius", dest="radius", type="int",
                      help="radius of a territory [default=%default]."  )

    parser.add_option("-f", "--flank", dest="flank", type="int",
                      help="size of the flanking region next to a gene [default=%default]."  )

    parser.add_option("-p", "--promotor", dest="promotor", type="int",
                      help="size of a promotor region [default=%default]."  )

    parser.add_option( "--merge-promotors", dest="merge_promotors", action="store_true",
                       help="merge promotors [default=%default]."  )

    parser.add_option( "--min-intron-length", dest="min_intron_length", type="int",
                      help="minimum intron length. If the distance between two consecutive exons is smaller, the region will be marked 'unknown' [default=%default]."  )

    parser.set_defaults(
        genome_file = None,
        flank = 1000,
        max_frameshift_length = 4,
        min_intron_length = 30,
        ignore_missing = False,
        restrict_source = None,
        method = "genome",
        radius = 50000,
        promotor = 5000,
        merge_promotors = False,
        )

    (options, args) = E.Start( parser )
    
    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )
    else:
        fasta = None

    if options.restrict_source:
        iterator = GTF.iterator_filtered( GTF.iterator(options.stdin), 
                                          source = options.restrict_source )
    elif options.method == "promotors" or options.method == "tts":
        iterator = GTF.iterator_filtered( GTF.iterator(options.stdin), source = "protein_coding" )
    else:
        iterator = GTF.iterator(options.stdin)

    if options.method == "full" or options.method == "genome":
        segmentor = fullSegmentation( iterator, fasta, options )
    elif options.method == "territories":
        segmentor = buildTerritories( iterator, fasta, options )
    elif options.method == "exons":
        segmentor = annotateExons( iterator, fasta, options )
    elif options.method == "promotors":
        segmentor = annotatePromoters( iterator, fasta, options )
    elif options.method == "regulons":
        segmentor = annotateRegulons( iterator, fasta, options )
    elif options.method == "tts":
        segmentor = annotateTTS( iterator, fasta, options )

    E.Stop()

