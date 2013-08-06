################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: snp2maf.py 2875 2010-03-27 17:42:04Z andreas $
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
snp2map.py - build maf formatted multiple genomic alignments from SNPs
======================================================================

:Author: Andreas Heger
:Release: $Id: snp2maf.py 2875 2010-03-27 17:42:04Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

construct alleles of a reference sequences based on
variants in an sqlite database.

Caveats
-------

* no phasing
* the two allelic sequences should be built per genome

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
import gzip
import sqlite3

import numpy
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Genomics as Genomics
import CGAT.GTF as GTF
import CGAT.Variants as Variants
import alignlib

def alignIndels( all_alleles, colcounts, extend_by = 0 ):
    '''align all indel-regions.'''

    aa = alignlib.makeAlignatorDPFull( alignlib.ALIGNMENT_LOCAL, 0, 0 )     
    alignator = alignlib.makeMultipleAlignatorSimple( aa)

    ids = all_alleles.keys()

    for x,c in enumerate(colcounts):
        if c <= 1: continue
        sequences = alignlib.StringVector()
        for sid in ids:
            for allele in all_alleles[sid]:
                sequences.append( allele[x] )

        mali = alignlib.makeMultAlignment()
        alignator.align( mali, sequences )
        realigned = []
        for line in str(alignlib.MultAlignmentFormatPlain( mali, sequences )).split("\n")[:-1]:
            data = line[:-1].split("\t")
            realigned.append( data[1] )
        assert len(realigned) == len(sequences)

        l = max( [len(r) for r in realigned] )
        i = 0
        for sid in ids:
            for allele in all_alleles[sid]:
                if realigned[i]: allele[x] = realigned[i] 
                else: allele[x] = "-" * l 
                i += 1
                
        colcounts[x] = l


def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: snp2maf.py 2875 2010-03-27 17:42:04Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default]."  )
    parser.add_option("-t", "--tracks", dest="tracks", type="string", action="append",
                      help="tracks (tablenames) to use in sqlite database [default=%default]."  )
    parser.add_option("-d", "--database", dest="database", type="string",
                      help="sqlite3 database [default=%default]."  )
    parser.add_option("-r", "--reference", dest="reference", type="string",
                      help="name of reference [default=%default]."  )
    parser.add_option("-i", "--is-gtf", dest="is_gtf", action="store_true",
                      help="if set, the gene_id will be added to the alignment header [default=%default]."  )
    parser.add_option("-z", "--compress", dest="compress", action="store_true",
                      help="compress output with gzip [default=%default]."  )
    parser.add_option("-p", "--pattern", dest="pattern_track", type="string",
                      help="regular expression pattern for track [default=%default].")

    parser.set_defaults(
            genome_file = None,
            tracks = [],
            database = "csvdb",
            output = [],
            border = 0,
            reference_name = "reference",
            pattern_track = "(\S+)",
            is_gtf = True,
            compress = False,
            )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

    ninput, nskipped, noutput = 0, 0, 0

    if not options.database or not options.tracks:
        raise ValueError("please supply both database and tracks")

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )
    else:
        fasta = None

    if options.is_gtf:
        infile_gff = GTF.iterator( options.stdin)
    else:
        infile_gff = GTF.iterator( options.stdin)

    dbhandle = sqlite3.connect( options.database )

    statement = '''SELECT pos, reference, genotype 
                   FROM %(track)s
                   WHERE contig = '%(contig)s' AND 
                   pos BETWEEN %(extended_start)s and %(extended_end)s
                '''

    counts = E.Counter()
    tracks = options.tracks
    try:
        translated_tracks = [ re.search( options.pattern_track, track).groups()[0] for track in tracks ]
    except AttributeError:
        raise AttributeError( "pattern `%s` does not match input tracks." % options.pattern_track )

    if options.compress:
        outfile = gzip.GzipFile( fileobj = options.stdout )
    else:
        outfile = options.stdout

    outfile.flush()
    outfile.write( "##maf version=1 program=snp2maf.py\n\n") 

    for gff in infile_gff:
        counts.input += 1
        
        contig = gff.contig
        strand = gff.strand
        lcontig = fasta.getLength( contig )
        region_start, region_end = gff.start, gff.end
        if contig.startswith("chr"): contig = contig[3:]
        extended_start = region_start - options.border
        extended_end = region_end + options.border
        is_positive = Genomics.IsPositiveStrand( strand )

        E.info("processing %s" % str(gff))

        # collect all variants
        all_variants = []
        for track in options.tracks:
            cc = dbhandle.cursor()
            cc.execute( statement % locals() )
            all_variants.append( map(Variants.Variant._make, cc.fetchall()) )
            cc.close()
            
        E.debug("%s:%i..%i collected %i variants for %i tracks" % (contig, 
                                                                   region_start, region_end, 
                                                                   sum( [len(x) for x in all_variants] ),
                                                                   len(all_variants)))

        reference_seq = fasta.getSequence( contig, "+", region_start, region_end ) 
        lseq = len(reference_seq)
        alleles = collections.defaultdict( list )
            
        # build allele sequences for track and count maximum chars per mali column
        colcounts = numpy.ones( lseq )
        for track, variants in zip(translated_tracks, all_variants):
            variants = Variants.updateVariants( variants, lcontig, "+" )
            a = Variants.buildAlleles(reference_seq, 
                                      variants, 
                                      reference_start = region_start )
            
            alleles[track] = a
            for allele in a:
                for pos, c in enumerate( allele):
                    colcounts[pos] = max(colcounts[pos], len(c))

        # realign gapped regions
        alignIndels( alleles, colcounts ) 
        
        if options.is_gtf:
            outfile.write( "a gene_id=%s\n" % gff.gene_id )
        else:
            outfile.write( "a\n" )

        maf_format = "s %(name)-30s %(pos)9i %(size)6i %(strand)s %(lcontig)9i %(seq)s\n" 
        def __addGaps( sequence, colcounts ):
            '''output gapped sequence.'''
            r = []
            for x,c in enumerate(sequence):
                r.append( c + "-" * (colcounts[x]-len(c)))
            return "".join(r)

        name = ".".join((options.reference, contig))
        if is_positive:
            pos = region_start
        else:
            pos = lcontig - region_start

        size = lseq
        seq = __addGaps( reference_seq, colcounts )
        outfile.write( maf_format % (locals()) )
        
        for track in translated_tracks:
            for aid, allele in enumerate(alleles[track]):
                seq = __addGaps( allele, colcounts )
                if not is_positive: Genomics.complement( seq )
                size = len(seq) - seq.count("-")
                name = ".".join( (track + "-%i" % aid, contig))
                outfile.write( maf_format % (locals()) )
                
        outfile.write("\n")

    E.info( "%s" % str(counts))

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

