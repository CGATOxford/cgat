'''
gff2gff.py - manipulate gff files
=================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics Intervals GFF Manipulation

Purpose
-------

This scripts reads a :term:`gff` formatted file, operates
on it and outputs the new intervals in :term:`gff` format.

Extension options:

``--extend``
   extend existing features

``--add-up-flank/--add-down-flank``
   add an upstream/downstearm flanking segment to first/last exon of a group.

Segment transformations:

``--crop``
   crop features according to features in a separate gff file.

``--crop-unique``
   remove non-unique features from gff file.

``--merge-features``
   merge consecutive features.

``--join-features``
   group consecutive features

Filtering options:

``--filter-range``
   extract features overlapping a chromosomal range.

``--sanitize``
   reconcile chromosome names between ENSEMBL/UCSC or with an indexed
   genomic fasta file (see :doc:`index_fasta`). Raises an exception if
   an unknown contig is found (unless ``--skip-missing`` is set).

``--remove-contigs``
   remove contig matching to a series of regular expressions.



Usage
-----

Make sure that a :term:`gff` formatted file contains only features on placed chromosomes::

   cat in.gff 
   | gff2gff.py --sanitize=genome --genome-file=hg19 --skip-missing 
   | gff2gff.py --remove-contigs="chrUn,_random" > gff.out

Type::

   python gff2gff.py --help

for command line help.

Command line options
--------------------

'''

import sys
import string
import re
import optparse


import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import CGAT.AGP as AGP
import CGAT.Genomics as Genomics
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Intervals as Intervals
import collections
import numpy
import bx.intervals.io
import bx.intervals.intersection

def combineGFF( gffs, options, merge = True ):
    """join intervals in gff file.

    Note: strandedness is ignored
    """

    if merge:
        min_distance, max_distance, min_features, max_features = map(int, options.merge_features.split(",") )
    else:
        min_distance, max_distance, min_features, max_features = map(int, options.join_features.split(",") )

    E.info( "joining features: min distance=%i, max_distance=%i, at least %i and at most %i features." %
            (min_distance, max_distance, min_features, max_features) )

    def iterate_chunks( gffs ):

        last = gffs.next()
        to_join = [last]

        for gff in gffs:
            d = gff.start - last.end
            if gff.contig == last.contig:
                assert gff.start >= last.start, "input file should be sorted by contig and position: d=%i:\n%s\n%s\n" % (d, last, gff)
            
            if gff.contig != last.contig or \
                    (max_distance and d > max_distance) or \
                    (min_distance and d < min_distance) or \
                    (max_features and len(to_join) >= max_features):
                
                if min_features or len(to_join) >= min_features: 
                    yield to_join
                to_join = []

            last = gff
            to_join.append( gff )
            
        if len(to_join) >= min_features: yield to_join
        raise StopIteration

    id = 1
    ninput, noutput, nfeatures= 0, 0, 0

    if merge:
        for to_join in iterate_chunks(gffs):

            ninput += 1
            y = GTF.Entry()
            t = options.output_format % id
            y.fromGTF( to_join[0], t, t )
            y.start = to_join[0].start
            y.end = to_join[-1].end
            
            options.stdout.write( "%s\n" % str(y ) )
            nfeatures += 1

            noutput += 1
            id += 1
    else:

        for to_join in iterate_chunks(gffs):

            ninput += 1
            for x in to_join:
                y = GTF.Entry()
                t = options.output_format % id
                y.fromGTF( x, t, t )
                options.stdout.write( "%s\n" % str(y ) )
                nfeatures += 1

            noutput += 1
            id += 1

    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nfeatures=%i\n" % (ninput, noutput, nfeatures) )
        
def cropGFFUnique( gffs, options ):
    """crop intervals in gff file.

    only unique regions are kept. This method ignores the feature field.

    If options.ignore_strand is set, strand information for cropping
    is ignored.
    """

    # read regions to crop with and convert intervals to intersectors
    E.info( "reading gff for cropping: started." )
    gffs = list(gffs)
    if len(gffs) == 0: return

    outf = options.stdout

    def _cmp_without_strand( this, last):
        return this.contig != last.contig

    def _cmp_with_strand( this, last):
        return this.contig != last.contig or this.strand != last.strand

    if options.ignore_strand:
        gffs.sort( key = lambda x: (x.contig, x.start ) )
        comp = _cmp_without_strand
    else:
        gffs.sort( key = lambda x: (x.contig, x.strand, x.start ) )
        comp = _cmp_with_strand

    last = gffs[0]
    c = E.Counter()

    c.input = len(gffs)

    for this in gffs[1:]:
        if comp( this, last) or last.end <= this.start:
            # no overlap
            if last.start < last.end: 
                c.output += 1
                outf.write( "%s\n" % str(last) )
            last = this
            
        elif this.end <= last.start:
            # this ends before last (happens in multiway overlaps)
            # nothing happens
            pass
            
        elif last.end <= this.end:
            # last ends before this
            l = last.end
            last.end = this.start
            this.start = l
            if last.start < last.end: 
                E.info( "overlap" )
                E.info( str(this) )
                E.info( str(last) )
                c.overlaps += 1
                c.output += 1
                outf.write( "%s\n" % str(last) )
            last = this

        elif last.end > this.end:
            # this within last - split last
            l = last.end
            last.end = this.start
            if last.start < last.end: 
                c.output += 1
                c.splits += 1
                outf.write( "%s\n" % str(last) )
            last.start = this.end
            last.end = l

    if last.start < last.end: 
        c.output += 1
        outf.write( "%s\n" % str(last) )

    E.info( "cropping finished: %s" % str(c) )

def cropGFF( gffs, options ):
    """crop intervals in gff file."""

    # read regions to crop with and convert intervals to intersectors
    E.info( "reading gff for cropping: started." )

    other_gffs = GTF.iterator( IOTools.openFile( options.crop, "r") )
    cropper = GTF.readAsIntervals( other_gffs )
    ntotal = 0
    for contig in cropper.keys():
        intersector = bx.intervals.intersection.Intersecter()
        for start, end in cropper[contig]:
            intersector.add_interval( bx.intervals.Interval(start,end) )
            ntotal += 1
        cropper[contig] = intersector

    E.info( "reading gff for cropping: finished." )
    E.info( "reading gff for cropping: %i contigs with %i intervals." % (len(cropper), ntotal ) )

    ninput, noutput, ncropped, ndeleted = 0, 0, 0, 0

    # do the actual cropping
    for gff in gffs:
        
        ninput += 1

        if gff.contig in cropper:
            start, end = gff.start, gff.end
            overlaps = cropper[gff.contig].find(start, end)

            if overlaps:
                l = end - start
                a = numpy.ones( l )
                for i in overlaps:
                    s = max( 0, i.start - start)
                    e = min( l, i.end - start)
                    a[s:e] = 0
                    
                segments = Intervals.fromArray( a )
                
                if len(segments) == 0:
                    ndeleted += 1
                else:
                    ncropped += 1

                for s, e in segments:
                    gff.start, gff.end = s + start, e + start
                    noutput += 1
                    options.stdout.write( "%s\n" % gff )
                    
                continue

        noutput += 1
        options.stdout.write( "%s\n" % gff )
            
    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, ncropped=%i, ndeleted=%i\n" % (ninput, noutput, ncropped, ndeleted) )

##------------------------------------------------------------------------
def main( argv = None ):
    
    if argv == None: argv = sys.argv

    parser = E.OptionParser( version = "%prog version: $Id: gff2gff.py 2868 2010-03-03 10:19:52Z andreas $")

    parser.add_option( "-f", "--forward-coordinates", dest="forward_coordinates", 
                      help="translate to forward coordinates.", action="store_true"  )

    parser.add_option( "--forward-strand", dest="forward_strand", 
                      help="convert to forward strand.", action="store_true"  )

    parser.add_option( "--ignore-strand", dest="ignore_strand", 
                      help="ignore strand information.", action="store_true"  )

    parser.add_option( "--is-gtf", dest="is_gtf", action="store_true",
                       help="input will be treated as gtf [default=%default]." ) 

    parser.add_option( "--add-up-flank", dest="add_up_flank", type="int", 
                      help="add an upstream flanking segment to first exon of a group."  )

    parser.add_option( "--add-down-flank", dest="add_down_flank", type="int", 
                      help="add a downstream flanking segment to last segment of a group."  )

    parser.add_option( "--extend", dest="extend", 
                      help="extend the existing features.", action="store_true"  )

    parser.add_option( "-c", "--contigs", dest="input_filename_contigs", type="string",
                       help="filename with contig lenghts."  )

    parser.add_option( "--filename-agp", dest="input_filename_agp", type="string",
                       help="agp file to map coordinates from contigs to scaffolds." )
    
    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome."  )

    parser.add_option( "--complement-groups", dest="complement_groups", action="store_true",
                       help="""complement groups. Will write introns from exons [%default].""" )

    parser.add_option( "--group-field", dest="group_field", type="string",
                       help="""gff field/attribute to group by such as gene_id, transrcipt_id, ... [%default].""" )

    parser.add_option( "--combine-groups", dest="combine_groups", action="store_true",
                       help="""combine groups.""" )

    parser.add_option( "--filter-range", dest="filter_range", type="string",
                       help="""extract all elements overlapping a range. A range is specified by eithor 'contig:from..to', 'contig:+:from..to', or 'from,to' .""" )

    parser.add_option( "--join-features", dest="join_features", type="string",
                       help="join features into a single transcript. Consecutive features are grouped "
                       " into the same transcript/gene. This metdo expects a string of for numbers ``a,b,c,d`` " 
                       " as input with:"
                       " a,b=minimum/maximum distance between features, "
                       " c,d=minimum,maximum number of features.""" )

    parser.add_option( "--merge-features", dest="merge_features", type="string",
                       help="merge features. Consecutive features are merged into a single feature. "
                       "This method expects a string of four numbers ``a,b,c,d`` as input; "
                       "a,b=minimum/maximum distance between features, "
                       "c,d=minimum,maximum number of features." )

    parser.add_option( "--crop-unique", dest="crop_unique", action="store_true",
                       help = "crop overlapping intervals, keeping only intervals that are unique [default=%default]" )

    parser.add_option( "--crop", dest="crop", type="string",
                       help="""crop features in gff file with features in another file. If a feature falls in the middle of another, two entries will be output.""" )

    parser.add_option( "--sanitize", dest="sanitize", type="choice",
                       choices=( "ucsc", "ensembl", "genome" ),
                       help="sanitize chr names for ucsc or ensembl or use the genome translator [%default]." )

    parser.add_option( "--skip-missing", dest="skip_missing", action="store_true",
                       help="skip entries on missing contigs. Otherwise an exception is raised [%default]." )

    parser.add_option( "--remove-contigs", dest="remove_contigs", type="string", action="store",
                       help="a comma separated list of regular expressions specifying contigs to be removed when runnnig sanitize [%default]." )

    parser.set_defaults(
        forward_coordinates = False,
        forward_strand = False,
        input_filename_contigs = False,
        input_filename_agp = False,
        genome_file = None,
        sanitize = None,
        add_up_flank = None,
        add_down_flank = None,
        extend = False,
        complement_groups=False,
        combine_groups=False,
        crop = None,
        crop_unique = False,
        ignore_strand = False,
        filter_range = None,
        join_features = None,
        merge_features = None,
        output_format = "%06i",
        skip_missing = False,
        remove_contigs = None,
        is_gtf = False,
        group_field = None,
        )

    (options, args) = E.Start( parser, argv = argv )

    if options.input_filename_contigs:
        contigs = Genomics.ReadContigSizes( IOTools.openFile( options.input_filename_contigs, "r") )

    if options.genome_file:
        genome_fasta = IndexedFasta.IndexedFasta( options.genome_file )
        contigs = genome_fasta.getContigSizes()
    else:
        genome_fasta = None

    if (options.forward_coordinates or options.forward_strand) and not contigs: 
        raise ValueError( "inverting coordinates requires genome file")

    if options.input_filename_agp:
        agp = AGP.AGP()
        agp.readFromFile( IOTools.openFile( options.input_filename_agp, "r" ) )
    else:
        agp = None
        
    gffs = GTF.iterator( options.stdin )

    if options.add_up_flank or options.add_down_flank:
        
        if options.is_gtf:
            iterator = GTF.flat_gene_iterator( gffs )
        else:
            iterator = GTF.joined_iterator( gffs, options.group_field )
        
        for chunk in iterator:
            is_positive = Genomics.IsPositiveStrand( chunk[0].strand )
            chunk.sort( lambda x, y: cmp( x.start, y.start) )
            lcontig = contigs[chunk[0].contig]

            if options.extend:
                if options.add_up_flank:
                    if is_positive:
                        chunk[0].start = max( 0, chunk[0].start - options.add_up_flank )
                    else:
                        chunk[-1].end = min( lcontig, chunk[-1].end + options.add_up_flank )
                if options.add_down_flank:
                    if is_positive:
                        chunk[-1].end = min( lcontig, chunk[-1].end + options.add_down_flank )
                    else:
                        chunk[0].start = max( 0, chunk[0].start - options.add_down_flank )
            else:
                if options.add_up_flank:
                    gff = GTF.Entry()
                    if is_positive:
                        gff.copy( chunk[0] )
                        gff.end = gff.start
                        gff.start = max( 0, gff.start - options.add_up_flank )
                        chunk.insert(0, gff)
                    else:
                        gff.copy( chunk[-1] )
                        gff.start = gff.end
                        gff.end = min( lcontig, gff.end + options.add_up_flank )
                        chunk.append( gff )
                    gff.feature = "5-Flank"
                    gff.mMethod = "gff2gff"
                if options.add_down_flank:
                    gff = GTF.Entry()
                    if is_positive:
                        gff.copy( chunk[-1] )
                        gff.start = gff.end
                        gff.end = min( lcontig, gff.end + options.add_up_flank )
                        chunk.append( gff )
                    else:
                        gff.copy( chunk[0] )
                        gff.end = gff.start
                        gff.start = max( 0, gff.start - options.add_up_flank )
                        chunk.insert(0, gff)
                    gff.feature = "3-Flank"
                    gff.mMethod = "gff2gff"

            if not is_positive: chunk.reverse()

            for gff in chunk:
                options.stdout.write( str(gff) + "\n" )

    elif options.complement_groups:
        
        iterator = GTF.joined_iterator( gffs, group_field = options.group_field )

        for chunk in iterator:
            if options.is_gtf:
                chunk = [ x for x in chunk if x.feature == "exon" ]
                if len(chunk) == 0: continue
            chunk.sort()
            x = GTF.Entry()
            x.copy( chunk[0] )
            x.start = x.end
            x.feature = "intron"
            for c in chunk[1:]:
                x.end = c.start
                options.stdout.write( str(x) + "\n" )
                x.start = c.end
            
    elif options.combine_groups:
        
        iterator = GTF.joined_iterator( gffs )

        for chunk in iterator:
            chunk.sort()
            x = GTF.Entry()
            x.copy( chunk[0] )
            x.end = chunk[-1].end
            x.feature = "segment"
            options.stdout.write( str(x) + "\n" )

    elif options.join_features:

        combineGFF( gffs, options, merge = False )

    elif options.merge_features:

        combineGFF( gffs, options, merge = True )

    elif options.crop:
        
        cropGFF( gffs, options )

    elif options.crop_unique:
        
        cropGFFUnique( gffs, options )

    elif options.filter_range:

        contig, strand, interval = None, None, None
        try:
            contig, strand, start, sep, end = re.match( "(\S+):(\S+):(\d+)(\.\.|-)(\d+)", options.filter_range ).groups()
        except AttributeError:
            pass
        
        if not contig:
            try:
                contig, start, sep, end = re.match( "(\S+):(\d+)(\.\.|-)(\d+)", options.filter_range ).groups()
                strand = None
            except AttributeError:
                pass
            
        if not contig:
            try:
                start,end = re.match( "(\d+)(\.\.|\,|\-)(\d+)", options.filter_range).groups()
            except AttributeError:
                raise "can not parse range %s" % options.filter_range
            contig = None
            strand = None
            
        if start:
            interval = (int(start), int(end) )
        else:
            interval = None

        if options.loglevel >= 2:
            options.stdlog.write("# filter: contig=%s, strand=%s, interval=%s\n" % (str(contig), str(strand), str(interval)))
            options.stdlog.flush()

        for gff in GTF.iterator_filtered( gffs, contig = contig, strand = strand, interval = interval ):
            options.stdout.write( str(gff) + "\n" )

    elif options.sanitize:

        def toUCSC(id):
            if not id.startswith("contig") and not id.startswith("chr"):
                id = "chr%s" % id
            return id

        def toEnsembl(id):
            if id.startswith("contig"): return id[len("contig"):]
            if id.startswith("chr"): return id[len("chr"):]
            return id

        if options.sanitize == "genome":
            if genome_fasta == None:
                raise ValueError("please specify --genome-file= when using --sanitize=genome")
            f = genome_fasta.getToken
        elif options.sanitize == "ucsc":
            f = toUCSC
        elif options.sanitize == "ensembl":
            f = toEnsembl

        skipped_contigs = collections.defaultdict( int )
        outofrange_contigs = collections.defaultdict( int )
        filtered_contigs = collections.defaultdict( int )
        
        for gff in gffs:
            try:
                gff.contig = f( gff.contig )
            except KeyError, msg:
                if options.skip_missing:
                    skipped_contigs[gff.contig] += 1
                    continue
                else:
                    raise

            if genome_fasta:
                lcontig = genome_fasta.getLength( gff.contig )
                if lcontig < gff.end:
                    outofrange_contigs[gff.contig] += 1
                    continue

            if options.remove_contigs:
                to_remove = [re.compile(x) for x in options.remove_contigs.split(",")]
                if any([x.match(gff.contig) for x in to_remove]):
                    filtered_contigs[gff.contig] += 1 
                    continue
                
            options.stdout.write( str(gff) + "\n" )

        if skipped_contigs:
            E.info( "skipped %i entries on %i contigs: %s" % (sum(skipped_contigs.values()),
                                                              len(skipped_contigs.keys()),
                                                              str(skipped_contigs) ) )
        if outofrange_contigs:
            E.warn( "skipped %i entries on %i contigs because they are out of range: %s" %\
                    (sum(outofrange_contigs.values()),
                     len(outofrange_contigs.keys()),
                     str(outofrange_contigs) ) )

        if filtered_contigs:
            E.info( "filtered out %i entries on %i contigs: %s" %\
                        (sum(filtered_contigs.values()), 
                         len(filtered_contigs.keys()), 
                         str(filtered_contigs) ) ) 

    else:

        
        for gff in gffs:

            if options.forward_coordinates:
                gff.invert( contigs[gff.contig] )

            if options.forward_strand:
                gff.invert( contigs[gff.contig] )
                gff.strand = "+"

            if agp:
                # note: this works only with forward coordinates
                gff.contig, gff.start, gff.end = agp.mapLocation( gff.contig, gff.start, gff.end )
        
            options.stdout.write( str(gff) + "\n" )

    E.Stop()
    
if __name__ == "__main__":
    sys.exit( main( sys.argv ))
