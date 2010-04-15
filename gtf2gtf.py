################################################################################
#
#   Gene prediction pipeline 
#
#   $Id: gtf2gtf.py 2861 2010-02-23 17:36:32Z andreas $
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
import os, sys, string, re, optparse, types, random, collections

USAGE="""python %s [OPTIONS] < in.gtf

convert a gtf merging all overlapping exons per gene.

This script expects the gtf file to be sorted by genes by contig and then by 
position.
"""

import GTF
import Experiment as E
import IndexedFasta
import Genomics
import Intervals
import IOTools
import Components

##------------------------------------------------------------
if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id: gtf2gtf.py 2861 2010-02-23 17:36:32Z andreas $", usage = USAGE)

    parser.add_option("-m", "--merge-exons", dest="merge_exons", action="store_true",
                      help="merge overlapping exons of all transcripts within a gene [default=%default]."  )

    parser.add_option( "--merge-genes", dest="merge_genes", action="store_true",
                      help="merge overlapping genes if their exons overlap. This ignores the strand [default=%default]."  )

    parser.add_option("-u", "--with-utr", dest="with_utr", action="store_true",
                      help="include utr in merged transcripts [default=%default]."  )

    parser.add_option("-t", "--merge-transcripts", dest="merge_transcripts", action="store_true",
                      help="merge all transcripts within a gene. The entry will span the whole gene (exons and introns). The transcript does not include the UTR unless --with-utr is set. [default=%default]."  )

    parser.add_option("-i", "--merge-introns", dest="merge_introns", action="store_true",
                      help="merge all introns within a gene [default=%default]."  )

    parser.add_option("-g", "--set-transcript-to-gene", "--set-transcript2gene", dest="set_transcript2gene", action="store_true",
                      help="set the transcript_id to the gene_id [default=%default]."  )

    parser.add_option("-G", "--set-gene-to-transcript", "--set-gene2transcript", dest="set_gene2transcript", action="store_true",
                      help="set the gene_id to the transcript_id [default=%default]."  )

    parser.add_option("-d", "--set-score2distance", dest="set_score2distance", action="store_true",
                      help="set the score field to the distance to transcription start [default=%default]."  )

    parser.add_option( "--exons2introns", dest="exons2introns", action="store_true",
                       help="convert exons to introns for a gene. Introns are labelled with the feature 'intron'." )

    parser.add_option("-f", "--filter", dest="filter", type="choice",
                      choices=("gene", "transcript","longest-gene"),
                      help="filter by gene-id or transcript-id. If `longest-gene` is chosen, the longest gene is chosen case of overlapping genes [default=%default]."  )

    parser.add_option("-r", "--rename", dest="rename", type="choice",
                      choices=("gene", "transcript","longest-gene"),
                      help="rename genes or transcripts with a map given in apply. Those that can not be renamed are removed [default=%default]."  )

    parser.add_option("-a", "--apply", dest="filename_filter", type="string",
                      help="filename to filter with [default=%default]."  )

    parser.add_option( "--invert-filter", dest="invert_filter", action="store_true",
                      help="invert filter (like grep -v) [default=%default]."  )

    parser.add_option("--sample-size", dest="sample_size", type="int",
                      help="extract a random sample of size # if the option --filter is set[default=%default]." )

    parser.add_option("--intron-min-length", dest="intron_min_length", type="int",
                      help="minimum length for introns [default=%default]."  )

    parser.add_option("--min-exons-length", dest="min_exons_length", type="int",
                      help="minimum length for gene (sum of exons) [default=%default]."  )

    parser.add_option("--intron-border", dest="intron_border", type="int",
                      help="number of residues to exclude at intron at either end [default=%default]."  )

    parser.add_option( "--transcripts2genes", dest="transcripts2genes", action="store_true",
                       help="cluster overlapping transcripts into genes." )

    parser.add_option( "--reset-strand", dest="reset_strand", action="store_true",
                       help="remove strandedness of features." )

    parser.add_option( "--remove-duplicates", dest="remove_duplicates", type="choice",
                       choices=("gene", "transcript"),
                       help="remove duplicates by gene/transcript [%default]" )



    parser.set_defaults(
        merge_exons = False,
        merge_transcripts = False,
        set_score2distance = False,
        set_gene2transcript = False,
        set_transcript2gene = False,
        filename_filter = None,
        filter = None,
        exons2introns = None,
        merge_genes = False,
        intron_border = None,
        intron_min_length = None,
        sample_size = 0,
        min_exons_length = 0,
        transripts2genes = False,
        reset_strand = False,
        with_utr = False,
        invert_filter = False,
        remove_duplications = None,
        )

    (options, args) = E.Start( parser )
    
    ninput, noutput, nfeatures, ndiscarded = 0, 0, 0, 0

    if options.set_transcript2gene:

        for gff in GTF.iterator(options.stdin):

            ninput += 1

            gff.setAttribute( "transcript_id", gff.gene_id)
            options.stdout.write( "%s\n" % str(gff) )                

            noutput += 1
            nfeatures += 1

    elif options.remove_duplicates:
        counts = collections.defaultdict(int)
        
        if options.remove_duplicates == "gene":
            gffs = GTF.gene_iterator(GTF.iterator(options.stdin), strict = False )
            f = lambda x: x[0][0].gene_id
            outf = lambda x: "\n".join( [ "\n".join( [ str(y) for y in xx] ) for xx in x] )
        elif options.remove_duplicates == "transcript":
            gffs = GTF.transcript_iterator(GTF.iterator(options.stdin), strict = False)
            f = lambda x: x[0].transcript_id
            outf = lambda x: "\n".join( [ str(y) for y in x] )

        store = []

        for entry in gffs:
            ninput += 1
            store.append(entry)
            id = f(entry)
            counts[id] += 1

        for entry in store:
            id = f(entry)
            if counts[id] == 1:
                options.stdout.write( outf(entry) + "\n" )
                noutput += 1
            else:
                ndiscarded += 1
                E.info("discarded duplicates for %s: %i" % (id, counts[id]))
            
    elif options.set_gene2transcript:

        for gff in GTF.iterator(options.stdin):

            ninput += 1

            gff.gene_id = gff.transcript_id
            options.stdout.write( "%s\n" % str(gff) )                
            noutput += 1
            nfeatures += 1

    elif options.merge_genes:

        gffs = GTF.iterator_sorted_chunks( GTF.flat_gene_iterator(GTF.iterator(options.stdin)) )
        
        def iterate_chunks( gff_chunks ):

            last = gff_chunks.next()
            to_join = [last]

            for gffs in gff_chunks:
                d = gffs[0].start - last[-1].end
                if gffs[0].contig == last[0].contig:
                    assert gffs[0].start >= last[0].start, "input file should be sorted by contig and position: d=%i:\n%s\n%s\n" % (d, str(last), str(gffs))

                if gffs[0].contig != last[0].contig or d > 0:
                    yield to_join
                    to_join = []

                last = gffs
                to_join.append( gffs )

            yield to_join
            raise StopIteration
        
        for chunks in iterate_chunks( gffs ):
            ninput += 1
            if len(chunks) > 1:
                gene_id = "merged_%s" % chunks[0][0].gene_id
                transcript_id = "merged_%s" % chunks[0][0].transcript_id
                info = ",".join([ x[0].gene_id for x in chunks ])
            else:
                gene_id = chunks[0][0].gene_id
                transcript_id = chunks[0][0].transcript_id
                info = None

            intervals = []
            for c in chunks: intervals += [ (x.start, x.end) for x in c ]
            intervals = Intervals.combine( intervals )

            for start, end in intervals:
                y = GTF.Entry()
                y.fromGFF( chunks[0][0], gene_id, transcript_id )
                y.start = start
                y.end = end
                y.strand = "."
                if info: y.addAttribute( "merged", info )
                options.stdout.write( "%s\n" % str(y ) )
                nfeatures += 1

            noutput += 1
        
    elif options.transcripts2genes:

        transcripts = set()
        genes = set()
        reset_strand = options.reset_strand
        for gtfs in GTF.iterator_transcripts2genes( GTF.iterator(options.stdin) ):

            ninput += 1
            for gtf in gtfs:
                if reset_strand: gtf.strand = "."
                options.stdout.write( "%s\n" % str(gtf) )                
                transcripts.add( gtf.transcript_id )
                genes.add( gtf.gene_id )
                nfeatures += 1
            noutput += 1

        if options.loglevel >= 1:
            options.stdlog.write("# transcripts2genes: transcripts=%i, genes=%i\n" % \
                                     (len(transcripts), len(genes)))

    elif options.rename:

        map_old2new = IOTools.readMap( open(options.filename_filter, "r") )

        if options.rename == "transcript":
            is_gene_id = False
        elif options.rename == "gene":
            is_gene_id = True
            
        for gff in GTF.iterator( options.stdin ):
            ninput += 1

            if is_gene_id and gff.gene_id in map_old2new:
                gff.gene_id = map_old2new[gff.gene_id]
            elif not is_gene_id and gff.transcript_id in map_old2new:
                gff.transcript_id = map_old2new[gff.transcript_id]
            else:
                ndiscarded += 1
                continue
            noutput += 1
            options.stdout.write("%s\n" % str(gff))
                
    elif options.filter:

        keep_genes = set()
        if options.filter == "longest-gene":
            iterator = GTF.flat_gene_iterator( GTF.iterator(options.stdin) )
            coords = []
            gffs = []
            for gff in iterator:
                gff.sort( key = lambda x: x.start )
                coords.append( (gff[0].contig,
                                min([x.start for x in gff]), 
                                max([x.end for x in gff]),
                                gff[0].gene_id ) )
                gffs.append( gff )
            coords.sort()
            
            last_contig = None
            max_end = 0
            longest_gene_id = None
            longest_length = None

            for contig, start, end, gene_id in coords:
                ninput += 1
                if contig != last_contig or start >= max_end:
                    if longest_gene_id: keep_genes.add( longest_gene_id )
                    longest_gene_id = gene_id
                    longest_length = end - start
                    max_end = end
                else:
                    if end-start > longest_length:
                        longest_length, longest_gene_id = end-start, gene_id
                last_contig = contig
                max_end = max(max_end, end)

            keep_genes.add(longest_gene_id)
            invert = options.invert_filter
            for gff in gffs:
                keep = gff[0].gene_id in keep_genes

                if (keep and not invert) or (not keep and invert):
                    noutput += 1
                    for g in gff:
                        nfeatures += 1
                        options.stdout.write("%s\n" % g )
                else:
                    ndiscarded += 1

        elif options.filter in ("gene", "transcript"):
            if options.filename_filter:

                ids, nerrors = IOTools.ReadList( open(options.filename_filter, "r") )

                ids = set(ids)
                by_gene = options.filter == "gene"
                by_transcript = options.filter == "transcript"
                invert = options.invert_filter

                reset_strand = options.reset_strand
                for gff in GTF.iterator(options.stdin):

                    ninput += 1

                    keep = False
                    if by_gene: keep = gff.gene_id in ids
                    if by_transcript: keep = gff.transcript_id in ids
                    if (invert and keep) or (not invert and not keep): continue

                    if reset_strand: gff.strand = "."

                    options.stdout.write( "%s\n" % str(gff) )                
                    nfeatures += 1
                    noutput += 1

            elif options.sample_size:

                if options.filter == "gene":
                    iterator = GTF.flat_gene_iterator( GTF.iterator(options.stdin) )
                elif options.filter == "transcript":
                    iterator = GTF.transcript_iterator( GTF.iterator(options.stdin) )
                if options.min_exons_length:
                    iterator = GTF.iterator_min_feature_length( iterator, 
                                                                min_length = options.min_exons_length,
                                                                feature = "exon" )

                data = [ x for x in iterator ]
                ninput = len(data)
                if len(data) > options.sample_size:
                    data = random.sample( data, options.sample_size )

                for d in data:
                    noutput += 1
                    for dd in d:
                        nfeatures += 1
                        options.stdout.write( str(dd) + "\n" )

            else:
                assert False, "please supply either a filename with ids to filter with (--apply) or a sample-size."
        
    elif options.exons2introns:

        for gffs in GTF.flat_gene_iterator(GTF.iterator(options.stdin)):

            ninput += 1

            cds_ranges = GTF.asRanges( gffs, "CDS" )
            exon_ranges = GTF.asRanges( gffs, "exon" )
            input_ranges = Intervals.combine( cds_ranges + exon_ranges )
            
            if len(input_ranges) > 1:
                last = input_ranges[0][1]
                output_ranges = []
                for start, end in input_ranges[1:]:
                    output_ranges.append( (last, start) )
                    last = end

                for start, end in output_ranges:
                    
                    entry = GTF.Entry()
                    entry.copy( gffs[0] )
                    entry.clearAttributes()
                    entry.transcript_id = "merged"
                    entry.feature = "intron"
                    entry.start = start
                    entry.end = end
                    options.stdout.write( "%s\n" % str( entry ) )
                    nfeatures += 1                    
                noutput += 1
            else:
                ndiscarded += 1

    elif options.set_score2distance:

        for gffs in GTF.transcript_iterator(GTF.iterator(options.stdin)):
            ninput += 1
            strand = Genomics.convertStrand( gffs[0].strand )
            all_start, all_end = min( [ x.start for x in gffs ] ), max( [x.end for x in gffs ] ) 

            if strand != ".":
                t = 0            
                if strand == "-": gffs.reverse()
                for gff in gffs:
                    gff.score = t
                    t += gff.end - gff.start

                if strand == "-": gffs.reverse()
            for gff in gffs:
                options.stdout.write( "%s\n" % str( gff ) )
                nfeatures += 1
            noutput += 1
    else:
        for gffs in GTF.flat_gene_iterator(GTF.iterator(options.stdin)):

            ninput += 1

            cds_ranges = GTF.asRanges( gffs, "CDS" )
            exon_ranges = GTF.asRanges( gffs, "exon" )
            strand = Genomics.convertStrand( gffs[0].strand )

            if cds_ranges and options.with_utr:
                cds_start, cds_end = cds_ranges[0][0], cds_ranges[-1][1]
                midpoint = ( cds_end - cds_start ) / 2 + cds_start

                utr_ranges = []
                for start, end in Intervals.truncate( exon_ranges, cds_ranges ):
                    if end - start > 3:
                        if strand == ".":
                            feature = "UTR"
                        elif strand == "+":
                            if start < midpoint:
                                feature = "UTR5"
                            else:
                                feature = "UTR3"
                        elif strand == "-":
                            if start < midpoint:
                                feature = "UTR3"
                            else:
                                feature = "UTR5"
                        utr_ranges.append( (feature, start,end) )
                output_feature = "CDS"
                output_ranges = cds_ranges
            else:
                output_feature = "exon"
                output_ranges = exon_ranges
                utr_ranges = []

            result = []

            if options.merge_exons:

                for feature, start, end in utr_ranges:
                    entry = GTF.Entry()
                    entry.copy( gffs[0] )
                    entry.clearAttributes()
                    entry.feature = feature
                    entry.transcript_id = "merged"
                    entry.start = start
                    entry.end = end
                    result.append( entry )

                for start, end in output_ranges:

                    entry = GTF.Entry()
                    entry.copy( gffs[0] )
                    entry.clearAttributes()
                    entry.transcript_id = "merged"
                    entry.feature = output_feature
                    entry.start = start
                    entry.end = end
                    result.append( entry )

            elif options.merge_transcripts:

                entry = GTF.Entry()
                entry.copy( gffs[0] )
                entry.clearAttributes()
                entry.transcript_id = entry.gene_id
                entry.start = output_ranges[0][0]
                entry.end = output_ranges[-1][1]
                result.append( entry )

            elif options.merge_introns:

                if len( output_ranges ) >= 2:
                    entry = GTF.Entry()
                    entry.copy( gffs[0] )
                    entry.clearAttributes()
                    entry.transcript_id = entry.gene_id
                    entry.start = output_ranges[0][1]
                    entry.end = output_ranges[-1][0]
                    result.append( entry )
                else:
                    ndiscarded += 1
                    continue

            result.sort( key=lambda x: x.start )

            for x in result:
                options.stdout.write( "%s\n" % str(x) )
                nfeatures += 1
            noutput += 1

    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nfeatures=%i, ndiscarded=%i\n" % (ninput, noutput, nfeatures, ndiscarded) )
    E.Stop()
