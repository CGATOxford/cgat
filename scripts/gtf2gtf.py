'''
gtf2gtf.py - manipulate transcript models
=========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics Genesets GTF Manipulation

Purpose
-------

This script reads a gene set in :term:`gtf` format from stdin, applies some 
transformation, and outputs a new gene set in :term:`gtf` format to stdout. 

Options
-------

Options available for use in this script can broadly be classified into four
categories:

1. sorting gene sets
2. manipulating gene models
3. filtering gene sets
4. setting/resetting fields within a gtf file

Further options for working with gtf files are available in gff2gff.py, 
which can be run with the specification --is-gtf


Sort gene sets
++++++++++++++

``--sort``
   Sorts entries in gtf file by one or more fields      

      +---------------+---------------------------------------+
      | option        | order in which fields are sorted      | 
      +---------------|---------------------------------------+
      | gene          | gene_id, transcript_id, contig, start |
      +---------------+---------------------------------------+
      | contig+gene   | contig, gene_id, transcript_id, start | 
      +---------------+---------------------------------------+
      | transcript    | transcript_id, contig, start          |
      +---------------+---------------------------------------+
      | position      | contig, start                         |
      +---------------+---------------------------------------+
      | position+gene | contig( gene_id, start )              |
      +---------------+---------------------------------------+
      | gene+position | gene_id, contig, start                | 
      +---------------+---------------------------------------+
   
   N.B. position+gene sorts by gene_id, start, then subsequently sorts
   flattened gene lists by contig, start


Manipulate gene-models
++++++++++++++++++++++

Options that can be used to alter the features represented in a :term:`gtf`
file. For further detail see command line options.

Input gtfs need to be sorted so that features for a gene or transcript 
appear consecutively within the file. This can be achevied using ``--sort``.

``--merge-exons``
    Merges overlapping exons for all transcripts of a gene, outputting the 
    merged exons. Can be used in conjunction with ``--merge-exons-distance`` 
    to set the minimum distance that may appear between two exons before 
    they are merged.If ``--with-utr`` is set, the output interval will also
    contain UTR.

``--merge-transcripts``
    Merges all transcripts of a gene. Outputs contains a single interval that
    spans the original gene (both introns and exons). If ``--with-utr`` is 
    set, the output interval will also contain UTR.

``--merge-genes``
    Merges genes that have overlapping exons, outputting a single gene_id and 
    transcript_id for all exons of overlapping genes. 
    (Ignores strand information.)
    
``--join-exons``
    Joins together all exons of a transcript, outputting a single interval that 
    spans the original transcript (both introns and exons).

``--intersect-transcripts``
    Finds regions representing the intersect of all transcripts of a gene.
    Output will contain intervals spanning only those bases covered by all 
    transcripts. If ``--with-utr`` is set, the UTR will also be included in the 
    intersect.

``--merge-introns``
    Merges the region spanned by introns for all transcripts of a gene.
    Outputs a single interval that spans the region between the start and end of
    the first and last intron, respectively.

``--exons2introns``
    Merges overlapping introns for all transcripts of a gene, outputting the 
    merged introns. Use ``--intron-min-length`` to ignore merged introns below a
    specified length. Use ``--intron-border`` to specify a number of residues to
    remove at either end of output introns (residues are removed prior to 
    filtering on size when used in conjunction with ``--intron-min-length``).
    
``--transcripts2genes``    
    Currently it doesn't work... 
    pysam.TabProxies.GTFProxy has no attribute 'gene_id'    
    May be used in conjunction with ``--reset-strand``    

The option ``--permit-duplicates`` may be specified in order to
allow gene-ids to be duplicated within the input :term:`gtf` file
(i.e. for the same gene-id to appear non-consecutively within the
input file). However, this option currently only works for
``--merge-exons``, ``--merge-transcripts``, ``--merge-introns``, and
``--intersect-transcripts``. It DOES NOT work for ``--merge-genes``,
``--join-exons``, or ``--exons2introns``.

Filter gene sets
++++++++++++++++

Options that can be used to filter :term:`gtf` files. For further
detail see command line options.

Input gtfs need to be sorted so that features for a gene or transcript 
 appear consecutively within the file. This can be achevied using ``--sort``.

``--filter``
    When filtering on the basis of 'gene-id' or 'transcript-id' a filename 
    containing ids to be removed may provided using ``--apply``. Alternatively,
    a random subsample of genes/transcripts may be retained using 
    ``--sample-size``. Use ``--min-exons-length`` in conjunction with 
    ``--sample-size`` to specify a minimum length for genes/transcripts to be 
    retained. Use ``--reset-strand`` to set strand to '.' in output.

    Other filter options include longest-gene, longest-transcript, 
    or representative-transcript.

    When filtering on the basis of gene-id, transcript-id or longest-gene,
    ``--invert-filter`` may be used to invert the selection.

``--remove-overlapping``
    Takes as argument a :term:`gff` formatted file. Any transcripts that 
    intersect intervals in the supplied file are removed. 
    (Does not account for strand.)

``--remove-duplicates``
    Setting to 'gene', 'transcript', or 'coordinates' will remove any interval for
    which non-consecutive occurrances of specified term appear in input :term:`gtf`
    file.
    Setting to 'ucsc', will remove any interval for which transcript-id contains
    '_dup'.


Set/reset fields
++++++++++++++++

Options for altering fields within :term:`gtf`. For further details see command 
 line options.

``--rename``
    When a mapping file is provided using ``--apply``, renames either 
    gene-id or transcript-id. Outputs a :term:`gtf` file with field renamed. Any
    entry in input :term:`gtf` not appearing in mapping file is discarded. 

``--add-protein-id``
    Takes as argument a file that maps transcript-id to protein-id and appends the
    protein-id provided to the attributes field.
    Any entry ininput :term:`gtf` not appearing in mapping file is discarded. 

``--renumber-genes``
    Renumber genes from 1 using the pattern provided, e.g. ENSG%s

``--renumber-transcripts``
    Renumber transcripts from 1 using the pattern provided, e.g. ENST%s

``--unset-genes``
    Renumber genes from 1 using the pattern provided, e.g. ENSG%s. Transcripts 
    with the same gene-id in input :term:`gtf` file will have different gene-ids
    in the output :term:`gtf` file.

``--set-transcript-to-gene``
    Will reset the transcript-id to the gene-id

``--set-gene-to-transcript``
    Will reset the gene-id to the transcript-id

``--set-protein-to-transcript``
    Will append transcript to attributes field as 'protein_id'

``--set-score-to-distance``
    Will reset the score field (field 6) of each feature in input :term:`gtf` to be
    the distance from transcription start site to the start of the feature. 
    (Assumes input file is sorted by transcript-id)


Usage
-----

Example::

    python gtf2gtf.py --sort=gene | python gtf2gtf.py --intersect-transcripts --with-utr --renumber-transcripts= MERGED%s

Type::

    python gtf2gtf.py --help

Command line Options
--------------------

'''
import os
import sys
import string
import re
import optparse
import types
import random
import collections
import itertools

import CGAT.GTF as GTF
import CGAT.Experiment as E
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Genomics as Genomics
import CGAT.Intervals as Intervals
import CGAT.IOTools as IOTools
import CGAT.Components as Components

##------------------------------------------------------------
## This script needs some attention.
##------------------------------------------------------------
def main( argv = None ):

    if not argv: argv = sys.argv

    parser = E.OptionParser( version = "%prog version: $Id: gtf2gtf.py 2861 2010-02-23 17:36:32Z andreas $", 
                             usage = globals()["__doc__"])

    parser.add_option("-m", "--merge-exons", dest="merge_exons", action="store_true",
                      help="merge overlapping exons of all transcripts within a gene. "
                           "The merged exons will be output. "
                           "Input needs to sorted by gene [default=%default]."  )

    parser.add_option("-t", "--merge-transcripts", dest="merge_transcripts", action="store_true",
                      help="merge all transcripts within a gene. The entry will span the whole gene (exons and introns). "
                           "The transcript does not include the UTR unless --with-utr is set. [default=%default]."  )

    parser.add_option( "--merge-genes", dest="merge_genes", action="store_true",
                       help="merge overlapping genes if their exons overlap. "
                       "A gene with a single transcript containing all exons "
                       "of the overlapping transcripts will be output. "
                       "This operation ignores strand information "
                       "The input needs te sorted by transcript [default=%default]."  )

    parser.add_option( "--merge-exons-distance", dest="merge_exons_distance", type="int",
                      help="distance in nucleotides between exons to be merged [default=%default]."  )

    parser.add_option("-j", "--join-exons", dest="join_exons", action="store_true",
                      help="join all exons per transcript. A new transcript will be "
                           "output that spans a whole transcript. "
                           "Input needs to be sorted by transcript [default=%default]."  )

    parser.add_option( "--unset-genes", dest="unset_genes", type="string",
                       help="unset gene identifiers, keeping transcripts intact. "
                       "New gene identifiers are set to the pattern given. For example, "
                       "'--unset-genes=%06i' [default=%default]."  )

    parser.add_option( "--sort", dest="sort", type="choice",
                       choices=("gene", 
                                "transcript", 
                                "position", 
                                "contig+gene", 
                                "position+gene", 
                                "gene+position" ),
                       help="sort input data [default=%default]."  )

    parser.add_option("-u", "--with-utr", dest="with_utr", action="store_true",
                      help="include utr in merged transcripts [default=%default]."  )


    parser.add_option( "--intersect-transcripts", dest="intersect_transcripts", action="store_true",
                      help="intersect all transcripts within a gene. The entry will only span those bases "
                      "that are covered by all transcrips."
                      "The transcript does not include the UTR unless --with-utr is set. This method"
                      "will remove all other features (stop_codon, etc.) "
                      "The input needs to be sorted by gene. "
                      "[default=%default]."  )

    parser.add_option("-i", "--merge-introns", dest="merge_introns", action="store_true",
                      help="merge and output all introns within a gene. The output will contain "
                      "all intronic regions within a gene. Single exon genes "
                      "are skipped. "
                      "The input needs to be sorted by gene. "
                      "[default=%default]."  )

    parser.add_option("-g", "--set-transcript-to-gene", "--set-transcript2gene", dest="set_transcript2gene", action="store_true",
                      help="set the transcript_id to the gene_id [default=%default]."  )

    parser.add_option( "--set-protein-to-transcript", dest="set_protein2transcript", action="store_true",
                      help="set the protein_id to the transcript_id [default=%default]."  )

    parser.add_option( "--add-protein-id", dest="add_protein_id", type="string",
                      help="add a protein_id for each transcript_id. "
                           "The argument is a filename containing a mapping between "
                           "transcript_id to protein_id [default=%default]."  )

    parser.add_option("-G", "--set-gene-to-transcript", "--set-gene2transcript", dest="set_gene2transcript", action="store_true",
                      help="set the gene_id to the transcript_id [default=%default]."  )

    parser.add_option("-d", "--set-score2distance", dest="set_score2distance", action="store_true",
                      help="set the score field for each feature to the distance to "
                           "transcription start site [default=%default]."  )

    parser.add_option( "--exons2introns", dest="exons2introns", action="store_true",
                       help="for each gene build an 'intronic' transcript "
                       "containing the union of all intronic regions "
                       "of all transcripts in a gene."
                       "The features are labeled as 'intron'." 
                       "The input needs to be sorted by gene. "
                      "[default=%default]."  )

    parser.add_option("-f", "--filter", dest="filter", type="choice",
                      choices=("gene", "transcript", "longest-gene", "longest-transcript", "representative-transcript" ),
                      help="apply a filter to the input file. Available filters are: "
                      "'gene': filter by gene_id, "
                      "'transcript': filter by transcript_id, "
                      "'longest-gene': output the longest gene for overlapping genes ," 
                      "'longest-transcript': output the longest transcript per gene," 
                      "'representative-transcript': output the representative transcript per gene. "
                      "The representative transcript is the transcript that shares most exons with "
                      "the other transcripts in a gene. "
                      "The input needs to be sorted by gene. "
                      "[default=%default]."  )
                      

    parser.add_option("-r", "--rename", dest="rename", type="choice",
                      choices=("gene", "transcript" ),
                      help="rename genes or transcripts with a map given by the option `--apply`. "
                      "Those that can not be renamed are removed "
                      "[default=%default]."  )

    parser.add_option( "--renumber-genes", dest="renumber_genes", type="string", 
                       help="renumber genes according to the given pattern. "
                       "[default=%default]."  )
    
    parser.add_option( "--renumber-transcripts", dest="renumber_transcripts", type="string", 
                       help="renumber transcripts according to the given pattern. "
                       "[default=%default]."  )

    parser.add_option("-a", "--apply", dest="filename_filter", type="string",
                      metavar="tsv",
                      help="filename of ids to map/filter [default=%default]." )
    
    parser.add_option( "--invert-filter", dest="invert_filter", action="store_true",
                       help="when using --filter, invert selection (like grep -v). "
                       "[default=%default]."  )

    parser.add_option("--sample-size", dest="sample_size", type="int",
                      help="extract a random sample of size # if the option "
                      "'--filter' is set[default=%default]." )

    parser.add_option("--intron-min-length", dest="intron_min_length", type="int",
                      help="minimum length for introns (for --exons2introns) "
                      "[default=%default]."  )

    parser.add_option("--min-exons-length", dest="min_exons_length", type="int",
                      help="minimum length for gene (sum of exons) (--sample-size) [default=%default]."  )

    parser.add_option("--intron-border", dest="intron_border", type="int",
                      help="number of residues to exclude at intron at either end "
                      "(--exons2introns) [default=%default]."  )

    parser.add_option( "--transcripts2genes", dest="transcripts2genes", action="store_true",
                       help="cluster overlapping transcripts into genes." )

    parser.add_option( "--reset-strand", dest="reset_strand", action="store_true",
                       help="remove strandedness of features (set to '.') when using --transcripts2genes"
                       "[default=%default]."  )

    parser.add_option( "--remove-overlapping", dest="remove_overlapping", type="string",
                       metavar = "gff",
                       help="remove all transcripts that overlap intervals "
                       "in a gff-formatted file."
                       "The comparison ignores strand "
                       "[default=%default]." )

    parser.add_option( "--permit-duplicates", dest="strict", action="store_false",
                       help="permit duplicate genes. "
                       "[default=%default]" )

    parser.add_option( "--remove-duplicates", dest="remove_duplicates", type="choice",
                       choices=("gene", "transcript", "ucsc", "coordinates"),
                       help="remove duplicates by gene/transcript. "
                       "If ``ucsc`` is chosen, transcripts ending on _dup# are "
                       "removed. This is necessary to remove duplicate entries "
                       "that are next to each other in the sort order [%default]" )

    parser.set_defaults(
        sort = None,
        merge_exons = False,
        join_exons = False,
        merge_exons_distance = 0,
        merge_transcripts = False,
        set_score2distance = False,
        set_gene2transcript = False,
        set_transcript2gene = False,
        set_protein2transcript = False,
        add_protein_id = None,
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
        remove_duplicates = None,
        remove_overlapping = None,
        renumber_genes = None,
        unset_genes = None,
        renumber_transcripts = None,
        strict = True,
        intersect_transcripts = False,
        )

    (options, args) = E.Start( parser, argv = argv )
    
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

        if options.remove_duplicates == "ucsc":
            store = []
            remove = set()
            f = lambda x: x[0].transcript_id

            gffs = GTF.transcript_iterator(GTF.iterator(options.stdin), strict = False)
            outf = lambda x: "\n".join( [ str(y) for y in x] )

            for entry in gffs:
                ninput += 1
                store.append(entry)
                id = f(entry)
                if "_dup" in id: 
                    remove.add(re.sub("_dup\d+","", id) )
                    remove.add( id )

            for entry in store:
                id = f(entry)
                if id not in remove: 
                    options.stdout.write( outf(entry) + "\n" )
                    noutput += 1
                else:
                    ndiscarded += 1
                    E.info("discarded duplicates for %s" % (id))
        else:
        
            if options.remove_duplicates == "gene":
                gffs = GTF.gene_iterator(GTF.iterator(options.stdin), strict = False )
                f = lambda x: x[0][0].gene_id
                outf = lambda x: "\n".join( [ "\n".join( [ str(y) for y in xx] ) for xx in x] )
            elif options.remove_duplicates == "transcript":
                gffs = GTF.transcript_iterator(GTF.iterator(options.stdin), strict = False)
                f = lambda x: x[0].transcript_id
                outf = lambda x: "\n".join( [ str(y) for y in x] )
            elif options.remove_duplicates == "coordinates":
                gffs = GTF.chunk_iterator(GTF.iterator(options.stdin) )
                f = lambda x: x[0].contig+"_"+str(x[0].start)+"-"+str(x[0].end)
                outf = lambda x: "\n".join( [ str(y) for y in x] )    

            store = []

            for entry in gffs:
                ninput += 1
                store.append(entry)
                id = f(entry)
                counts[id] += 1

            # Assumes GTF file sorted by contig then start
            last_id = ""
            if options.remove_duplicates == "coordinates":
                for entry in store:
                    id = f(entry)
                    if id == last_id:
                        ndiscarded += 1
                        E.info("discarded duplicates for %s: %i" % (id, counts[id]))
                    else:
                        options.stdout.write( outf(entry) + "\n" )
                        noutput += 1
                    last_id = id
                        
            else:
                for entry in store:
                    id = f(entry)
                    if counts[id] == 1:
                        options.stdout.write( outf(entry) + "\n" )
                        noutput += 1
                    else:
                        ndiscarded += 1
                        E.info("discarded duplicates for %s: %i" % (id, counts[id]))

    elif options.sort:

        for gff in GTF.iterator_sorted( GTF.iterator( options.stdin ), sort_order = options.sort ):
            ninput += 1
            options.stdout.write( "%s\n" % str(gff) )                
            noutput += 1
            nfeatures += 1
            
    elif options.set_gene2transcript:

        for gff in GTF.iterator(options.stdin):

            ninput += 1

            gff.gene_id = gff.transcript_id
            options.stdout.write( "%s\n" % str(gff) )                
            noutput += 1
            nfeatures += 1

    elif options.set_protein2transcript:

        for gff in GTF.iterator(options.stdin):
            ninput += 1
            gff.setAttribute( "protein_id", gff.transcript_id)
            options.stdout.write( "%s\n" % str(gff) )                
            noutput += 1
            nfeatures += 1

    elif options.add_protein_id:

        transcript2protein = IOTools.readMap( open( options.add_protein_id, "r") )

        missing = set()
        for gff in GTF.iterator(options.stdin):
            ninput += 1
            if gff.transcript_id not in transcript2protein:
                if gff.transcript_id not in missing:
                    E.debug( "removing transcript '%s' due to missing protein id" % gff.transcript_id)
                    missing.add( gff.transcript_id)
                ndiscarded += 1
                continue
            
            gff.setAttribute( "protein_id", transcript2protein[gff.transcript_id])
            options.stdout.write( "%s\n" % str(gff) )                
            noutput += 1
            nfeatures += 1
        
        E.info("transcripts removed due to missing protein ids: %i" % len(missing))

    elif options.join_exons:

        for exons in GTF.transcript_iterator( GTF.iterator(options.stdin) ):
            ninput += 1
            strand = Genomics.convertStrand( exons[0].strand )
            contig = exons[0].contig
            transid = exons[0].transcript_id
            geneid = exons[0].gene_id
            biotype = exons[0].source
            all_start, all_end = min( [ x.start for x in exons ] ), max( [x.end for x in exons ] )
            y = GTF.Entry()
            y.contig = contig
            y.source = biotype
            y.feature = "transcript"
            y.start = all_start
            y.end = all_end
            y.strand = strand
            y.transcript_id = transid
            y.gene_id = geneid
            options.stdout.write( "%s\n" % str(y ) )

    elif options.merge_genes:
        # merges overlapping genes
        # 
        gffs = GTF.iterator_sorted_chunks( 
            GTF.flat_gene_iterator(GTF.iterator(options.stdin)),
            sort_by = "contig-strand-start" )
        
        def iterate_chunks( gff_chunks ):

            last = gff_chunks.next()
            to_join = [last]

            for gffs in gff_chunks:
                d = gffs[0].start - last[-1].end

                if gffs[0].contig == last[0].contig and gffs[0].strand == last[0].strand:
                    assert gffs[0].start >= last[0].start, \
                        "input file should be sorted by contig, strand and position: d=%i:\nlast=\n%s\nthis=\n%s\n" % \
                        (d, 
                         "\n".join( [str(x) for x in last] ),
                         "\n".join( [str(x) for x in gffs] ) )

                if gffs[0].contig != last[0].contig or \
                        gffs[0].strand != last[0].strand or \
                        d > 0:
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
            # take single strand
            strand = chunks[0][0].strand

            for start, end in intervals:
                y = GTF.Entry()
                y.fromGTF( chunks[0][0], gene_id, transcript_id )
                y.start = start
                y.end = end
                y.strand = strand

                if info: y.addAttribute( "merged", info )
                options.stdout.write( "%s\n" % str(y ) )
                nfeatures += 1

            noutput += 1

    elif options.renumber_genes:
        
        map_old2new = {}
        for gtf in GTF.iterator(options.stdin):
            ninput += 1
            if gtf.gene_id not in map_old2new:
                map_old2new[gtf.gene_id ] = options.renumber_genes % (len(map_old2new) + 1)
            gtf.setAttribute("gene_id", map_old2new[gtf.gene_id ] )
            options.stdout.write( "%s\n" % str(gtf) )
            noutput += 1

    elif options.unset_genes:
        
        map_old2new = {}
        for gtf in GTF.iterator(options.stdin):
            ninput += 1
            key = gtf.transcript_id
            if key not in map_old2new:
                map_old2new[key] = options.unset_genes % (len(map_old2new) + 1)
            gtf.setAttribute( "gene_id",  map_old2new[key] )
            options.stdout.write( "%s\n" % str(gtf) )
            noutput += 1

    elif options.renumber_transcripts:
        
        map_old2new = {}
        for gtf in GTF.iterator(options.stdin):
            ninput += 1
            key = (gtf.gene_id, gtf.transcript_id)
            if key not in map_old2new:
                map_old2new[key] = options.renumber_transcripts % (len(map_old2new) + 1)
            gtf.setAttribute( "transcript_id",  map_old2new[key] )
            options.stdout.write( "%s\n" % str(gtf) )
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

            if is_gene_id:
                if gff.gene_id in map_old2new:
                    gff.setAttribute("gene_id", map_old2new[gff.gene_id])
                else:
                    E.debug( "removing missing gene_id %s" % gff.gene_id)
                    ndiscarded += 1
                    continue

            else:
                if gff.transcript_id in map_old2new:
                    gff.setAttribute("transcript_id", map_old2new[gff.transcript_id])
                else:
                    E.debug( "removing missing transcript_id %s" % gff.transcript_id)
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
        elif options.filter in ("longest-transcript", "representative-transcript" ):

            iterator = GTF.gene_iterator( GTF.iterator(options.stdin) )

            def selectLongestTranscript( gene ):
                r = []
                for transcript in gene:
                    transcript.sort( key = lambda x: x.start )
                    length = transcript[-1].end - transcript[0].start
                    r.append( (length, transcript) )
                r.sort()
                return r[-1][1]
            
            def selectRepresentativeTranscript( gene ):
                '''select a representative transcript.

                The representative transcript represent the largest number
                of exons over all transcripts.
                '''
                all_exons = []
                for transcript in gene:
                    all_exons.extend( [ (x.start, x.end) for x in transcript if x.feature == "exon"] )
                exon_counts = {}
                for key, exons in itertools.groupby( all_exons ):
                    exon_counts[key] = len(list(exons))
                transcript_counts = []
                for transcript in gene:
                    count = sum( [ exon_counts[(x.start, x.end)] for x in transcript if x.feature == "exon" ] )
                    transcript_counts.append( (count, transcript) )
                transcript_counts.sort()
                return transcript_counts[-1][1]

            if options.filter == "longest-transcript":
                _select = selectLongestTranscript
            elif options.filter == "representative-transcript":
                _select = selectRepresentativeTranscript

            for gene in iterator:
                ninput += 1
                transcript = _select( gene )
                noutput += 1
                for g in transcript:
                    nfeatures += 1
                    options.stdout.write("%s\n" % g )

        elif options.filter in ("gene", "transcript"):

            if options.filename_filter:
                
                ids, nerrors = IOTools.ReadList( open(options.filename_filter, "r") )
                E.info( "read %i ids" % len(ids))

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

                if options.intron_border:
                    b = options.intron_border
                    output_ranges = [ (x[0]+b, x[1]-b ) for x in output_ranges ]

                if options.intron_min_length:
                    l = options.intron_min_length
                    output_ranges = [ x for x in output_ranges if x[1]-x[0] > l ]
            
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

    elif options.remove_overlapping:
        
        index = GTF.readAndIndex( GTF.iterator( IOTools.openFile( options.remove_overlapping, "r" ) ) )
        
        for gffs in GTF.transcript_iterator(GTF.iterator(options.stdin)):
            ninput += 1
            found = False
            for e in gffs:
                if index.contains( e.contig, e.start, e.end ):
                    found = True
                    break
            
            if found:
                ndiscarded += 1    
            else:
                noutput += 1
                for e in gffs: 
                    nfeatures += 1
                    options.stdout.write( "%s\n" % str(e ) )

    elif options.intersect_transcripts:
        
        for gffs in GTF.gene_iterator(GTF.iterator(options.stdin), strict=options.strict ):
            
            ninput += 1
            r = []
            for g in gffs:
                if options.with_utr:
                    ranges = GTF.asRanges( g, "exon" )
                else:
                    ranges = GTF.asRanges( g, "CDS" )
                r.append( ranges )

            result = r[0]
            for x in r[1:]:
                result = Intervals.intersect( result, x )
            
            entry = GTF.Entry()
            entry.copy( gffs[0][0] )
            entry.clearAttributes()
            entry.transcript_id = "merged"
            entry.feature = "exon"
            for start, end in result:
                entry.start = start
                entry.end = end
                options.stdout.write( "%s\n" % str( entry ) )
                nfeatures += 1
                
            noutput += 1
    else:
        for gffs in GTF.flat_gene_iterator(GTF.iterator(options.stdin), strict=options.strict ):

            ninput += 1

            cds_ranges = GTF.asRanges( gffs, "CDS" )
            exon_ranges = GTF.asRanges( gffs, "exon" )

            # sanity checks
            strands = set( [x.strand for x in gffs ] )
            contigs = set( [x.contig for x in gffs ] )
            if len(strands) > 1:
                raise ValueError( "can not merge gene '%s' on multiple strands: %s" % (gffs[0].gene_id, str(strands)))

            if len(contigs) > 1:
                raise ValueError( "can not merge gene '%s' on multiple contigs: %s" % (gffs[0].gene_id, str(contigs)))

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
                # need to combine per feature - skip
                # utr_ranges = Intervals.combineAtDistance( utr_ranges, options.merge_exons_distance )

                output_ranges = Intervals.combineAtDistance( output_ranges, options.merge_exons_distance )

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

    E.info("ninput=%i, noutput=%i, nfeatures=%i, ndiscarded=%i" % (ninput, noutput, nfeatures, ndiscarded) )
    E.Stop()


if __name__ == "__main__":
    sys.exit( main( sys.argv ))
