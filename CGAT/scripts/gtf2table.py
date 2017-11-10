'''gtf2table.py - annotate genes/transcripts
============================================

:Tags: Genomics Genesets GTF Annotation

Purpose
-------

annotate genes or transcripts given in a :term:`gtf` formatted file
and output them in tabular format.

Annotations can be either computed per gene (all exons across all
transcripts) or per transcript. The input needs to be sorted
accordingly.

The script iterates over each gene or transcript model in turn and
outputs a row in a table for each. Multiple annotations can be
computed at the same time adding additional columns to the table.

For example, to output information about exons in genes::

    zcat in.gtf.gz | cgat gtf2table -v 0 --counter=length

+---------------+----+----+----+---------+------+--------+----+----+----+
|gene_id        |nval|min |max |mean     |median|stddev  |sum |q1  |q3  |
+---------------+----+----+----+---------+------+--------+----+----+----+
|ENSG00000225373|4   |39  |688 |253.0000 |142.5 |261.5922|1012|58  |688 |
+---------------+----+----+----+---------+------+--------+----+----+----+
|ENSG00000267111|1   |744 |744 |744.0000 |744.0 |0.0000  |744 |744 |744 |
+---------------+----+----+----+---------+------+--------+----+----+----+
|ENSG00000267588|3   |135 |489 |370.0000 |486.0 |166.1746|1110|135 |489 |
+---------------+----+----+----+---------+------+--------+----+----+----+
|ENSG00000220978|1   |138 |138 |138.0000 |138.0 |0.0000  |138 |138 |138 |
+---------------+----+----+----+---------+------+--------+----+----+----+

The table contains exon length statistics of each gene, the number (nval),
min, max, and total length(sum) of exons per gene. To count per transcript::

    zcat in.gtf.gz | cgat gtf2table -v 0 --reporter=transcripts --counter=length

+---------------+----+---+---+--------+------+--------+----+---+---+
|transcript_id  |nval|min|max|mean    |median|stddev  |sum |q1 |q3 |
+---------------+----+---+---+--------+------+--------+----+---+---+
|ENST00000592209|3   |58 |227|149.6667|164.0 |69.7344 |449 |58 |227|
+---------------+----+---+---+--------+------+--------+----+---+---+
|ENST00000589741|2   |71 |312|191.5000|191.5 |120.5000|383 |71 |312|
+---------------+----+---+---+--------+------+--------+----+---+---+
|ENST00000391654|2   |39 |583|311.0000|311.0 |272.0000|622 |39 |583|
+---------------+----+---+---+--------+------+--------+----+---+---+
|ENST00000587045|1   |173|173|173.0000|173.0 |0.0000  |173 |173|173|
+---------------+----+---+---+--------+------+--------+----+---+---+

To add also information about cpg-composition, add another counter::

    zcat in.gtf.gz | cgat gtf2table -v 0 --genome=hg19 --reporter=transcripts --counter=length --counter=composition-cpg

+---------------+----+---+---+--------+------+--------+----+---+---+----------+------------+-----------+
|transcript_id  |nval|min|max|mean    |median|stddev  |sum |q1 |q3 |CpG_count |CpG_density |CpG_ObsExp |
+---------------+----+---+---+--------+------+--------+----+---+---+----------+------------+-----------+
|ENST00000592209|3   |58 |227|149.6667|164.0 |69.7344 |449 |58 |227|4         |0.01781     |0.13946    |
+---------------+----+---+---+--------+------+--------+----+---+---+----------+------------+-----------+
|ENST00000589741|2   |71 |312|191.5000|191.5 |120.5000|383 |71 |312|5         |0.02610     |0.17067    |
+---------------+----+---+---+--------+------+--------+----+---+---+----------+------------+-----------+
|ENST00000391654|2   |39 |583|311.0000|311.0 |272.0000|622 |39 |583|4         |0.01286     |0.09402    |
+---------------+----+---+---+--------+------+--------+----+---+---+----------+------------+-----------+
|ENST00000587045|1   |173|173|173.0000|173.0 |0.0000  |173 |173|173|1         |0.01156     |0.11396    |
+---------------+----+---+---+--------+------+--------+----+---+---+----------+------------+-----------+

Note that we had to use the ``--genome-file`` option to supply the
genomic sequence.

Additional switches permit counting introns instead of exons (option
``--section``).

Annotations
-----------

The annotations might be derived from a property of the transcript or
gene model itself (such as length, number of exons, etc), the genomic
sequence (such as composition, etc) or require additional data sets.
For example, to compute the coverage with reads, a :term:`bam` formatted
file is required (see option ``--bam-file``). Or to compute the overlap
with genomic densities, a :term:`bigwig` formatted file is required (see
option ``--bigwig-file``).

This section lists the counters available. They are grouped by data sources
they require.

Generic annotations
+++++++++++++++++++

Generic annotations require no additional data source to annotate a
transcript or gene.

length
   output exon length summary of transcript/gene. This counter outputs the
   number of exons in each transcript/gene together with exon length summary
   statistics (minimal exon length, maximal exon length, total exon length).

position
   output genomic coordinates of transcript/gene (chromosome, start, end).

Sequence derived annotations
++++++++++++++++++++++++++++

Sequence derived annotations require the genomic sequence to compute
properties of the gene/transcript model (see option ``--genome-file``).

composition-na
   output nucleotide composition of transcript/gene

composition-cgp
   output CpG count, CpG density and CpG observed / expected for
   each transcript/gene.

splice
   output splicing summary of transcript/gene. Outputs the number of
   canonical and non-canonical splice sites.

quality
   output base-quality information summary of gene. Needs quality scores.

Interval derived annotations
++++++++++++++++++++++++++++

Annotations in the list below relate a gene or transcript to a set of
intervals given as a second file. The second set of intervals is given
by the option ``filename-gff``. By default, the intervals are expected
to be given as a :term:`gff` formatted file, but alternative formats
(:term:`gtf` and :term:`bed`) are possible (see option
``--filename-format``)

overlap
    compute overlap of genes in input with features in second stream.
    Requires a :term:`gff` formatted file with gene territories.

overlap-stranded
    count overlap with genomic features in second another file. Outputs
    the number of overlapping exons. Records the direction of overlap
    (sense/antisense). Requires a :term:`gff` formatted file with
    features.

territories
    compute overlap of transcripts/genes with territories. Territories
    are genomic, non-overlapping intervals. For each transcript the
    counter outputs the status (U=unique match, A=ambiguous match
    (overlap with multiple territories, 0=no match), the number of
    territories it overlaps and the territories it overlaps.

coverage

   compute the coverage - per nucleotide - of the gene/transcript
   models with intervals given in ``--gff-file``.  Coverage values
   are output in 5' to 3' together with summary statistics (bases
   covered, minimum, maximum coverage, etc.). By using the options
   ``--restrict-feature`` or ``--restrict-source`` the counting can be
   rescricted to particular features or sources in the :term:`gff` file.

distance

   compute distance of genes to features in a second file. Requires a
   second :term:`gff` formatted file. The strand information of the
   features is ignored.

binding-pattern

   given a list of intervals, determine the binding pattern within and
   surrounding the gene. For each gene, intervals overlapping the CDS,
   introns, UTRs and the flank are collected and recorded. The binding
   is summarized with a binding pattern, a binary pattern indicating
   overlap/no overlap with 5' flank, 5' UTR, CDS, Introns, 3' UTR, 3'
   flank. This method is useful to check where transcription factor binding
   sites are located aronud a gene/transcript model.

proximity

   report summary stats (lengths, values) of features in proximity to
   genes input gene set. Requires a :term:`gff` formatted file with
   genomic features. This feature is useful when aiming to normalize a
   value, such as a substitution rate of a transcript model, by
   substitution rates of segments in the neighbourhood such as
   ancestral repeats. The values are given in the ``score`` field of
   the :term:`gff` formatted file. The radius for proximity is controlled
   by the option ``--proximal-distance``.

proximity-exclusive

   as proximity, but exclude any ranges in the :term:`gff` formatted file that
   overlap the transcript/gene model.

proximity-lengthmatched

   as proximity-exclusive, but length-match features with
   genes. Segments are declared equal in length if they are within
   10% of the original segments length.

neighbours
    output features in second stream that are in proximity to genes
    in input. This is similar to the ``proximity`` counters, but also
    outputs the features that are in proximity.

Gene set derived annotations
++++++++++++++++++++++++++++

distance-genes
   compute distance of genes to genes in a second file. Requires a
   second :term:`gtf` formatted file with genes. The counter distinguishes
   a variety of cases (closest upstream/downstream).

distance-tss
   compute distance of genes to transcription start sites. Requires a
   second :term:`gtf` formatted file with genes.

overlap-transcripts
    count overlap of genes with transcripts in another set.
    Requires a :term:`gtf` formatted file.

overrun
   output intron overrun, exons in the input gene set extending
   into the introns of a reference gene set. Requries a :term:`gtf`
   formatted file with a reference gene set.

splice-comparison
   Compare how splice site usage compares between a gene/transcript
   with transcripts in a reference gene set. Outputs found, missed,
   perfect, partial, incomplete splice sites and exon-skipping events.


Short-read derived annotations
++++++++++++++++++++++++++++++

Short-read derived annotations count the overlap of reads given
in a :term:`bam` formatted file (see option ``--bam-file``) with
gene or transcript models.

read-overlap

   output number of reads overlapping a transcript/gene model.
   Outputs the number of reads overlapping the transcript/gene model.
   Reads are counted separately for sense and antisense overlap.

read-coverage

   output read coverage summary statistics of transcript/gene model.
   Outputs the number of reads overlapping the transcript/gene model,
   the number of bases overlapped by reads and summary statistics of
   the minimimum, maximum, etc. coverage per pase. Reads are counted
   separately for sense and antisense overlap. The counter does not
   take into account splice sites. As it does per-base counting,
   it is slower than ``read-overlap`` and there is no need to
   use both at the same time.

read-extension

   Counter of special interest. This counter outputs the read density
   in bins upstream, within and downstream of transcript models. The
   counter can be used to predict the length of the 3' and 5' UTR.

read-fullcounts

   count number of reads overlapping a gene or transcript. Reads
   are classified according to the type of overlap (sense/antisense),
   full or partial match to exon, splice status (spliced/unspliced).
   See below for more information.

   Unique and non-unique matches are counted (by alignment start
   position).

read-counts

   count number of reads overlapping a gene or transcript. Summarizes
   the output of read-fullcounts.

readpair-fullcounts

   count number of read pairs overlapping a gene or transcript. Pairs
   are counted according to a variety of attributes (exonic
   overlap/read pair status/splice status/...). See below for more
   information.

readpair-counts

   count number of read pairs overlapping a gene or
   transcript. Summarizes the output of readpair-fullcounts.


Classifiers
-----------

Classifiers not only annotate the transcripts or gene model, but also
aim to provide some classification based on these annotations. They
require a secondary file (see option ``--gff-file``) for the
classification.

classifier

   classify transcripts according to genomic annotation.  Requires a
   :term:`gff` file with genomic annotations (see :doc:`gtf2gff`)
   delineating genomic regions as intronic, intergenic, exonic,
   etc. This is useful for a rough classification of transcripts as
   intergenic, intronic, etc.

   Best used for ChIP-Seq data sets (each peak is a "transcript"). The
   method is a little out of place in this script, but is here as it
   uses much of the code implemented here.

classifier-rnaseq

   classify transcripts with respect to a reference geneset. The
   classifiers aims to match a transcript up with the "most similar"
   transcript in the reference gene set and dependending on
   attributes, classifies it as a good match, alternative transcript,
   fragment, etc. The classifier takes into account strand and prefers
   sense matches to antisense matches.

classifier-polii

   classify according to PolII transcripts. A gene/transcript is
   transcribed, if it is covered by large PolII intervals over 80% of
   its length. A gene/transript is primed if its promotor/UTR is
   covered by 50% of its length, while the rest of the gene body
   isn't.

Other annotations
++++++++++++++++++

The following methods (see option ``--counter``) are available:

bigwig-counts

   collect density values from a :term:`bigwig` formatted file and output
   summary statistics and percentage of bases covered (``pcovered``)
   by value in bigwig file. Requires option ``--bigwig-file``.


Read counting
-------------

The methods ``read-counts`` and ``readpair-counts`` count the number
of reads inside a :term:`bam` formatted file that map inside a
transcript or gene.

These counters proceed on a per-gene or per-transcript basis depending
on the ``--reporter`` option. All reads overlapping the exons or
introns of a transcript are collected and counted either individually
(``read-counts``) or as pairs (``readpair-counts``).

For paired-read counting, each pair is classified and counted
according to the following for axes:

1. Pair status: proper pair, unmapped read in pair, ...
2. Direction: for stranded protocols, sense, antisense, ...
3. Overlap status: Pair overlaps exons only, introns only, ...
4. Splice status: Reads are unspliced, spliced consistently with
   transcript model, or inconsistently with transcript model

Counts are then collated into a smaller set of summaries for convenience:

+--------------------+----------------------------------------------------+
|Column              |Number of read pairs                                |
+--------------------+----------------------------------------------------+
|counted_all         |considered to be correct: they are sense, overlap   |
|                    |with exons fully.                                   |
+--------------------+----------------------------------------------------+
|counted_spliced     |in sense direction and exonic, spliced and spliced  |
|                    |correctly.                                          |
|                    |                                                    |
+--------------------+----------------------------------------------------+
|counted_unspliced   |in sense direction and exonic, not spliced.         |
+--------------------+----------------------------------------------------+
|sense_intronic      |intronic and sense direction.                       |
+--------------------+----------------------------------------------------+
|sense_inconsistent  |in sense direction, but overlap both introns and    |
|                    |exons or extend beyond the transcript modelse       |
+--------------------+----------------------------------------------------+
|sense_other         |in sense direction, but not counted, intronic or    |
|                    |inconsistent.                                       |
+--------------------+----------------------------------------------------+
|antisense           |in antisense direction.                             |
+--------------------+----------------------------------------------------+
|nonsense            |in unexpected orientation                           |
+--------------------+----------------------------------------------------+
|notproper           |that are not proper pairs                           |
+--------------------+----------------------------------------------------+
|quality_pairs       |affected by a the removal of a low-quality read     |
+--------------------+----------------------------------------------------+
|quality_reads       |Number of low quality reads removed                 |
+--------------------+----------------------------------------------------+
|total               |Total number of pairs considered                    |
+--------------------+----------------------------------------------------+

The Counter ``readpair-fullcounts`` provides the detailed counts
for each transcript model according to the four axes. The output
can be used to implement custom counting schemes.

Reads below a minimum quality score will be ignored
(``--min-read-quality``). By default, all reads will be counted.

Reads mapping to multiple locations will be downweighted if
the ``--weight-multi-mapping`` option is set. This requires
the presence of the ``NH`` flag in the :term:`bam` file.

For paired read counting, the library type can be specified with the
``--library-type`` option to make use of strand information. Library
types are labelled according to the tophat_ and cufflinks_
convention. A summary is `here
<http://www.nature.com/nprot/journal/v7/n3/fig_tab/nprot.2012.016_T1.html>`_.

+--------------------+--------------------+------------------------------+
|Library type        |RNA-seq protocol    |Explanation                   |
+--------------------+--------------------+------------------------------+
|fr-unstranded       |Illumina TruSeq     |Reads from the leftmost end of|
| (default)          |                    |the fragment (in transcript   |
|                    |                    |coordinates) map to the       |
|                    |                    |transcript strand, and the    |
|                    |                    |rightmost end maps to the     |
|                    |                    |opposite strand left          |
|                    |                    |                              |
+--------------------+--------------------+------------------------------+
|fr-firststrand      |dUTP, NSR, NNSR39   |Same as above except we       |
|                    |                    |enforce the rule that the     |
|                    |                    |rightmost end of the fragment |
|                    |                    |(in transcript coordinates) is|
|                    |                    |the first sequenced (or only  |
|                    |                    |sequenced for single-end      |
|                    |                    |reads). Equivalently, it is   |
|                    |                    |assumed that only the strand  |
|                    |                    |generated during first strand |
|                    |                    |synthesis is sequenced        |
+--------------------+--------------------+------------------------------+
|fr-secondstrand     |Directional Illumina|Same as above                 |
|                    |(Ligation)          |except TopHat/Cufflinks       |
|                    |Standard SOLiD      |enforce the rule that the     |
|                    |                    |leftmost end of the fragment  |
|                    |                    |(in transcript coordinates) is|
|                    |                    |the first sequenced (or only  |
|                    |                    |sequenced for single-end      |
|                    |                    |reads). Equivalently, it is   |
|                    |                    |assumed that only the strand  |
|                    |                    |generated during second strand|
|                    |                    |synthesis is sequenced        |
+--------------------+--------------------+------------------------------+

For unstranded protocols, all reads and pairs are considered to matching
in the sense direction.

Usage
-----

Type::

   python gtf2table.py --help

for command line help.

Command line options
--------------------

'''

import sys
import pysam

import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.IndexedFasta as IndexedFasta
import CGAT.GeneModelAnalysis as GeneModelAnalysis

import pyBigWig


def main(argv=None):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default].")

    parser.add_option("-q", "--quality-file",
                      dest="quality_file",
                      type="string",
                      help="filename with genomic base quality "
                      "information [default=%default].")

    parser.add_option("-b", "--bam-file", dest="bam_files",
                      type="string", metavar="bam",
                      help="filename with read mapping information. "
                      "Multiple files can be submitted in a "
                      "comma-separated list [default=%default].")

    parser.add_option("-i", "--bigwig-file", dest="bigwig_file",
                      type="string", metavar="bigwig",
                      help="filename with bigwig information "
                      "[default=%default].")

    parser.add_option("-f", "--gff-file", dest="filename_gff",
                      type="string", action="append", metavar='bed',
                      help="filename with extra gff files. The order "
                      "is important [default=%default].")

    parser.add_option("--filename-format", dest="filename_format",
                      type="choice",
                      choices=("bed", "gff", "gtf"),
                      help="format of secondary stream [default=%default].")

    parser.add_option("--restrict-source", dest="gff_sources", type="string",
                      action="append",
                      help="restrict input to this 'source' in extra "
                      "gff file (for counter: overlap) [default=%default].")

    parser.add_option("--restrict-feature", dest="gff_features", type="string",
                      action="append",
                      help="restrict input to this 'feature' in extra gff "
                      "file (for counter: overlap) [default=%default].")

    parser.add_option("-r", "--reporter", dest="reporter", type="choice",
                      choices=("genes", "transcripts"),
                      help="report results for 'genes' or 'transcripts' "
                      "[default=%default].")

    parser.add_option("-s", "--section", dest="sections",
                      type="choice",
                      action="append",
                      choices=("exons", "introns"),
                      help="select range on which counters will operate "
                      "[default=%default].")

    parser.add_option("-c", "--counter", dest="counters",
                      type="choice",
                      action="append",
                      choices=(	"bigwig-counts",
                                "binding-pattern",
                                "classifier",
                                "classifier-rnaseq",
                                "classifier-rnaseq-splicing",
                                "classifier-polii",
                                "composition-na",
                                "composition-cpg",
                                "coverage",
                                "distance",
                                "distance-genes",
                                "distance-tss",
                                "length",
                                'neighbours',
                                "overlap",
                                "overlap-stranded",
                                "overlap-transcripts",
                                "overrun",
                                "position",
                                "proximity",
                                "proximity-exclusive",
                                "proximity-lengthmatched",
                                "quality",
                                "read-coverage",
                                "read-extension",
                                "read-overlap",
                                "read-counts",
                                "read-fullcounts",
                                "readpair-counts",
                                "readpair-fullcounts",
                                "splice",
                                "splice-comparison",
                                "territories"),
                      help="select counters to apply to input "
                      "[default=%default].")

    parser.add_option("--add-gtf-source", dest="add_gtf_source",
                      action="store_true",
                      help="add gtf field of source to output "
                      "[default=%default].")

    parser.add_option("--proximal-distance", dest="proximal_distance",
                      type="int",
                      help="distance to be considered proximal to "
                      "an interval [default=%default].")

    parser.add_option("--multi-mapping-method",
                      dest="multi_mapping",
                      type="choice",
                      choices=('all', 'ignore', 'weight'),
                      help="how to treat multi-mapping reads in "
                      "bam-files. Requires "
                      "the NH flag to be set by the mapper "
                      "[default=%default].")

    parser.add_option("--use-barcodes",
                      dest="use_barcodes",
                      action="store_true",
                      help="Use barcodes to count unique umi's. "
                      "UMI's are specified in the read identifier "
                      "as the last field, where fields are separated "
                      "by underscores, e.g. "
                      "@READ:ILLUMINA:STUFF_NAMINGSTUFF_UMI. "
                      "When true, unique counts are returned. "
                      "Currently only compatible with count-reads")

    parser.add_option("--sample-probability",
                      dest="sample_probability",
                      type="float",
                      help="Specify the probability of whether any"
                      "given read or read pair in a file bam is counted"
                      "Currently only compatible with count-reads")

    parser.add_option("--column-prefix", dest="prefixes",
                      type="string",
                      action="append",
                      help="add prefix to column headers - prefixes "
                      "are used in the same order as the counters "
                      "[default=%default].")

    parser.add_option("--library-type",
                      dest="library_type",
                      type="choice",
                      choices=("unstranded",
                               "firststrand",
                               "secondstrand",
                               "fr-unstranded",
                               "fr-firststrand",
                               "fr-secondstrand"),
                      help="library type of reads in bam file. "
                      "[default=%default]")

    parser.add_option("--min-mapping-quality",
                      dest="minimum_mapping_quality",
                      type="float",
                      help="minimum mapping quality. Reads with a quality "
                      "score of less will be ignored. "
                      "[default=%default]")

    parser.set_defaults(
        genome_file=None,
        reporter="genes",
        with_values=True,
        sections=[],
        counters=[],
        filename_gff=[],
        filename_format=None,
        gff_features=[],
        gff_sources=[],
        add_gtf_source=False,
        proximal_distance=10000,
        bam_files=None,
        multi_mapping='all',
        library_type='fr-unstranded',
        prefixes=[],
        minimum_mapping_quality=0,
        use_barcodes=False,
        sample_probability=1.0
    )

    if not argv:
        argv = sys.argv

    (options, args) = E.Start(parser, add_output_options=True, argv=argv)

    if options.prefixes:
        if len(options.prefixes) != len(options.counters):
            raise ValueError(
                "if any prefix is given, the number of prefixes "
                "must be the same as the number of counters")

    # get files
    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta(options.genome_file)
    else:
        fasta = None

    if options.quality_file:
        quality = IndexedFasta.IndexedFasta(options.quality_file)
        quality.setTranslator(IndexedFasta.TranslatorBytes())
    else:
        quality = None

    if options.bam_files:
        bam_files = []
        for bamfile in options.bam_files.split(","):
            bam_files.append(pysam.AlignmentFile(bamfile, "rb"))
    else:
        bam_files = None

    if options.bigwig_file:
        bigwig_file = pyBigWig.open(options.bigwig_file)
    else:
        bigwig_file = None

    counters = []

    if not options.sections:
        E.info("counters will use the default section (exons)")
        options.sections.append(None)

    if not options.gff_sources:
        options.gff_sources.append(None)
    if not options.gff_features:
        options.gff_features.append(None)

    cc = E.Counter()

    for n, c in enumerate(options.counters):
        if options.prefixes:
            prefix = options.prefixes[n]
        else:
            prefix = None

        if c == "position":
            for section in options.sections:
                counters.append(
                    GeneModelAnalysis.CounterPosition(
                        section=section,
                        options=options,
                        prefix=prefix))
        elif c == "length":
            for section in options.sections:
                counters.append(
                    GeneModelAnalysis.CounterLengths(
                        section=section,
                        options=options,
                        prefix=prefix))
        elif c == "splice":
            if fasta is None:
                raise ValueError('splice requires a genomic sequence')
            counters.append(GeneModelAnalysis.CounterSpliceSites(fasta=fasta, prefix=prefix))
        elif c == "quality":
            if fasta is None:
                raise ValueError('quality requires a quality score sequence')
            counters.append(GeneModelAnalysis.CounterQuality(fasta=quality, prefix=prefix))
        elif c == "overrun":
            counters.append(GeneModelAnalysis.CounterOverrun(
                filename_gff=options.filename_gff,
                options=options,
                prefix=prefix))
        elif c == "read-coverage":
            counters.append(GeneModelAnalysis.CounterReadCoverage(
                bam_files,
                options=options,
                prefix=prefix))
        elif c == "read-extension":
            counters.append(GeneModelAnalysis.CounterReadExtension(
                bam_files,
                filename_gff=options.filename_gff,
                options=options,
                prefix=prefix))
        elif c == "read-overlap":
            counters.append(GeneModelAnalysis.CounterReadOverlap(
                bam_files,
                multi_mapping=options.multi_mapping,
                minimum_mapping_quality=options.minimum_mapping_quality,
                options=options,
                prefix=prefix))
        elif c == "read-counts":
            counters.append(GeneModelAnalysis.CounterReadCounts(
                bam_files,
                multi_mapping=options.multi_mapping,
                use_barcodes=options.use_barcodes,
                sample_probability=options.sample_probability,
                minimum_mapping_quality=options.minimum_mapping_quality,
                options=options,
                prefix=prefix))
        elif c == "read-fullcounts":
            counters.append(GeneModelAnalysis.CounterReadCountsFull(
                bam_files,
                multi_mapping=options.multi_mapping,
                sample_probability=options.sample_probability,
                minimum_mapping_quality=options.minimum_mapping_quality,
                options=options,
                prefix=prefix))
        elif c == "readpair-counts":
            counters.append(GeneModelAnalysis.CounterReadPairCounts(
                bam_files,
                multi_mapping=options.multi_mapping,
                sample_probability=options.sample_probability,
                library_type=options.library_type,
                minimum_mapping_quality=options.minimum_mapping_quality,
                options=options,
                prefix=prefix))
        elif c == "readpair-fullcounts":
            counters.append(GeneModelAnalysis.CounterReadPairCountsFull(
                bam_files,
                multi_mapping=options.multi_mapping,
                sample_probability=options.sample_probability,
                minimum_mapping_quality=options.minimum_mapping_quality,
                options=options,
                prefix=prefix))
        elif c == "bigwig-counts":
            counters.append(GeneModelAnalysis.CounterBigwigCounts(
                bigwig_file,
                options=options, prefix=prefix))
        elif c == "splice-comparison":
            if fasta is None:
                raise ValueError('splice-comparison requires a genomic '
                                 'sequence')
            counters.append(GeneModelAnalysis.CounterSpliceSiteComparison(
                fasta=fasta,
                filename_gff=options.filename_gff,
                feature=None,
                source=None,
                options=options, prefix=prefix))
        elif c == "composition-na":
            if fasta is None:
                raise ValueError('composition-na requires a genomic sequence')
            for section in options.sections:
                counters.append(GeneModelAnalysis.CounterCompositionNucleotides(
                    fasta=fasta,
                    section=section,
                    options=options,
                    prefix=prefix))
        elif c == "composition-cpg":
            if fasta is None:
                raise ValueError('composition-cpg requires a genomic sequence')
            for section in options.sections:
                counters.append(GeneModelAnalysis.CounterCompositionCpG(
                    fasta=fasta,
                    section=section,
                    options=options, prefix=prefix))

        elif c in ("overlap",
                   "overlap-stranded",
                   "overlap-transcripts",
                   "proximity",
                   "proximity-exclusive",
                   "proximity-lengthmatched",
                   "neighbours",
                   "territories",
                   "distance",
                   "distance-genes",
                   "distance-tss",
                   "binding-pattern",
                   "coverage"):
            if c == "overlap":
                template = GeneModelAnalysis.CounterOverlap
            if c == "overlap-stranded":
                template = GeneModelAnalysis.CounterOverlapStranded
            elif c == "overlap-transcripts":
                template = GeneModelAnalysis.CounterOverlapTranscripts
            elif c == "proximity":
                template = GeneModelAnalysis.CounterProximity
            elif c == "neighbours":
                template = GeneModelAnalysis.CounterNeighbours
            elif c == "proximity-exclusive":
                template = GeneModelAnalysis.CounterProximityExclusive
            elif c == "proximity-lengthmatched":
                template = GeneModelAnalysis.CounterProximityLengthMatched
            elif c == "territories":
                template = GeneModelAnalysis.CounterTerritories
            elif c == "distance":
                template = GeneModelAnalysis.CounterDistance
            elif c == "distance-genes":
                template = GeneModelAnalysis.CounterDistanceGenes
            elif c == "distance-tss":
                template = GeneModelAnalysis.CounterDistanceTranscriptionStartSites
            elif c == "coverage":
                template = GeneModelAnalysis.CounterCoverage
            elif c == "binding-pattern":
                template = GeneModelAnalysis.CounterBindingPattern

            for section in options.sections:
                for source in options.gff_sources:
                    for feature in options.gff_features:
                        counters.append(template(
                            filename_gff=options.filename_gff,
                            feature=feature,
                            source=source,
                            fasta=fasta,
                            section=section,
                            options=options,
                            prefix=prefix))

        elif c == "classifier":
            counters.append(GeneModelAnalysis.Classifier(
                filename_gff=options.filename_gff,
                fasta=fasta,
                options=options, prefix=prefix))

        elif c == "classifier-rnaseq":
            counters.append(GeneModelAnalysis.ClassifierRNASeq(
                filename_gff=options.filename_gff,
                fasta=fasta,
                options=options, prefix=prefix))
        elif c == "classifier-rnaseq-splicing":
            counters.append(GeneModelAnalysis.ClassifierRNASeqSplicing(
                filename_gff=options.filename_gff,
                fasta=fasta,
                options=options,
                prefix=prefix))
        elif c == "classifier-polii":
            counters.append(GeneModelAnalysis.ClassifierPolII(
                filename_gff=options.filename_gff,
                feature=None,
                source=None,
                fasta=fasta,
                options=options,
                prefix=prefix))
        elif c == "binding-pattern":
            counters.append(GeneModelAnalysis.CounterBindingPattern(
                filename_gff=options.filename_gff,
                feature=None,
                source=None,
                fasta=fasta,
                options=options,
                prefix=prefix))

    if options.reporter == "genes":
        iterator = GTF.flat_gene_iterator
        header = ["gene_id"]
        fheader = lambda x: [x[0].gene_id]
    elif options.reporter == "transcripts":
        iterator = GTF.transcript_iterator
        header = ["transcript_id"]
        fheader = lambda x: [x[0].transcript_id]

    if options.add_gtf_source:
        header.append("source")
        ffields = lambda x: [x[0].source]
    else:
        ffields = lambda x: []

    options.stdout.write("\t".join(
        header + [x.getHeader() for x in counters]) + "\n")

    for gffs in iterator(GTF.iterator(options.stdin)):
        cc.input += 1

        for counter in counters:
            counter.update(gffs)

        skip = len([x for x in counters if x.skip]) == len(counters)
        if skip:
            cc.skipped += 1
            continue

        options.stdout.write("\t".join(
            fheader(gffs) +
            ffields(gffs) +
            [str(counter) for counter in counters]) + "\n")

        cc.output += 1

    E.info("%s" % str(cc))
    for counter in counters:
        E.info("%s\t%s" % (repr(counter), str(counter.counter)))
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
