"""========================================
RNA-Seq Differential expression pipeline
========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

The RNA-Seq differential expression pipeline performs differential
expression analysis. It requires three inputs:

   1 one or more genesets in :term:`gtf` formatted files
   2 mapped reads in :term:`bam` formatted files
   3 design files as :term:`tsv`-separated format

This pipeline works on a single genome.

Overview
========

The pipeline performs the following:

   * compute tag counts on an exon, transcript and gene level for each
     of the genesets. The following tag counting methods are implemented:

      * featureCounts_
      * gtf2table

   * estimate FPKM using cufflinks for each of the gene sets

   * perform differential expression analysis. The methods currently
     implemented are:

      * DESeq
      * EdgeR
      * Cuffdiff

Background
============

The quality of the rnaseq data (read-length, paired-end) determines
the quality of transcript models. For instance, if reads are short
(35bp) and/or reads are not paired-ended, transcript models will be
short and truncated.  In these cases it might be better to concentrate
the analysis on only previously known transcript models.

Transcripts are the natural choice to measure expression of. However
other quantities might be of interest. Some quantities are biological
meaningful, for example differential expression from a promotor shared
by several trancripts. Other quantities might no biologically
meaningful but are necessary as a technical comprise.  For example,
the overlapping transcripts might be hard to resolve and thus might
need to be aggregated per gene. Furthermore, functional annotation is
primarily associated with genes and not individual transcripts. The
pipeline attempts to measure transcription and differential expression
for a variety of entities following the classification laid down by
:term:`cuffdiff`:

isoform
   Transcript level
gene
   Gene level, aggregates several isoform/transcripts
tss
   Transcription start site. Aggregate all isoforms starting from the
   same :term:`tss`.
cds
   Coding sequence expression. Ignore reads overlapping non-coding
   parts of transcripts (UTRs, etc.). Requires annotation of the cds
   and thus only available for :file:`reference.gtf.gz`.

Methods differ in their ability to measure transcription on all
levels.

Overprediction of differential expression for low-level expressed
transcripts with :term:`cuffdiff` is a `known problem
<http://seqanswers.com/forums/showthread.php?t=6283&highlight=fpkm>`_.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.

The sphinxreport report requires a :file:`conf.py` and
:file:`sphinxreport.ini` file (see :ref:`PipelineReporting`). To start
with, use the files supplied with the Example_ data.

Input
-----

Reads
+++++

Reads are imported by placing :term:`bam` formatted files are linking
to files in the :term:`working directory`.

The default file format assumes the following convention::

   <samplename>.bam

Genesets
++++++++

Genesets are imported by placing :term:`gtf` formatted files in the
:term:`working directory`.

Design matrices
+++++++++++++++

Design matrices are imported by placing :term:`tsv` formatted files
into the :term:`working directory`. A design matrix describes the
experimental design to test. Each design file has four columns::

      track   include group   pair
      CW-CD14-R1      0       CD14    1
      CW-CD14-R2      0       CD14    1
      CW-CD14-R3      1       CD14    1
      CW-CD4-R1       1       CD4     1
      FM-CD14-R1      1       CD14    2
      FM-CD4-R2       0       CD4     2
      FM-CD4-R3       0       CD4     2
      FM-CD4-R4       0       CD4     2

track
     name of track - should correspond to a sample name.
include
     flag to indicate whether or not to include this data
group
     group indicator - experimental group
pair
     pair that sample belongs to (for paired tests) - set to 0 if the
     design is not paired.

Optional inputs
+++++++++++++++

Requirements
------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

+-----------+----------+----------------------------+
|*Program*  |*Version* |*Purpose*                   |
+-----------+----------+----------------------------+
|cufflinks_ |>=1.3.0   |transcription levels        |
+-----------+----------+----------------------------+
|samtools   |>=0.1.16  |bam/sam files               |
+-----------+----------+----------------------------+
|bedtools   |          |working with intervals      |
+-----------+----------+----------------------------+
|R/DESeq    |          |differential expression     |
+-----------+----------+----------------------------+

Pipeline output
===============

FPKM values
-----------

The directory :file:`fpkm.dir` contains the output of running
cufflinks_ against each set of reads for each geneset.

The syntax of the output files is ``<geneset>_<track>.cufflinks``. The
FPKM values are uploaded into tables as ``<geneset>_<track>_fpkm`` for
per-transcript FPKM values and ``<geneset>_<track>_genefpkm`` for
per-gene FPKM values.

Differential gene expression results
-------------------------------------

Differential expression is estimated for different genesets with a
variety of methods. Results are stored per method in subdirectories
such as :file:`deseq.dir`, :file:`edger.dir` or :file:`cuffdiff.dir`.

The syntax of the output files is ``<design>_<geneset>_...``.

:file:`.diff` files:
    The .diff files are :term:`tsv` formatted tables with the result
    of the differential expression analysis.

    The column names are:

    test_id
       gene_id/transcript_id/tss_id, etc
    treatment_name
       name of the treatment
    control_name
       name of the control
    !=_mean
       mean expression value of treatment/control
    !=_std
       standard deviation of treatment/control
    pvalue
       pvalue of test
    qvalue
       qvalue of test (FDR if the threshold was set at pvalue)
    l2fold
       log2 of fold change
    fold
       fold change, treatment / control
    significant
       flag: 1 if called significant
    status
       status of test. Values are
       OK: test ok, FAIL: test failed,
       NOCALL: no test (insufficient data, etc.).

The results are uploaded into the database as tables called
``<design>_<geneset>_<method>_<level>_<section>``.

Level denotes the biological entity that the differential analysis was
performed on.  Level can be one of ``cds``,``isoform``,``tss`` and
``gene``. The later should always be present.

Section is ``diff`` for differential expression results
(:file:`.diff`` files) and ``levels`` for expression levels.


Example
=======

Example data is available at
http://www.cgat.org/~andreas/sample_data/pipeline_rnaseqdiffexpression.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_rnaseqdiffexpression.tgz
   tar -xvzf pipeline_rnaseq.tgz
   cd pipeline_rnaseq
   python <srcdir>/pipeline_rnaseq.py make full

.. note::

   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::

   cufflinks
      cufflinks_ - transcriptome analysis

   tophat
      tophat_ - a read mapper to detect splice-junctions

   deseq
      deseq_ - differential expression analysis

   cuffdiff
      find differentially expressed transcripts. Part of cufflinks_.

   cuffcompare
      compare transcriptomes. Part of cufflinks_.

.. _cufflinks: http://cufflinks.cbcb.umd.edu/index.html
.. _tophat: http://tophat.cbcb.umd.edu/
.. _bamstats: http://www.agf.liv.ac.uk/454/sabkea/samStats_13-01-2011
.. _deseq: http://www-huber.embl.de/users/anders/DESeq/
.. _featurecounts: http://bioinf.wehi.edu.au/featureCounts/

ChangeLog
=========

28.03.2014  Andreas Heger
            added automated selection of paired counting to featureCounts
            and gtf2table counting.

11.4.2014   Andreas Heger
            changed workflow. Multiple counters are applied and
            differential expression is computed on all.


Code
====

"""

# load modules
from ruffus import *
from ruffus.combinatorics import *

import CGAT.Experiment as E
import CGAT.Database as Database

import sys
import os
import re
import itertools
import glob
import collections
import sqlite3
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
from rpy2.robjects import r as R
import rpy2.robjects as ro

import CGAT.Expression as Expression
import CGAT.BamTools as BamTools
import CGATPipelines.PipelineGeneset as PipelineGeneset
import CGATPipelines.PipelineRnaseq as PipelineRnaseq

# levels of cuffdiff analysis
# (no promotor and splice -> no lfold column)
CUFFDIFF_LEVELS = ("gene", "cds", "isoform", "tss")

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# load options from the config file
import CGAT.Pipeline as P
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters(PARAMS["annotations_dir"],
                                      "pipeline_annotations.py", on_error_raise=__name__ == "__main__")

###################################################################
###################################################################
# Helper functions mapping tracks to conditions, etc
###################################################################
import CGATPipelines.PipelineTracks as PipelineTracks

Sample = PipelineTracks.AutoSample

# collect sra nd fastq.gz tracks
TRACKS = PipelineTracks.Tracks(Sample).loadFromDirectory(
    glob.glob("*.bam"), "(\S+).bam")

# group by experiment (assume that last field is a replicate identifier)
EXPERIMENTS = PipelineTracks.Aggregate(TRACKS, labels=("condition", "tissue"))

GENESETS = PipelineTracks.Tracks(Sample).loadFromDirectory(
    glob.glob("*.gtf.gz"), "(\S+).gtf.gz")

###################################################################
###################################################################
###################################################################


def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh

#########################################################################
#########################################################################
#########################################################################
# Definition of targets


# Expression levels: geneset vs bam files
TARGETS_FPKM = [(("%s.gtf.gz" % x.asFile(), "%s.bam" % y.asFile()),
                 "%s_%s.cufflinks" % (x.asFile(), y.asFile()))
                for x, y in itertools.product(GENESETS, TRACKS)]

#########################################################################
#########################################################################
#########################################################################
# preparation targets

#########################################################################
#########################################################################
#########################################################################


@files(os.path.join(PARAMS["annotations_dir"], "geneset_all.gtf.gz"),
       "geneset_mask.gtf")
def buildMaskGtf(infile, outfile):
    '''This takes ensembl annotations (geneset_all.gtf.gz) and writes out
    all entries that have a 'source' match to "rRNA" or 'contig' match
    to "chrM". for use with cufflinks
    '''
    geneset = IOTools.openFile(infile)
    outf = open(outfile, "wb")
    for entry in GTF.iterator(geneset):
        if re.findall("rRNA", entry.source) or re.findall(
                "chrM", entry.contig):
            outf.write("\t".join((map(
                str,
                [entry.contig, entry.source, entry.feature,
                 entry.start, entry.end, ".", entry.strand,
                 ".", "transcript_id" + " " + '"' +
                 entry.transcript_id + '"' + ";" + " " +
                 "gene_id" + " " + '"' + entry.gene_id + '"'])))
                + "\n")

    outf.close()

#########################################################################
#########################################################################
#########################################################################


@transform("*.gtf.gz",
           suffix(".gtf.gz"),
           "_geneinfo.load")
def loadGeneSetGeneInformation(infile, outfile):
    PipelineGeneset.loadGeneStats(infile, outfile)

#########################################################################


@follows(mkdir("fpkm.dir"))
@files([(x, os.path.join("fpkm.dir", y)) for x, y in TARGETS_FPKM])
def runCufflinks(infiles, outfile):
    '''estimate expression levels in each set using cufflinks.'''
    PipelineRnaseq.runCufflinks(infiles, outfile)


@transform(runCufflinks,
           suffix(".cufflinks"),
           ".load")
def loadCufflinks(infile, outfile):
    '''load expression level measurements.'''
    PipelineRnaseq.loadCufflinks(infile, outfile)


@collate(runCufflinks,
         regex("fpkm.dir/(.*)_(.*).cufflinks"),
         r"fpkm.dir/\1_fpkm_genes.tsv.gz")
def mergeCufflinksGeneFPKM(infiles, outfile):
    '''build aggregate table with cufflinks FPKM values.'''

    prefix = os.path.basename(outfile)
    prefix = prefix[:prefix.index("_")]

    headers = ",".join(
        [re.match("fpkm.dir/.*_(.*).cufflinks", x).groups()[0]
         for x in infiles])

    statement = '''
    python %(scriptsdir)s/combine_tables.py
        --columns=1
        --skip-titles
        --headers=%(headers)s
        --take=FPKM fpkm.dir/%(prefix)s_*.genes_tracking.gz
    | perl -p -e "s/tracking_id/gene_id/"
    | gzip
    > %(outfile)s
    '''
    P.run()


@collate(runCufflinks,
         regex("fpkm.dir/(.*)_(.*).cufflinks"),
         r"fpkm.dir/\1_fpkm_isoforms.tsv.gz")
def mergeCufflinksIsoformFPKM(infiles, outfile):
    '''build aggregate table with cufflinks FPKM values.'''

    prefix = os.path.basename(outfile)
    prefix = prefix[:prefix.index("_")]

    headers = ",".join(
        [re.match("fpkm.dir/.*_(.*).cufflinks", x).groups()[0]
         for x in infiles])

    statement = '''
    python %(scriptsdir)s/combine_tables.py
        --columns=1
        --skip-titles
        --headers=%(headers)s
        --take=FPKM fpkm.dir/%(prefix)s_*.fpkm_tracking.gz
    | perl -p -e "s/tracking_id/transcript_id/"
    | gzip
    > %(outfile)s
    '''
    P.run()


@transform((mergeCufflinksGeneFPKM, mergeCufflinksIsoformFPKM),
           suffix(".tsv.gz"),
           ".load")
def loadCufflinksFPKM(infile, outfile):
    '''load fkpm data into table.'''

    P.load(infile, outfile,
           "--index=gene_id --index=transcript_id")


#########################################################################
#########################################################################
#########################################################################
CUFFDIFF_TARGETS_DE = [
    (
        (x, y, glob.glob("*.bam")),
        "%s_%s.tsv.gz" % (P.snip(x, ".tsv"),
                          P.snip(y, ".gtf.gz")))
    for x, y in itertools.product(glob.glob("design*.tsv"),
                                  glob.glob("*.gtf.gz"))
]


@follows(mkdir("cuffdiff.dir"), buildMaskGtf)
@files([(x, os.path.join("cuffdiff.dir", y))
        for x, y in CUFFDIFF_TARGETS_DE])
def runCuffdiff(infiles, outfile):
    '''perform differential expression analysis using cuffdiff.'''

    design_file = infiles[0]
    geneset_file = infiles[1]
    bamfiles = infiles[2]

    if PARAMS["cuffdiff_include_mask"]:
        mask_file = os.path.abspath("geneset_mask.gtf")
    else:
        mask_file = None

    options = PARAMS["cuffdiff_options"] + \
        " --library-type %s" % PARAMS["cufflinks_library_type"]

    Expression.runCuffdiff(bamfiles,
                           design_file,
                           geneset_file,
                           outfile,
                           threads=PARAMS.get("cuffdiff_threads", 4),
                           cuffdiff_options=options,
                           fdr=PARAMS["cuffdiff_fdr"],
                           mask_file=mask_file)

#########################################################################
#########################################################################
#########################################################################


@transform(runCuffdiff,
           suffix(".tsv.gz"),
           "_cuffdiff.load")
def loadCuffdiff(infile, outfile):
    '''load results from differential expression analysis and produce
    summary plots.

    Note: converts from ln(fold change) to log2 fold change.

    The cuffdiff output is parsed.

    Pairwise comparisons in which one gene is not expressed (fpkm <
    fpkm_silent) are set to status 'NOCALL'. These transcripts might
    nevertheless be significant.

    '''

    Expression.loadCuffdiff(infile, outfile)

#########################################################################
#########################################################################
#########################################################################


def buildExpressionStats(tables, method, outfile, outdir):
    '''build expression summary statistics.

    Creates also diagnostic plots in

    <exportdir>/<method> directory.
    '''
    dbhandle = sqlite3.connect(PARAMS["database"])

    def _split(tablename):
        # this would be much easier, if feature_counts/gene_counts/etc.
        # would not contain an underscore.
        try:
            design, geneset, counting_method = re.match(
                "([^_]+)_vs_([^_]+)_(.*)_%s" % method,
                tablename).groups()
        except AttributeError:
            try:
                design, geneset = re.match(
                    "([^_]+)_([^_]+)_%s" % method,
                    tablename).groups()
                counting_method = "na"
            except AttributeError:
                raise ValueError("can't parse tablename %s" % tablename)

        return design, geneset, counting_method

        # return re.match("([^_]+)_", tablename ).groups()[0]

    keys_status = "OK", "NOTEST", "FAIL", "NOCALL"

    outf = IOTools.openFile(outfile, "w")
    outf.write("\t".join(
        ("design",
         "geneset",
         "level",
         "treatment_name",
         "counting_method",
         "control_name",
         "tested",
         "\t".join(["status_%s" % x for x in keys_status]),
         "significant",
         "twofold")) + "\n")

    all_tables = set(Database.getTables(dbhandle))

    for level in CUFFDIFF_LEVELS:

        for tablename in tables:

            tablename_diff = "%s_%s_diff" % (tablename, level)
            tablename_levels = "%s_%s_diff" % (tablename, level)
            design, geneset, counting_method = _split(tablename_diff)
            if tablename_diff not in all_tables:
                continue

            def toDict(vals, l=2):
                return collections.defaultdict(
                    int,
                    [(tuple(x[:l]), x[l]) for x in vals])

            tested = toDict(
                Database.executewait(
                    dbhandle,
                    "SELECT treatment_name, control_name, "
                    "COUNT(*) FROM %(tablename_diff)s "
                    "GROUP BY treatment_name,control_name" % locals()
                    ).fetchall())
            status = toDict(Database.executewait(
                dbhandle,
                "SELECT treatment_name, control_name, status, "
                "COUNT(*) FROM %(tablename_diff)s "
                "GROUP BY treatment_name,control_name,status"
                % locals()).fetchall(), 3)
            signif = toDict(Database.executewait(
                dbhandle,
                "SELECT treatment_name, control_name, "
                "COUNT(*) FROM %(tablename_diff)s "
                "WHERE significant "
                "GROUP BY treatment_name,control_name" % locals()
                ).fetchall())

            fold2 = toDict(Database.executewait(
                dbhandle,
                "SELECT treatment_name, control_name, "
                "COUNT(*) FROM %(tablename_diff)s "
                "WHERE (l2fold >= 1 or l2fold <= -1) AND significant "
                "GROUP BY treatment_name,control_name,significant"
                % locals()).fetchall())

            for treatment_name, control_name in tested.keys():
                outf.write("\t".join(map(str, (
                    design,
                    geneset,
                    level,
                    counting_method,
                    treatment_name,
                    control_name,
                    tested[(treatment_name, control_name)],
                    "\t".join(
                        [str(status[(treatment_name, control_name, x)])
                         for x in keys_status]),
                    signif[(treatment_name, control_name)],
                    fold2[(treatment_name, control_name)]))) + "\n")

            ###########################################
            ###########################################
            ###########################################
            # plot length versus P-Value
            data = Database.executewait(
                dbhandle,
                "SELECT i.sum, pvalue "
                "FROM %(tablename_diff)s, "
                "%(geneset)s_geneinfo as i "
                "WHERE i.gene_id = test_id AND "
                "significant" % locals()).fetchall()

            # require at least 10 datapoints - otherwise smooth scatter fails
            if len(data) > 10:
                data = zip(*data)

                pngfile = "%(outdir)s/%(design)s_%(geneset)s_%(level)s_pvalue_vs_length.png" % locals()
                R.png(pngfile)
                R.smoothScatter(R.log10(ro.FloatVector(data[0])),
                                R.log10(ro.FloatVector(data[1])),
                                xlab='log10( length )',
                                ylab='log10( pvalue )',
                                log="x", pch=20, cex=.1)

                R['dev.off']()

    outf.close()

#########################################################################
#########################################################################
#########################################################################


@follows(mkdir(os.path.join(PARAMS["exportdir"], "cuffdiff")))
@transform(loadCuffdiff,
           suffix(".load"),
           ".plots")
def buildCuffdiffPlots(infile, outfile):
    '''create summaries of cufflinks results (including some diagnostic plots)

    Plots are created in the <exportdir>/cuffdiff directory.

    Plots are:

    <geneset>_<method>_<level>_<track1>_vs_<track2>_significance.png
        fold change against expression level
    '''
    ###########################################
    ###########################################
    # create diagnostic plots
    ###########################################
    outdir = os.path.join(PARAMS["exportdir"], "cuffdiff")

    dbhandle = sqlite3.connect(PARAMS["database"])

    prefix = P.snip(infile, ".load")

    geneset, method = prefix.split("_")

    for level in CUFFDIFF_LEVELS:
        tablename_diff = prefix + "_%s_diff" % level
        tablename_levels = prefix + "_%s_levels" % level

        # note that the ordering of EXPERIMENTS and the _diff table
        # needs to be the same as only one triangle is stored of the
        # pairwise results.  do not plot "undefined" lfold values
        # (where treatment_mean or control_mean = 0) do not plot lfold
        # values where the confidence bounds contain 0.
        for track1, track2 in itertools.combinations(EXPERIMENTS, 2):
            statement = """
            SELECT CASE WHEN d.treatment_mean < d.control_mean
            THEN d.treatment_mean
            ELSE d.control_mean END,
            d.l2fold, d.significant
            FROM %(tablename_diff)s AS d
            WHERE treatment_name = '%(track1)s' AND
            control_name = '%(track2)s' AND
            status = 'OK' AND
            treatment_mean > 0 AND
            control_mean > 0
            """ % locals()

            data = zip(*Database.executewait(dbhandle, statement))

            pngfile = "%(outdir)s/%(geneset)s_%(method)s_%(level)s_%(track1)s_vs_%(track2)s_significance.png" % locals()

            # ian: Bug fix: moved R.png to after data check so that no
            #     plot is started if there is no data this was leading
            #     to R falling over from too many open devices

            if len(data) == 0:
                E.warn("no plot for %s - %s -%s vs %s" %
                       (pngfile, level, track1, track2))
                continue

            R.png(pngfile)
            R.plot(ro.FloatVector(data[0]),
                   ro.FloatVector(data[1]),
                   xlab='min(FPKM)',
                   ylab='log2fold',
                   log="x", pch=20, cex=.1,
                   col=R.ifelse(ro.IntVector(data[2]), "red", "black"))

            R['dev.off']()

    P.touch(outfile)

#########################################################################
#########################################################################
#########################################################################


@follows(loadGeneSetGeneInformation)
@merge(loadCuffdiff,
       "cuffdiff_stats.tsv")
def buildCuffdiffStats(infiles, outfile):
    tablenames = [P.toTable(x) for x in infiles]
    outdir = os.path.dirname(infiles[0])
    buildExpressionStats(tablenames, "cuffdiff", outfile, outdir)

#########################################################################
#########################################################################
#########################################################################


@transform(buildCuffdiffStats,
           suffix(".tsv"),
           ".load")
def loadCuffdiffStats(infile, outfile):
    '''import cuffdiff results.'''
    P.load(infile, outfile)

#########################################################################
#########################################################################
#########################################################################


@merge(os.path.join(PARAMS["annotations_dir"],
                    PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
       "coding_exons.gtf.gz")
def buildCodingExons(infile, outfile):
    '''compile set of protein coding exons.

    This set is used for splice-site validation
    '''

    statement = '''
    zcat %(infile)s
    | awk '$2 == "protein_coding" && $3 == "CDS"'
    | perl -p -e "s/CDS/exon/"
    | python %(scriptsdir)s/gtf2gtf.py --merge-exons --log=%(outfile)s.log
    | gzip
    > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################


@transform("*.gtf.gz",
           suffix(".gtf.gz"),
           ".unionintersection.bed.gz")
def buildUnionIntersectionExons(infile, outfile):
    '''build union/intersection genes according to
    Bullard et al. (2010) BMC Bioinformatics.

    Builds a single-segment bed file.
    '''

    statement = '''
    gunzip < %(infile)s
    | python %(scriptsdir)s/gtf2gtf.py --intersect-transcripts
    --with-utr --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2gff.py --is-gtf
    --crop-unique  --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2bed.py --is-gtf --log=%(outfile)s.log
    | sort -k1,1 -k2,2n
    | gzip
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################


@transform("*.gtf.gz",
           suffix(".gtf.gz"),
           ".union.bed.gz")
def buildUnionExons(infile, outfile):
    '''build union genes.

    Exons across all transcripts of a gene are merged.
    They are then intersected between genes to remove any overlap.

    Builds a single-segment bed file.
    '''

    statement = '''
    gunzip < %(infile)s
    | python %(scriptsdir)s/gtf2gtf.py
         --merge-exons --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2gff.py
         --is-gtf --crop-unique  --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2bed.py
         --is-gtf --log=%(outfile)s.log
    | sort -k1,1 -k2,2n
    | gzip
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################


@follows(mkdir("gene_counts.dir"))
@files([(("%s.bam" % x.asFile(), "%s.gtf.gz" % y.asFile()),
         ("gene_counts.dir/%s_vs_%s.tsv.gz" % (x.asFile(), y.asFile())))
        for x, y in itertools.product(TRACKS, GENESETS)])
def buildGeneLevelReadCounts(infiles, outfile):
    '''compute read counts and coverage of exons with reads.
    '''

    bamfile, exons = infiles

    if BamTools.isPaired(bamfile):
        counter = 'readpair-counts'
    else:
        counter = 'read-counts'

    # ignore multi-mapping reads
    statement = '''
    zcat %(exons)s
    | python %(scriptsdir)s/gtf2table.py
          --reporter=genes
          --bam-file=%(bamfile)s
          --counter=length
          --prefix="exons_"
          --counter=%(counter)s
          --prefix=""
          --counter=read-coverage
          --prefix=coverage_
          --min-mapping-quality=%(counting_min_mapping_quality)i
          --multi-mapping=ignore
          --log=%(outfile)s.log
    | gzip
    > %(outfile)s
    '''

    P.run()


@transform(buildGeneLevelReadCounts,
           suffix(".tsv.gz"),
           "_gene_counts.load")
def loadGeneLevelReadCounts(infile, outfile):
    P.load(infile, outfile, options="--index=gene_id")

#########################################################################
#########################################################################
#########################################################################


@collate(buildGeneLevelReadCounts,
         regex("gene_counts.dir/(.+)_vs_(.+).tsv.gz"),
         r"gene_counts.dir/\2.gene_counts.tsv.gz")
def aggregateGeneLevelReadCounts(infiles, outfile):
    ''' build a matrix of counts with genes and tracks dimensions '''

    infiles = " ".join(infiles)
    # use anysense unique counts, needs to parameterized
    # for stranded/unstranded rnaseq data
    statement = '''python %(scriptsdir)s/combine_tables.py
    --columns=1
    --take=%(counting_type)s
    --use-file-prefix
    --regex-filename='(.+)_vs.+.tsv.gz'
    --log=%(outfile)s.log
    %(infiles)s
    | sed 's/geneid/gene_id/'
    | gzip > %(outfile)s '''

    P.run()


#########################################################################
#########################################################################
#########################################################################
@follows(mkdir("extension_counts.dir"))
@transform("*.bam",
           regex(r"(\S+).bam"),
           r"extension_counts.dir/\1.extension_counts.tsv.gz")
def buildGeneLevelReadExtension(infile, outfile):
    '''compute extension of cds.

    Known UTRs are counted as well.
    '''

    cds = os.path.join(PARAMS["annotations_dir"],
                       PARAMS_ANNOTATIONS["interface_geneset_cds_gtf"])

    territories = os.path.join(PARAMS["annotations_dir"],
                               PARAMS_ANNOTATIONS["interface_territories_gff"])

    utrs = os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_annotation_gff"])

    if "geneset_remove_contigs" in PARAMS:
        remove_contigs = '''| awk '$1 !~ /%s/' ''' % PARAMS[
            "geneset_remove_contigs"]
    else:
        remove_contigs = ""

    statement = '''
    zcat %(cds)s
    %(remove_contigs)s
    | python %(scriptsdir)s/gtf2table.py
          --reporter=genes
          --bam-file=%(infile)s
          --counter=position
          --counter=read-extension
          --min-mapping-quality=%(counting_min_mapping_quality)i
          --output-filename-pattern=%(outfile)s.%%s.tsv.gz
          --filename-gff=%(territories)s
          --filename-gff=%(utrs)s
          --log=%(outfile)s.log
    | gzip
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################


@follows(mkdir("transcript_counts.dir"))
@files([(("%s.bam" % x.asFile(), "%s.gtf.gz" % y.asFile()),
         ("transcript_counts.dir/%s_vs_%s.tsv.gz" % (x.asFile(), y.asFile())))
        for x, y in itertools.product(TRACKS, GENESETS)])
def buildTranscriptLevelReadCounts(infiles, outfile):
    '''count reads falling into transcripts of protein coding gene models.

    .. note::
       In paired-end data sets each mate will be counted. Thus
       the actual read counts are approximately twice the fragment
       counts.

    '''
    bamfile, geneset = infiles

    if BamTools.isPaired(bamfile):
        counter = 'readpair-counts'
    else:
        counter = 'read-counts'

    statement = '''
    zcat %(geneset)s
    | python %(scriptsdir)s/gtf2table.py
          --reporter=transcripts
          --bam-file=%(bamfile)s
          --counter=length
          --prefix="exons_"
          --counter=%(counter)s
          --prefix=""
          --counter=read-coverage
          --prefix=coverage_
          --min-mapping-quality=%(counting_min_mapping_quality)i
          --multi-mapping=ignore
          --log=%(outfile)s.log
    | gzip
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################


@transform(buildTranscriptLevelReadCounts,
           suffix(".tsv.gz"),
           ".load")
def loadTranscriptLevelReadCounts(infile, outfile):
    P.load(infile, outfile, options="--index=transcript_id")


#########################################################################
#########################################################################
#########################################################################
@follows(mkdir("feature_counts.dir"))
@files([(("%s.bam" % x.asFile(), "%s.gtf.gz" % y.asFile()),
         ("feature_counts.dir/%s_vs_%s.tsv.gz" % (x.asFile(), y.asFile())))
        for x, y in itertools.product(TRACKS, GENESETS)])
def buildFeatureCounts(infiles, outfile):
    '''counts reads falling into "features", which by default are genes.

    A read overlaps if at least one bp overlaps.

    Pairs and strandedness can be used to resolve reads falling into
    more than one feature. Reads that cannot be resolved to a single
    feature are ignored.

    '''
    bamfile, annotations = infiles
    PipelineRnaseq.runFeatureCounts(
        annotations,
        bamfile,
        outfile,
        nthreads=PARAMS['featurecounts_threads'],
        strand=PARAMS['featurecounts_strand'],
        options=PARAMS['featurecounts_options'])

#########################################################################
#########################################################################
#########################################################################


@collate(buildFeatureCounts,
         regex("feature_counts.dir/(.+)_vs_(.+).tsv.gz"),
         r"feature_counts.dir/\2.feature_counts.tsv.gz")
def aggregateFeatureCounts(infiles, outfile):
    ''' build a matrix of counts with genes and tracks dimensions.
    '''

    # Use column 7 as counts This is a possible source of bugs, the
    # column position has changed before.

    infiles = " ".join(infiles)
    statement = '''python %(scriptsdir)s/combine_tables.py
                                            --columns=1
                                            --take=7
                                            --use-file-prefix
                                            --regex-filename='(.+)_vs.+.tsv.gz'
                                            --log=%(outfile)s.log
                                             %(infiles)s
                  | sed 's/geneid/gene_id/'
                  | gzip > %(outfile)s '''

    P.run()


@transform(aggregateFeatureCounts,
           suffix(".tsv.gz"),
           ".load")
def loadFeatureCounts(infile, outfile):
    '''load individual feature counts into database'''
    P.load(infile, outfile, "--index=gene_id")


@merge(buildFeatureCounts,
       "feature_counts_summary.load")
def loadFeatureCountsSummary(infiles, outfile):
    '''load feature counts summary data into table.'''
    infiles = [P.snip(x, ".gz") + ".summary" for x in infiles]
    P.mergeAndLoad(infiles, outfile, options="--index=track")


@transform((aggregateGeneLevelReadCounts,
            aggregateFeatureCounts),
           suffix(".tsv.gz"),
           "_stats.tsv")
def summarizeCounts(infile, outfile):
    '''perform summarization of read counts'''

    prefix = P.snip(outfile, ".tsv")
    job_options = "-l mem_free=32G"
    statement = '''python %(scriptsdir)s/runExpression.py
              --method=summary
              --filename-tags=%(infile)s
              --output-filename-pattern=%(prefix)s_
              --log=%(outfile)s.log
              > %(outfile)s'''
    P.run()


@follows(mkdir("designs.dir"))
@product("design*.tsv",
         formatter(".*/(?P<PART1>.*).tsv$"),
         (aggregateGeneLevelReadCounts,
          aggregateFeatureCounts),
         formatter(".*/(?P<PART2>.*).tsv.gz$"),
         "designs.dir/{PART1[0][0]}_vs_{PART2[1][0]}_stats.tsv")
def summarizeCountsPerDesign(infiles, outfile):
    '''perform summarization of read counts within experiments.
    '''

    design_file, counts_file = infiles
    prefix = P.snip(outfile, ".tsv")
    statement = '''python %(scriptsdir)s/runExpression.py
              --method=summary
              --filename-design=%(design_file)s
              --filename-tags=%(counts_file)s
              --output-filename-pattern=%(prefix)s_
              --log=%(outfile)s.log
              > %(outfile)s'''
    P.run()

#########################################################################
#########################################################################
#########################################################################


@transform((summarizeCounts,
            summarizeCountsPerDesign),
           suffix("_stats.tsv"), "_stats.load")
def loadTagCountSummary(infile, outfile):
    '''load windows summary.'''
    P.load(infile, outfile)
    P.load(P.snip(infile, ".tsv") + "_correlation.tsv",
           P.snip(outfile, "_stats.load") + "_correlation.load",
           options="--first-column=track")


@follows(loadTagCountSummary,
         loadFeatureCountsSummary,
         aggregateGeneLevelReadCounts,
         aggregateFeatureCounts)
def counting():
    pass


@follows(mkdir("deseq.dir"), counting)
@product("design*.tsv",
         formatter("(.*).tsv$"),
         (aggregateGeneLevelReadCounts,
          aggregateFeatureCounts),
         formatter("(.*).tsv.gz$"),
         "deseq.dir/{basename[0][0]}_vs_{basename[1][0]}.gz")
def runDESeq(infiles, outfile):
    '''perform differential expression analysis using deseq.'''

    design_file, count_file = infiles

    track = P.snip(outfile, ".tsv.gz")

    statement = '''python %(scriptsdir)s/runExpression.py
    --method=deseq
    --filename-tags=%(count_file)s
    --filename-design=%(design_file)s
    --output-filename-pattern=%(track)s_
    --outfile=%(outfile)s
    --fdr=%(deseq_fdr)f
    --deseq-fit-type=%(deseq_fit_type)s
    --deseq-dispersion-method=%(deseq_dispersion_method)s
    --deseq-sharing-mode=%(deseq_sharing_mode)s
    --filter-min-counts-per-row=%(tags_filter_min_counts_per_row)i
    --filter-min-counts-per-sample=%(tags_filter_min_counts_per_sample)i
    --filter-percentile-rowsums=%(tags_filter_percentile_rowsums)i
    > %(outfile)s.log '''

    P.run()


@transform(runDESeq, suffix(".tsv.gz"), "_deseq.load")
def loadDESeq(infile, outfile):
    '''load differential expression results.
    '''
    # add gene level follow convention "<level>_diff"
    tablename = P.toTable(outfile) + "_gene_diff"
    statement = '''zcat %(infile)s
            | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --index=test_id
              --table=%(tablename)s
            > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################


@follows(loadGeneSetGeneInformation)
@merge(loadDESeq, "deseq_stats.tsv")
def buildDESeqStats(infiles, outfile):
    tablenames = [P.toTable(x) for x in infiles]
    outdir = os.path.dirname(infiles[0])
    buildExpressionStats(tablenames, "deseq", outfile, outdir)

#########################################################################
#########################################################################
#########################################################################


@transform(buildDESeqStats,
           suffix(".tsv"),
           ".load")
def loadDESeqStats(infile, outfile):
    P.load(infile, outfile)

#########################################################################
#########################################################################
#########################################################################


@follows(counting, mkdir("edger.dir"))
@product("design*.tsv",
         formatter("(.*).tsv$"),
         (aggregateGeneLevelReadCounts,
          aggregateFeatureCounts),
         formatter("(.*).tsv.gz$"),
         "edger.dir/{basename[0][0]}_vs_{basename[1][0]}.gz")
def runEdgeR(infiles, outfile):
    '''perform differential expression analysis using edger.'''

    design_file, count_file = infiles
    track = P.snip(outfile, ".tsv.gz")

    statement = '''python %(scriptsdir)s/runExpression.py
    --method=edger
    --filename-tags=%(count_file)s
    --filename-design=%(design_file)s
    --output-filename-pattern=%(track)s_
    --outfile=%(outfile)s
    --fdr=%(edger_fdr)f
    --filter-min-counts-per-row=%(tags_filter_min_counts_per_row)i
    --filter-min-counts-per-sample=%(tags_filter_min_counts_per_sample)i
    --filter-percentile-rowsums=%(tags_filter_percentile_rowsums)i
    > %(outfile)s.log '''

    P.run()

#########################################################################
#########################################################################
#########################################################################


@transform(runEdgeR, suffix(".tsv.gz"), "_edger.load")
def loadEdgeR(infile, outfile):
    '''load differential expression results.
    '''
    # add gene level follow convention "<level>_diff"
    tablename = P.toTable(outfile) + "_gene_diff"
    statement = '''zcat %(infile)s
    | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
    --allow-empty
    --index=test_id
    --table=%(tablename)s
    > %(outfile)s
    '''
    P.run()


@follows(loadGeneSetGeneInformation)
@merge(loadEdgeR, "edger_stats.tsv")
def buildEdgeRStats(infiles, outfile):
    tablenames = [P.toTable(x) for x in infiles]
    outdir = os.path.dirname(infiles[0])
    buildExpressionStats(tablenames, "edger", outfile, outdir)


@transform(buildEdgeRStats,
           suffix(".tsv"),
           ".load")
def loadEdgeRStats(infile, outfile):
    P.load(infile, outfile)


@follows(loadCufflinks,
         loadCufflinksFPKM,
         loadGeneLevelReadCounts)
def expression():
    pass

mapToTargets = {'cuffdiff': loadCuffdiffStats,
                'deseq': loadDESeqStats,
                'edger': loadEdgeRStats,
                }
TARGETS_DIFFEXPRESSION = [mapToTargets[x] for x in
                          P.asList(PARAMS["methods"])]


@follows(*TARGETS_DIFFEXPRESSION)
def diff_expression():
    pass


@follows(diff_expression)
@merge("*_stats.tsv", "de_stats.load")
def loadDEStats(infiles, outfile):
    '''load DE stats into table.'''
    P.concatenateAndLoad(infiles, outfile,
                         missing_value=0,
                         regex_filename="(.*)_stats.tsv")


@follows(mkdir("tagplots.dir"))
@product("design*.tsv",
         formatter("(.*).tsv$"),
         (aggregateGeneLevelReadCounts,
          aggregateFeatureCounts),
         formatter("(.*).tsv.gz$"),
         "tagplots.dir/{basename[0][0]}_vs_{basename[1][0]}")
def plotTagStats(infiles, outfile):
    '''perform differential expression analysis using deseq.'''

    design_file, counts_file = infiles

    statement = '''
    python %(scriptsdir)s/runExpression.py
    --log=%(outfile)s.log
    --filename-tags=%(counts_file)s
    --filename-design=%(design_file)s
    --method=plottagstats
    --output-filename-pattern=%(outfile)s
    > %(outfile)s
    '''
    P.run()

mapToQCTargets = {'cuffdiff': runCuffdiff,
                  'deseq': runDESeq,
                  'edger': runEdgeR,
                  }
QCTARGETS = [mapToQCTargets[x] for x in P.asList(PARAMS["methods"])]


@jobs_limit(1, "R")
@transform(QCTARGETS,
           suffix(".tsv.gz"),
           ".plots")
def plotDETagStats(infile, outfile):
    '''plot differential expression stats'''
    Expression.plotDETagStats(infile, outfile)
    P.touch(outfile)

###################################################################
###################################################################
###################################################################


@follows(plotTagStats,
         plotDETagStats,
         loadTagCountSummary,
         loadDEStats)
def qc():
    pass

###################################################################
###################################################################
###################################################################


@follows(expression, diff_expression, qc)
def full():
    pass

###################################################################
###################################################################
###################################################################


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''

    E.info("starting documentation build process from scratch")
    P.run_report(clean=True)

###################################################################
###################################################################
###################################################################


@follows(mkdir("report"))
def update_report():
    '''update report.'''

    E.info("updating documentation")
    P.run_report(clean=False)


@follows(update_report)
def publish():
    '''publish files.'''
    # publish web pages
    P.publish_report()


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
