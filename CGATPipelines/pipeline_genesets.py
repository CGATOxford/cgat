"""===========================
Geneset analysis
===========================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Overview
========

This pipeline performs gene set analysis of one or more
genesets.

Input data are two collections of files, genelists and pathways.

Genelists are tabular data with a gene for each row and associated
attributes in additional columns such as expression level, probability
of being called differentially expressed, etc.

Pathways are tabular data linking genes to pathways that they exist
in.


Generally, it performs the following tasks:

1. The pipeline merges separately prepared gene lists into
   a single gene list matrix. There is a continuous scale
   version (P-Values, expression values, ...) and a thresholded
   version (0 and 1 for genelist membership).

2. The pipeline builds a matrix of gene list annotations to
   test against. To this end, it collects:

   ENSEMBL GO annotations
   KEGG Pathways
   User supplied pathways
   GSEA database signatures

3. The pipeline performs various gene set enrichment analyses.
   These are:

   1. Hypergeometric GO analysis
   2. Gene set enrichment analysis

4. The pipeline creates various QC metrics. To this end it looks
   for biases in any of the gene lists supplied. Biases the pipeline
   looks at are:

   1. Gene length
   2. Nucleotide composition
   3. Gene intron/exon structure
   4. User supplied table with biases.


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

Optional inputs
+++++++++++++++

Requirements
------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

+----------+-----------+---------------------------+
|*Program* |*Version*  |*Purpose*                  |
+----------+-----------+---------------------------+
|          |           |                           |
+----------+-----------+---------------------------+

Pipeline output
===============

The major output is in the database file :file:`csvdb`.

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *

import sys
import os
import sqlite3
import pandas

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Database as Database
import CGAT.SetTools as SetTools
import CGATPipelines.PipelineGO as PipelineGO

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# load options from the config file
import CGAT.Pipeline as P
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

PARAMS.update(P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    on_error_raise=__name__ == "__main__",
    prefix="annotations_",
    update_interface=True))


# Update the PARAMS dictionary in any PipelineModules
# e.g.:
# import CGATPipelines.PipelineGeneset as PipelineGeneset
# PipelineGeneset.PARAMS = PARAMS


def connect():
    '''connect to database.

    Use this method to connect to additional databases.

    Returns a database connection.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh


@transform('genelists.dir/*.tsv.gz',
           suffix(".tsv.gz"),
           ".load")
def loadGeneLists(infile, outfile):
    '''load gene list data into database.'''
    P.load(infile, outfile,
           tablename="genelist_%s" % P.toTable(outfile))


@merge('genelists.dir/*.tsv.gz', 'genelists.tsv.gz')
def buildGeneListMatrix(infiles, outfile):
    '''build a gene list matrix for simple pathway analysis
    based on hypergeometric test.

    A gene list is derived from a gene set by
    applying thresholds to the input data set. The
    thresholds are defined in the configuration file.
    '''

    genesets = []
    backgrounds = []
    headers = []
    for infile in infiles:
        genelist = pandas.read_csv(
            IOTools.openFile(infile),
            index_col=0,
            sep='\t')

        track = P.snip(os.path.basename(infile), ".tsv.gz")
        headers.append(track)

        field = PARAMS[P.matchParameter("%s_foreground_field" % track)]
        min_threshold = PARAMS[P.matchParameter(
            "%s_foreground_min_threshold" % track)]
        max_threshold = PARAMS[P.matchParameter(
            "%s_foreground_max_threshold" % track)]
        genesets.append(set(genelist[
            (genelist[field] >= min_threshold) &
            (genelist[field] <= max_threshold)].index))

        E.info('%s: foreground: %f <= %s <= %f' % (track,
                                                   min_threshold,
                                                   field,
                                                   max_threshold))

        field = PARAMS[P.matchParameter("%s_background_field" % track)]
        min_threshold = PARAMS[P.matchParameter(
            "%s_background_min_threshold" % track)]
        max_threshold = PARAMS[P.matchParameter(
            "%s_background_max_threshold" % track)]

        E.info('%s: background: %f <= %s <= %f' % (track,
                                                   min_threshold,
                                                   field,
                                                   max_threshold))
        backgrounds.append(set(genelist[
            (genelist[field] >= min_threshold) &
            (genelist[field] <= max_threshold)].index))

        E.info("%s: fg=%i, bg=%i" % (track,
                                     len(genesets[-1]),
                                     len(backgrounds[-1])))

    E.info("writing gene list matrix")
    with IOTools.openFile(outfile, "w") as outf:
        SetTools.writeSets(outf, genesets, labels=headers)
    with IOTools.openFile(outfile + ".bg.tsv.gz", "w") as outf:
        SetTools.writeSets(outf, backgrounds, labels=headers)

    E.info("writing intersection/union matrix")
    # build set intersection matrix
    matrix = SetTools.unionIntersectionMatrix(genesets)
    with IOTools.openFile(outfile + ".matrix.gz", "w") as outf:
        IOTools.writeMatrix(outf, matrix, headers, headers)
    matrix = SetTools.unionIntersectionMatrix(backgrounds)
    with IOTools.openFile(outfile + ".bg.matrix.gz", "w") as outf:
        IOTools.writeMatrix(outf, matrix, headers, headers)


@transform(buildGeneListMatrix,
           suffix(".tsv.gz"),
           ".load")
def loadGeneListMatrix(infile, outfile):
    '''load fgene list matrix into table.'''
    track = P.snip(infile, ".tsv.gz")
    P.load(infile, outfile, tablename="%s_foreground" % track)
    P.load(infile + ".bg.tsv.gz", outfile, tablename="%s_background" % track)


@transform('pathways.dir/*.tsv.gz',
           regex('.*/(.*).tsv.gz'),
           r"pathways_\1.load")
def loadPathways(infile, outfile):
    '''load pathway information into database.'''
    P.load(infile, outfile, "--index=gene_id --index=go_id")


@follows(mkdir('hypergeometric.dir'))
@transform('pathways.dir/*.tsv.gz',
           regex('.*/(.*).tsv.gz'),
           add_inputs(buildGeneListMatrix),
           r'hypergeometric.dir/\1.tsv')
def runHypergeometricAnalysis(infiles, outfile):
    '''run pathway analysis on pathway files in
    the directory pathways.dir.
    '''
    infile_pathways, infile_genelist = infiles
    infile_background = infile_genelist + ".bg.tsv.gz"

    # TODO:
    # gene annotations
    # category annotations
    #
    # os.path.join(
    #        PARAMS["annotations_dir"],
    #        PARAMS_ANNOTATIONS["interface_go_obo"]),

    PipelineGO.runGOFromFiles(
        outfile=outfile,
        outdir=outfile + ".dir",
        fg_file=infile_genelist,
        bg_file=infile_background,
        go_file=infile_pathways,
        ontology_file=None,
        minimum_counts=PARAMS["hypergeometric_minimum_counts"],
        pairs=False,
        gene2name=None)


def computePathwayBiases(infile, outfile):
    pass


@transform(runHypergeometricAnalysis,
           suffix(".tsv"),
           r"\1.load")
def loadHypergeometricAnalysis(infile, outfile):
    '''load GO results.'''

    track = P.toTable(outfile)
    tablename = 'hypergeometric_%s_summary' % track
    P.load(infile, outfile, tablename=tablename)

    dbh = connect()
    ontologies = [x[0] for x in Database.executewait(
        dbh,
        '''SELECT DISTINCT ontology FROM %s''' % tablename).fetchall()]

    genelists = [x[0] for x in Database.executewait(
        dbh,
        '''SELECT DISTINCT genelist FROM %s''' % tablename).fetchall()]

    # output files from runGO.py
    sections = ('results', 'parameters', 'withgenes')

    for section in sections:
        tablename = 'hypergeometric_%s_%s' % (track, section)
        statement = '''
        python %(scriptsdir)s/combine_tables.py
        --cat=track
        --regex-filename="hypergeometric.dir/%(track)s.tsv.dir/(\S+).%(section)s"
        hypergeometric.dir/%(track)s.tsv.dir/*.%(section)s
        | python %(scriptsdir)s/csv2db.py
        %(csv2db_options)s
        --table=%(tablename)s
        >> %(outfile)s'''
        P.run()

    for ontology in ontologies:

        fn = os.path.join(infile + ".dir", "all_alldesc.%s.l2fold" % ontology)

        if not os.path.exists(fn):
            E.warn("file %s does not exist" % fn)
            continue

        P.load(fn,
               outfile,
               tablename='hypergeometric_%s_%s_l2fold' % (track, ontology),
               options='--allow-empty')

        fn = os.path.join(
            infile + ".dir", "all_alldesc.%s.l10pvalue" % ontology)

        P.load(fn,
               outfile,
               tablename='hypergeometric_%s_%s_l10pvalue' % (track, ontology),
               options='--allow-empty')

        fn = os.path.join(
            infile + ".dir", "all_alldesc.%s.l10qvalue" % ontology)

        P.load(fn,
               outfile,
               tablename='hypergeometric_%s_%s_l10qvalue' % (track, ontology),
               options='--allow-empty')


@merge(runHypergeometricAnalysis,
       "hypergeometric_summary.load")
def loadHypergeometricResultsSummary(infiles, outfile):
    '''load GO summary results.'''
    infiles = glob.glob("hypergeometric.dir/*/*.parameters")
    P.mergeAndLoad(infiles, outfile)


@collate("hypergeometric.dir/go.tsv.dir/*.results",
         regex(r"hypergeometric.dir/go.tsv.dir/(.*)\.(.*).results"),
         r"hypergeometric.go.dir/go.tsv.dir/\1.revigo")
def plotGOResults(infiles, outfile):
    '''.'''

    infiles = " ".join(infiles)

    track = P.snip(outfile, ".revigo")

    statement = '''
    cat %(infiles)s
    | python %(scriptsdir)s/revigo.py
      --filename-go=%(annotations_filename_go)s
      --output-filename-pattern=%(track)s.%%s
      --ontology=all
      --max-similarity=0.5
      --reverse-palette
      --force
      -v 2
    > %(outfile)s
    '''

    P.run()


@follows(loadPathways,
         loadGeneLists,
         loadGeneListMatrix,
         loadHypergeometricAnalysis)
def full():
    pass


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.'''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish_report():
    '''publish report.'''

    E.info("publishing report")
    P.publish_report()

if __name__ == "__main__":

    # P.checkFiles( ("genome.fasta", "genome.idx" ) )
    sys.exit(P.main(sys.argv))
