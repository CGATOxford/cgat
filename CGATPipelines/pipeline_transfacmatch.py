"""

========================
Transfac match pipeline
========================

:Author: Nick Ilott
:Release: $Id: pipeline_transfacmatch
:Date: |today|
:Tags: Python

The transfac match pipeline takes a several set of :term:`bed` or
:term:`gtf` formatted files and scans intervals for transcription
factor binding motifs.

It performs the following analyses:
   * transfac match motif analysis
      * If bed files are submitted then intervals must be specified
         with a unique name id
      * If gtf files are submitted then specification of interval
         locations should be provided
      * enrichment of motifs in a foreground set of ids vs.
        background set of ids
      * foreground and background sets must be specified

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.
The pipeline looks for a configuration file in several places:

   1. The default configuration in the :term:`code directory`.
   2. A shared configuration file :file:`../pipeline.ini`.
   3. A local configuration :file:`pipeline.ini`.

The order is as above. Thus, a local configuration setting will
override a shared configuration setting and a default configuration
setting.

Configuration files follow the ini format (see the python
`ConfigParser <http://docs.python.org/library/configparser.html>`
documentation).
The configuration file is organized by section and the variables
are documented within the file. In order to get a local configuration
file in the current directory, type::

    python <codedir>/pipeline_transfacmatch.py config

The sphinxreport report requires a :file:`conf.py` and
:file:`sphinxreport.ini` file (see :ref:`PipelineReporting`).
To start with, use the files supplied with the Example_ data.

Input
v-----

Intervals / GTF entries
++++++++++++++++++++++++

Input are :term:`bed`-formatted files of intervals or :term:`gtf`
formatted files. If intervals are submiited then they should be at
least bed4 formatted, i.e., each interval should be labelled
(uniquely).


Background and foreground sets
+++++++++++++++++++++++++++++++

Background and foreground sets are required to test for enrichments
of TF motifs. These are therefore subsets of the bed or gtf file.
It consists of a :term:`tsv` file:

+---------+
|   id    |
+---------+
|  id1    |
|  id2    |
+---------+

Background and foreground files must take the form
<name>.background.tsv <name>.foreground.tsv, repectively.

Motif library
++++++++++++++

The pipeline requires that a transfac matrix library be present.
Its location is specified in the pipeline.ini.


Motif profiles
+++++++++++++++

The pipeline requires a file that contains a selection of matrices
to search with defined cutoffs.
Transfac match will search only matrices that are specified in profile.
It is advised that all matrices in the transfac profile are used.

It will return only hits with scores that are
equal or higher than specified in profiles.


Requirements
------------

+--------------------+-----------+-----------------------------------+
|*Program*           |*Version*  |*Purpose*                          |
+--------------------+-----------+-----------------------------------+
|(transfac) match    |           |sequence scanning for TF motifs    |
+--------------------+-----------+-----------------------------------+

Pipline Output
==============

The results of the computation are all stored in an sqlite relational
database :file:`csvdb`.

Example
=======

Example data is available at ...
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_transfacmatch.tgz
   tar -xvzf pipeline_transfacmatch.tgz
   cd pipeline_transfacmatch.dir
   python <srcdir>/pipeline_transfacmatch.py make full

.. note::
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

"""


# load modules
from ruffus import *
import CGAT.Experiment as E
import logging as L
import CGAT.Database as Database
import CGAT.CSV as CSV
import numpy as np
import fnmatch
import sqlite3
import CGAT.FastaIterator as FastaIterator
import random
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import sys
import os
import re
import shutil
import itertools
import math
import glob
import time
import gzip
import collections
import random
import pandas
import pandas.io.sql as pdsql
from rpy2.robjects import r as R
from rpy2.robjects.packages import importr
from rpy2.robjects.numpy2ri import numpy2ri
import rpy2.robjects as robjects
import CGATPipelines.PipelineTracks as PipelineTracks
import CGATPipelines.PipelineTransfacMatch as PipelineTFM

###############################################################################
###############################################################################
###############################################################################
use_cluster = True

import CGAT.Pipeline as P
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])
PARAMS = P.PARAMS

###############################################################################
# helper functions mapping input files
###############################################################################

INPUT_FILE = None
# input intervals / gtf entries
# they need to be in a single file
INPUT_FORMATS = (r"*.gtf.gz", r"*.bed.gz")
REGEX_FORMATS = regex(r"(\S+).(gtf.gz|bed.gz)")
input_file = []
for x in INPUT_FORMATS:
    input_file += glob.glob(x)
    if input_file:
        assert len(input_file) <= 1, ("Multiple bed/gtf files in working "
                                      "directory please make sure all"
                                      " intervals are in one file")
        INPUT_FILE = input_file[0]
    else:
        pass


# foreground and background sets
INPUT_FOREGROUND = glob.glob("*.foreground.tsv")
INPUT_BACKGROUND = glob.glob("*.background.tsv")

TRACK_FOREGROUND = [P.snip(x, ".foreground.tsv") for x in INPUT_FOREGROUND]
TRACK_BACKGROUND = [P.snip(x, ".background.tsv") for x in INPUT_BACKGROUND]

###############################################################################
###############################################################################
###############################################################################


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

# P.toTable()
###############################################################################
###############################################################################
###############################################################################


def filenameToTablename(filename):
    '''
    converts filename containing "." to tablename where "." converted to "_"
    '''
    return filename.replace(".", "_")

# P.touch()
###############################################################################
###############################################################################
###############################################################################


def sentinelFile(filename):
    '''
    create empty file for updating purposes
    '''
    outf = open(filename, "w")
    outf.write("file created for ruffus update")
    outf.close()

###############################################################################
###############################################################################
# Section 1: Convert interval file to fasta file
###############################################################################
###############################################################################


@follows(mkdir("fasta.dir"))
@transform(INPUT_FILE, REGEX_FORMATS, r"fasta.dir/\1.fasta")
def buildIntervalsFasta(infile, outfile):
    '''
    build fasta file from intervals. Alternatively
    if a gtf file is specified this function will
    use parameters specified in the .ini file to
    use intervals upstream / downstream of tss
    '''

    # define upstream and downstream extensions
    upstream = PARAMS["intervals_extension_upstream"]
    downstream = PARAMS["intervals_extension_downstream"]

    assert len(str(upstream)), ("extension_upstream cannot be of %s type."
                                "If no extension is to be used specify 0" %
                                type(upstream))
    assert len(str(downstream)), ("downstream extension cannot be of %s type."
                                  "If no extension is to be used specify 0" %
                                  type(downstream))

    # if input is gtf then convert to bed
    # with intervals defined by .ini file
    temp = P.getTempFile("/ifs/scratch")
    if infile.endswith(".gtf.gz"):
        # the resulting temporary file will not be zipped
        concatenate = "cat"
        for gene in GTF.merged_gene_iterator(
                GTF.iterator(IOTools.openFile(infile))):
            if gene.strand == "+":
                temp.write("%s\t%s\t%s\t%s\t%s\t%s\n" %
                           (gene.contig,
                            str(gene.start - upstream),
                            str(gene.start + downstream),
                            gene.gene_id,
                            ".",
                            gene.strand))
            elif gene.strand == "-":
                temp.write("%s\t%s\t%s\t%s\t%s\t%s\n" %
                           (gene.contig,
                            str(gene.end - downstream),
                            str(gene.end + upstream),
                            gene.gene_id,
                            ".",
                            gene.strand))
        temp.close()
        inf = temp.name
    else:
        inf = infile
        concatenate = "zcat"

    to_cluster = True

    # define statement
    # option to specify strand in config file.
    statement = ("%(concatenate)s %(inf)s |"
                 " python %(scriptsdir)s/bed2fasta.py"
                 "  --genome=%(genomedir)s/%(genome)s")

    if PARAMS["intervals_stranded"]:
        statement += ("  --use-strand --log=%(outfile)s.log > %(outfile)s")
    else:
        statement += ("  --log=%(outfile)s.log > %(outfile)s")

    P.run()

    if infile.endswith(".gtf.gz"):
        os.remove(inf)

###########
# target
###########


@follows(buildIntervalsFasta)
def Fasta():
    pass

###############################################################################
###############################################################################
# Section 2: Select background set with equal GC to foreground set
###############################################################################
###############################################################################


@follows(mkdir("GC_content.dir"), buildIntervalsFasta)
@transform(INPUT_BACKGROUND + INPUT_FOREGROUND,
           regex(r"(\S+).tsv"),
           add_inputs(buildIntervalsFasta),
           r"GC_content.dir/\1.gc.tsv")
def calculateGCContent(infiles, outfile):
    '''
    calculate the GC content across foreground and background sets
    '''
    PipelineTFM.calculateSequenceComposition(infiles[0],
                                             infiles[1],
                                             outfile,
                                             PARAMS["genesets_header"])

###############################################################################
###############################################################################
###############################################################################


@transform(calculateGCContent, suffix(".tsv"), ".load")
def loadGCContent(infile, outfile):
    '''
    load the results the GC content for each background
    and foreground
    '''
    tablename = filenameToTablename(P.snip(os.path.basename(infile), ".tsv"))
    statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s
                   --log=%(outfile)s.log
                   --add-index=id
                   %(csv2db_options)s
                   < %(infile)s > %(outfile)s'''
    P.run()

###############################################################################
###############################################################################
###############################################################################
if PARAMS["background_match"]:
    @collate(loadGCContent, regex("GC_content.dir/"
                                  "(.+)\.(?:background|foreground)\.gc\.load"),
             r"GC_content.dir/\1.matched.background.tsv")
    def matchBackgroundForSequenceComposition(infiles, outfile):
        '''
        take the background set and subset it for intervals with the same
        sequence composition distribution that is the same as the foreground
        set (for the composition statistic specified in config file).
        This requires that the background set is sufficiently large.
        '''
        # get gene set name
        track = re.match("GC_content.dir/"
                         "(.+)\.(?:background|foreground)\.gc\.load",
                         infiles[0]).groups()[0]

        # get list of foreground genes
        input_background = "%s.background.tsv" % track

        # get list of backround genes
        input_foreground = "%s.foreground.tsv" % track

        # get name of fasta file containing intervals
        fasta_file = os.path.basename(INPUT_FILE)[:-len(".gtf.gz")]
        fasta_file = os.path.join("fasta.dir", fasta_file)

        PipelineTFM.matchBgSequenceComposition(infiles,
                                               input_background,
                                               input_foreground,
                                               fasta_file,
                                               outfile,
                                               PARAMS["database"],
                                               PARAMS["genesets_header"],
                                               PARAMS["background_match_stat"],
                                               PARAMS["sig_testing_method"])

    ###############################################
    ###############################################
    ###############################################
    @transform(matchBackgroundForSequenceComposition,
               regex(r"(\S+).tsv"),
               add_inputs(buildIntervalsFasta),
               r"\1.gc.tsv")
    def calculateMatchedGCComposition(infiles, outfile):
        '''
        calculate the GC content for the CpG matched data
        Should be the same as the CpG content of the foreground
        set
        '''
        PipelineTFM.calculateSequenceComposition(infiles[0],
                                                 infiles[1],
                                                 outfile)

    ###############################################
    ###############################################
    ###############################################
    @transform(calculateMatchedGCComposition, suffix(".tsv"), ".load")
    def loadMatchedGCComposition(infile, outfile):
        '''
        load the CpG compostion of matched background set
        '''

        to_cluster = True
        tablename = filenameToTablename(
            P.snip(os.path.basename(outfile), ".load"))
        statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s
                       --log=%(outfile)s.log
                       %(csv2db_options)s
                       < %(infile)s > %(outfile)s'''
        P.run()

###############
# target
##############
if PARAMS["background_match"]:
    @follows(loadMatchedGCComposition)
    def GC():
        pass
else:
    @follows(loadGCContent)
    def GC():
        pass

###############################################################################
###############################################################################
# Section 3: Run TransFac Match
###############################################################################
###############################################################################


@follows(mkdir("match.dir"))
@merge([buildIntervalsFasta,
        PARAMS["transfac_matrix"],
        PARAMS["transfac_profile"]],
       "match.dir/match.result")
def runMatch(infiles, outfile):
    '''
    run transfac(R) match
    '''
    to_cluster = True
    seq = infiles[0]
    mxlib = PARAMS["transfac_matrix"]
    mxprf = PARAMS["transfac_profile"]
    match_executable = "/ifs/data/biobase/transfac/match/bin/match_linux64"
    out = outfile

    statement = '''%(match_executable)s %(mxlib)s %(seq)s %(out)s %(mxprf)s'''
    P.run()

###############################################################################
###############################################################################
###############################################################################


@transform(runMatch, suffix(".result"), ".load")
def loadMatchResults(infile, outfile):
    '''
    load the results of the match analysis into sqlite database
    '''
    temp = P.getTempFile("./match.dir")
    temp.write("seq_id\tmatrix_id\tposition\tstrand\t"
               "core_score\tmatrix_score\tsequence\n")
    for details in PipelineTFM.match_iterator(infile):
        temp.write("\t".join(map(str, [details.seq_id,
                                       details.matrix_id,
                                       details.position,
                                       details.strand,
                                       details.core_score,
                                       details.matrix_score,
                                       details.sequence])) + "\n")
    temp.close()

    to_cluster = True
    job_options = "-l mem_free=64G"

    inf = temp.name
    tablename = filenameToTablename(os.path.basename(infile))
    statement = ("python %(scriptsdir)s/csv2db.py"
                 "  -t %(tablename)s"
                 "  --log=%(outfile)s.log"
                 "  --add-index=seq_id"
                 "  %(csv2db_options)s"
                 " < %(inf)s > %(outfile)s")
    P.run()
    os.unlink(temp.name)

###############################################################################
###############################################################################
###############################################################################


@transform(loadMatchResults, suffix(".load"), ".metrics")
def buildMatchMetrics(infile, outfile):
    '''
    match outputs transcription factors that are found in the supplied
    sequences. We are interested in the following metrics:

       * No. unique transcription factors found per sequence

       * Maximal number of TF motifs found per sequence

    '''
    tablename = filenameToTablename(
        os.path.basename(P.snip(infile, ".load"))) + "_result"
    PipelineTFM.frequencyMetrics(PARAMS["database"], tablename, outfile)

###############################################################################
###############################################################################
###############################################################################


@transform(buildMatchMetrics, suffix(""), ".load")
def loadMatchMetrics(infile, outfile):
    '''
    load match metrics
    '''
    tablename = filenameToTablename(os.path.basename(infile))

    to_cluster = True

    job_options = " -l mem_free=16G"

    statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s
                   --log=%(outfile)s.log
                   --add-index=seq_id
                   %(csv2db_options)s
                   < %(infile)s > %(outfile)s'''
    P.run()

###########
# target
###########


@follows(loadMatchMetrics)
def Match():
    pass


###############################################################################
###############################################################################
# Section 4: For each transcription factor test for significant enrichment
#            between foreground and background gene sets.
###############################################################################
###############################################################################
if PARAMS['sig_testing_method'] == "fisher":

    if PARAMS["background_match"]:
        @follows(loadMatchResults,
                 loadMatchedGCComposition,
                 mkdir("match_test.dir"))
        #@jobs_limit(1, "FisherTest")
        @collate([matchBackgroundForSequenceComposition, calculateGCContent],
                 regex(".+/(.+)\.(?:foreground.gc|matched.background)\.tsv"),
                 r"match_test.dir/\1.matched.significance")
        def estimateEnrichmentOfTFBS(infiles, outfile):
            '''
            Estimate the significance of transcription factors that are
            associated with a foreground set of intervals vs a background
            set matched for sequence composition.
            '''
            E.info("Running Fisher's exact test for TF enrichment between %s" %
                   " & ".join([os.path.basename(x) for x in infiles]))

            # required files
            match_table = "match_result"

            # we don't know which order the foreground and background
            # will come in
            background = [infile for infile in infiles if
                          re.search("background", infile)][0]
            foreground = ["%s.foreground.tsv" %
                          re.match(".+/(.+)\.foreground\.gc\.tsv",
                                   infile).groups()[0]
                          for infile in infiles if re.search("foreground",
                                                             infile)][0]

            # run significance testing
            # MM: added in directionality into FET - might only be looking for
            # enrichment OR depletion so don't want to hammer those p-value
            # too hard
            # MM: 23/12/14 - refactor to run on cluster

            job_options = "-l mem_free=1G"

            statement = '''
            python %(scriptsdir)s/tfbs2enrichment.py
            --foreground=%(foreground)s
            --background=%(background)s
            --database=%(database)s
            --log=%(outfile)s.log
            --match-table=%(match_table)s
            --outfile=%(outfile)s
            --geneset-header=%(genesets_header)s
            --direction=%(sig_testing_direction)s
            '''

            P.run()

            E.info("Completed Fisher's exact test for "
                   "TF enrichment between %s" %
                   " & ".join([os.path.basename(x) for x in infiles]))

    else:
        @follows(loadMatchResults, mkdir("match_test.dir"))
        @jobs_limit(1, "FisherTest")
        @collate(INPUT_BACKGROUND + INPUT_FOREGROUND,
                 regex("(.+)\.(?:foreground|background)\.tsv"),
                 r"match_test.dir/\1.significance")
        def estimateEnrichmentOfTFBS(infiles, outfile):
            '''
            Estimate the significance of trnascription factors that are
            associated with a foreground set of intervals vs a background set
            '''
            E.info("Running Fisher's exact test for TF enrichment between %s" %
                   " & ".join([os.path.basename(x) for x in infiles]))

            # required files
            match_table = "match_result"

            # we don't know which order the foreground and backgorund
            # will come in
            background = [infile for infile in infiles if
                          re.search("background", infile)][0]
            foreground = [infile for infile in infiles if
                          re.search("foreground", infile)][0]

            # run significance testing

            PipelineTFM.testSignificanceOfMatrices(background,
                                                   foreground,
                                                   PARAMS["database"],
                                                   match_table,
                                                   outfile)

            E.info("Completed Fisher's exact test for "
                   "TF enrichment between %s" %
                   " & ".join([os.path.basename(x) for x in infiles]))

elif PARAMS['sig_testing_method'] == "permutation":
    @follows(loadMatchMetrics,
             mkdir("match_test.dir"))
    @collate([matchBackgroundForSequenceComposition, calculateGCContent],
             regex(".+/(.+)\.(?:foreground.gc|background.gc)\.tsv"),
             r"match_test.dir/\1.matched.significance")
    def estimateEnrichmentOfTFBS(infiles, outfile):
        '''
        Test for enrichment of TFBS within a gene set by permutation.
        '''

        E.info("Running permutation testing for TFBS enrichment between %s" %
               " & ".join([os.path.basename(x) for x in infiles]))

        dbh = sqlite3.connect(PARAMS['database'])
        # table from sql db
        match_table = "match_result"
        tfbs_state = '''SELECT matrix_id, seq_id FROM %s;''' % match_table
        tfbs_table = pdsql.read_sql(sql=tfbs_state,
                                    con=dbh,
                                    index_col='matrix_id')

        # get foreground and background gene files
        # setup gc content dataframes

        background = [inf for inf in infiles if re.search("background",
                                                          inf)][0]
        foreground = [inf for inf in infiles if re.search("foreground",
                                                          inf)][0]

        back_gc = pandas.read_table(background,
                                    sep="\t",
                                    index_col=0,
                                    header=0)
        bg_gene_id = [x.split(" ")[0] for x in back_gc.index.tolist()]
        back_gc['gene_id'] = bg_gene_id
        back_gc.index = bg_gene_id

        fore_gc = pandas.read_table(foreground,
                                    sep="\t",
                                    index_col=0,
                                    header=0)
        fg_gene_id = [x.split(" ")[0] for x in fore_gc.index.tolist()]
        fore_gc['gene_id'] = fg_gene_id
        fore_gc.index = fg_gene_id

        # run permutation significance testing
        # check if there are sufficient genes in the background
        # to do a permutation test.  If not, do Fishers' exact.
        # if yes, but less than specific number of permutations,
        # limit to this number.

        perms = int(PARAMS['sig_testing_nperms'])
        poss_perms = PipelineTFM.nCr(n=len(back_gc.index),
                                     r=len(fore_gc.index))
        if poss_perms < 1000:
            E.warn("Insufficient background genes to perform"
                   "permutations.  Please use Fisher's Exact test")
            raise ValueError("Insufficient background size. "
                             "Cannot use permutation test")
        elif poss_perms > 1000 and poss_perms < perms:
            E.info("Maximum possible permutations with this background"
                   " set is %i.  Running %i permutations only" % (poss_perms,
                                                                  poss_perms))
            out_dict = PipelineTFM.permuteTFBSEnrich(tfbs_table=tfbs_table,
                                                     fg_gc=fore_gc,
                                                     bg_gc=back_gc,
                                                     nPerms=poss_perms,
                                                     bg_stat=PARAMS[""])
        else:
            bg_stat = PARAMS["background_match_stat"]
            out_dict = PipelineTFM.permuteTFBSEnrich(tfbs_table=tfbs_table,
                                                     fg_gc=fore_gc,
                                                     bg_gc=back_gc,
                                                     nPerms=perms,
                                                     bg_stat=bg_stat)

        out_frame = pandas.DataFrame(out_dict).T
        pyadjust = R['p.adjust']
        pvs = robjects.FloatVector([p for p in out_frame['pvalue']])
        out_frame['qvalue'] = pyadjust(pvs)

        out_frame.to_csv(outfile, sep="\t", index_label='matrix_id')
###############################################################################
###############################################################################
###############################################################################


@transform(estimateEnrichmentOfTFBS, suffix(".significance"), ".load")
def loadEnrichmentOfTFBS(infile, outfile):
    '''
    load the results of the enrichment
    '''

    to_cluster = True
    tablename = filenameToTablename(os.path.basename(infile))
    statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s
                  --log=%(outfile)s.log
                  --add-index=matrix_id
                  %(csv2db_options)s
                  < %(infile)s > %(outfile)s'''
    P.run()


@collate(loadEnrichmentOfTFBS,
         regex("(.+)/(.+)\.load"),
         r"\1/all_genesets.significance.tsv")
def collateEnrichmentOfTFBS(infiles, outfile):
    '''
    Concatenate the enrichment scores for all gene sets into a single table
    and recalculate qvalues.
    '''

    def _fetch(table_name):
        table_name = table_name.replace("-", "_")
        # retrieve table
        dbh = sqlite3.connect(PARAMS["database"])
        cc = dbh.cursor()
        data = cc.execute("SELECT * FROM %s" % table_name).fetchall()
        cc.close()

        # return pandas dataframe
        headers = [d[0] for d in cc.description]
        return pandas.DataFrame.from_records(data, columns=headers)

    first = True
    for infile in infiles:
        table_name = filenameToTablename(P.snip(os.path.basename(infile)))
        table_name = table_name + "_significance"
        geneset_id = table_name.split("_")[0]
        if first:
            # set up first dataframe
            first = False
            df = _fetch(table_name)
            # removing qvalue as it will be recalculated
            try:
                df.drop("qvalue", axis=1, inplace=True)
            except ValueError:
                # qvalue not contained in df
                pass
            # adding column for geneset_id
            geneset_id = [geneset_id, ]*len(df.index)
            df["geneset_id"] = geneset_id
            continue

        # append successive dataframes
        df_n = _fetch(table_name)
        df_n.drop("qvalue", axis=1, inplace=True)
        geneset_id = [geneset_id, ]*len(df_n.index)
        df_n["geneset_id"] = geneset_id
        df = df.append(df_n)

    # get globally adjusted pvalues
    p_vals = df["pvalue"].tolist()
    padjpy = R["p.adjust"]
    # returns a python FloatVector object
    q_vals = padjpy(p_vals)
    df["qvalue"] = [x for x in q_vals]

    df.to_csv(outfile, sep="\t", index=False)


@transform(collateEnrichmentOfTFBS,
           suffix(".significance.tsv"),
           ".load")
def loadCollatedEnrichmentOfTFBS(infile, outfile):
    tablename = P.snip(os.path.basename(infile), ".significance.tsv")
    statement = ("python %(scriptsdir)s/csv2db.py"
                 "  -t %(tablename)s"
                 "  --log=%(outfile)s.log"
                 "  --add-index=matrix_id"
                 "  %(csv2db_options)s"
                 " < %(infile)s > %(outfile)s")
    P.run()


##############
# target
##############


@posttask(touch_file("complete.flag"))
@follows(loadCollatedEnrichmentOfTFBS)
def Significance():
    pass


@follows(Fasta, Match, Significance, GC)
def full():
    pass

###############################################################################
###############################################################################
###############################################################################
if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
