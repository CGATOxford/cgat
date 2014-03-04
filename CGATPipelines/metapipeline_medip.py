"""
=====================
MeDIP pipeline - Meta
=====================

:Author: David Sims
:Release: $Id$
:Date: |today|
:Tags: Python

Meta pipeline for :doc:`pipeline_medip`. A meta pipeline analysis the output of multiple
runs of the same pipeline.

Methods
=======



Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

Input
-----

Source pipelines are added by linking them into the current directory. The pipeline expects the following
naming convention:

   medip_<name>

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

Pipeline output
===============

?

Example
=======

ToDo: make exome sequencing example


Code
====

"""

# load modules
from ruffus import *

import Experiment as E
import logging as L
import Database
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
import csv
import numpy
import sqlite3
import GTF
import IOTools
import IndexedFasta
import PipelineGeneset
import PipelineMapping
import Stats
import PipelineTracks
import PipelineMappingQC
import PipelineMedip
import Pipeline as P
import Expression

from rpy2.robjects import r as R
import rpy2.robjects as ro

#########################################################################
#########################################################################
#########################################################################
# load options from the config file
P.getParameters(["%s/pipeline.ini" % __file__[:-len(".py")],
                 "../pipeline.ini",
                 "pipeline.ini"])

PARAMS = P.PARAMS

PARAMS_ANNOTATIONS = P.peekParameters(PARAMS["annotations_dir"],
                                      "pipeline_annotations.py")

###################################################################
###################################################################
###################################################################
Sample = PipelineTracks.Sample
TRACKS = PipelineTracks.Tracks(Sample).loadFromDirectory(glob.glob("medip_*"),
                                                         "medip_(\S+)")


###################################################################
###################################################################
###################################################################
# if conf.py exists: execute to change the above assignmentsn
if os.path.exists("pipeline_conf.py"):
    L.info("reading additional configuration from pipeline_conf.py")
    execfile("pipeline_conf.py")

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

    statement = '''ATTACH DATABASE '%s' as medip_%s'''

    for track in TRACKS:
        cc = dbh.cursor()
        db = os.path.join("medip_%s" % track, "csvdb")
        cc.execute(statement % (db, track))
        cc.close()

    return dbh

###################################################################
###################################################################
###################################################################


@merge(None, "mapping.tsv")
def buildSummaryMapping(infiles, outfile):

    dbh = connect()
    cc = dbh.cursor()

    outf = IOTools.openFile(outfile, "w")

    table = "bam_stats"

    colnames = None
    for track in TRACKS:

        statement = """SELECT * 
                         FROM medip_%(track)s.%(table)s"""

        data = cc.execute(statement % locals()).fetchall()
        _colnames = [x[0] for x in cc.description]
        if not colnames:
            colnames = _colnames
            outf.write("\t".join(["metatrack"] + colnames,) + "\n")

        assert colnames == _colnames

        for row in data:
            outf.write("\t".join(map(str, (track,) + row)) + "\n")

    outf.close()

###################################################################
###################################################################
###################################################################


@merge(None, "called_dmrs.tsv")
def buildSummaryCalledDMRs(infiles, outfile):
    '''build summary of differentially methylated regions.'''

    dbh = connect()
    cc = dbh.cursor()

    outf = IOTools.openFile(outfile, "w")
    outf.write("metatrack\ttest\tntested\tnok\tnsignificant\tn2fold\n")

    for track in TRACKS:
        tables = [x[0] for x in cc.execute( """SELECT name FROM medip_%s.sqlite_master 
            WHERE type='table' and sql LIKE '%%control_mean%%' and sql LIKE '%%treatment_mean%%'""" % track
                                            ).fetchall()]

        for table in tables:

            statement = """SELECT 
                         COUNT(*) as ntested, 
                         SUM(CASE WHEN status='OK' THEN 1 ELSE 0 END) AS nok, 
                         SUM(CASE WHEN significant THEN 1 ELSE 0 END) AS nsignificant, 
                         SUM(CASE WHEN significant AND (l2fold < -1 OR l2fold > 1) THEN 1 ELSE 0 END) as n2fold 
                         FROM medip_%(track)s.%(table)s"""

            ntested, nok, nsignificant, n2fold = cc.execute(
                statement % locals()).fetchone()

            outf.write(
                "\t".join(map(str, (track, table, ntested, nok, nsignificant, n2fold))) + "\n")

    outf.close()

###################################################################
###################################################################
###################################################################


@merge(None, "cpg_coverage.tsv")
def buildSummaryCpGCoverage(infiles, outfile):
    '''build summary of differentially methylated regions.'''

    dbh = connect()
    cc = dbh.cursor()

    outf = IOTools.openFile(outfile, "w")
    outf.write("metatrack\ttrack\tcoverage\tncovered\tpcovered\n")

    for track in TRACKS:

        tables = [x[0] for x in cc.execute( """SELECT name FROM medip_%s.sqlite_master 
            WHERE type='table' and name LIKE '%%coveredpos%%' """ % track
                                            ).fetchall()]

        for table in tables:

            statement = """SELECT '%(track)s' as metatrack,
                         '%(table)s' as track,
                         coverage, ncovered, pcovered FROM medip_%(track)s.%(table)s"""

            for x in cc.execute(statement % locals()):
                outf.write("\t".join(map(str, x)) + "\n")

    outf.close()

###################################################################
###################################################################
###################################################################


@transform((buildSummaryCalledDMRs,
            buildSummaryMapping,
            buildSummaryCpGCoverage),
           suffix(".tsv"), ".load")
def loadSummary(infile, outfile):
    '''load all summary tables.'''
    P.load(infile, outfile)

###################################################################
###################################################################
###################################################################


@follows(loadSummary)
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

###################################################################
###################################################################
###################################################################


@follows(mkdir("%s/bamfiles" % PARAMS["web_dir"]),
         mkdir("%s/medips" % PARAMS["web_dir"]),
         )
def publish():
    '''publish files.'''

    # publish web pages
    P.publish_report()

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
