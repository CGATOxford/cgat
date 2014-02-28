"""=================================================
pipeline_testing - automated testing of pipelines
=================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

This pipeline executes other pipelines for testing purposes.

Overview
========

This pipeline implements automated testing of CGAT pipelines. This
pipeline is run regularly on the code in the master repository in
order to check if all pipelines are functional.

This pipeline can also be run from within a cloned repository before
checking in changes. This is highly recommended, as changes in shared
scripts might break existing workflows.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

In order to run all tests, simply enter an empty directory and type::

   python <srcdir>/pipeline_testing.py make full
   python <srcdir>/pipeline_testing.py make build_report

The first command will run all pipelines while the second will build a summary
report.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineReporting`). To start with, use the files supplied with the
Example_ data.

Input
-----

As input, the pipeline requires a directory with test cases. The
directory is specified in the ``data_dir`` configuration
variable. Each data directory needs to start with the name of the
pipeline that is to be tested. Valid names for test cases are thus
``pipeline_annotations``, ``pipeline_annotations_lizard``,
``pipeline_rnaseq_illumina``, ``pipeline_rnaseq_solid``, etc.

A test case directory contains all the data that is required for
running a test within a directory, such as the genomic sequence,
various data files and configuration files. For example, for testing
:doc:`pipeline_rnaseq`, the following files are given:

   data files
      Brain-F1-R1.fastq.gz,Brain-F1-R2.fastq.gz,Brain-F2-R1.fastq.gz,Brain-F2-R2.fastq.gz,UHR-F1-R1.fastq.gz,UHR-F1-R2.fastq.gz
   
   genomic sequence
      hg19.fasta,hg19.idx

   bowtie indices
      hg19.1.ebwt,hg19.2.ebwt,hg19.3.ebwt,hg19.4.ebwt,hg19.fa,hg19.fa.fai,hg19.rev.1.ebwt,hg19.rev.2.ebwt
      hg19_cs.1.ebwt,hg19_cs.2.ebwt,hg19_cs.3.ebwt,hg19_cs.4.ebwt,hg19_cs.rev.1.ebwt,hg19_cs.rev.2.ebwt

   configuration files
      conf.py,pipeline.ini,sphinxreport.ini

Dependencies between pipelines are currently implemented within the
script itself. The dependencies that are currently implemented are:

* :doc:`pipeline_annotations` is run before any other pipelines.

Optional inputs
+++++++++++++++

Requirements
------------

All the test cases need to be properly configured. In particular, all
the required software for all the pipelines need to be in place in
order for all tests to pass.

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
import glob
import os
import shutil
import CGAT.Experiment as E

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import CGAT.Pipeline as P
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

PARAMS = P.PARAMS


def getPipelinesDir():
    '''return directory with pipelines.'''
    return os.path.dirname(__file__)

###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
def splitTestName(outfile):
    '''split the test name into pipeline and data part.

    The name of a test starts with the name of a pipeline followed by
    some description.

    '''

    pipeline_components = outfile.split("_")

    g = os.path.join(getPipelinesDir(), "pipeline_*.py")

    pipelines = set([os.path.basename(x) for x in glob.glob(g)])

    p = []
    for x in range(len(pipeline_components)):
        p.append(pipeline_components[x])
        if ("_".join(p) + ".py") in pipelines:
            pipeline_name = "_".join(p)
            test_description = "_".join(
                pipeline_components[x + 1:])
            break
    else:
        raise ValueError("pipeline for test %s not found" % outfile)

    return pipeline_name, test_description

###################################################################
###################################################################
###################################################################
def runTest(infile, outfile, update = False):
    '''run a test.'''

    test_name = P.snip(outfile, ".log")
    pipeline_name, test_description = splitTestName(test_name)

    # do not run on cluster, mirro
    # that a pipeline is set up on
    # the head node
    to_cluster = False

    pipelines_dir = getPipelinesDir()

    if not update:
        E.info("building test directory")
        try:
            shutil.rmtree("%s.dir" % test_name)
        except OSError:
            pass

        os.mkdir("%s.dir" % test_name)

        # create default config files
        statement = '''
        (cd %(test_name)s.dir;
        python %(pipelines_dir)s/%(pipeline_name)s.py config) > create.log
        '''
        P.run()

        # link to new files, overwriting default
        # configuration files if necessary.
        statement = '''ln -fs %(infile)s.dir/* %(test_name)s.dir/'''
        P.run()

    statement = '''
    (cd %(test_name)s.dir;
    python %(pipelines_dir)s/%(pipeline_name)s.py
    %(pipeline_options)s make full) >& %(outfile)s
    ''' 
    
    P.run()


###################################################################
###################################################################
###################################################################
## general tests
###################################################################
@files([(os.path.join(PARAMS["data_dir"], x + ".dir"), x + ".log")
        for x in P.asList(PARAMS["prerequisites"])])
def runPreparationTests(infile, outfile):
    '''run pre-requisite pipelines.'''
    runTest(infile, outfile)

###################################################################
###################################################################
###################################################################
## run a test
###################################################################
@follows(runPreparationTests)
@files([(x, os.path.basename(x) + ".log")
        for x in glob.glob(
            os.path.join(PARAMS["data_dir"], "pipeline_*"))
        if x not in P.asList(PARAMS["prerequisites"])])
def runTests(infile, outfile):
    '''run a pipeline with test data.'''
    runTest(infile, outfile)

###################################################################
###################################################################
###################################################################
## update a test
###################################################################
@follows(runPreparationTests)
@files([(x, os.path.basename(x) + ".log") for x in
        glob.glob(os.path.join(PARAMS["data_dir"], "pipeline_*"))])
def updateTests(infile, outfile):
    '''run a pipeline with test data.'''
    runTest(infile, outfile, update=True)

###################################################################
###################################################################
###################################################################
## build reports
###################################################################
@transform(runTests, suffix(".log"), ".report")
def runReports(infile, outfile):
    '''run a pipeline report.'''
    
    test_name = P.snip(outfile, ".report")

    pipeline_name, test_description = splitTestName(test_name)

    statement = '''
    (cd %(test_name)s.dir; python %(scriptsdir)s/%(pipeline_name)s.py
    %(pipeline_options)s make build_report) >& %(outfile)s
    ''' 
    
    P.run()

###################################################################
###################################################################
###################################################################
## primary targets
###################################################################
@follows(runTests, runReports)
def full():
    pass

###################################################################
###################################################################
###################################################################
## primary targets
###################################################################
@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''

    E.info("starting report build process from scratch")
    P.run_report(clean = True)

@follows(mkdir("report"))
def update_report():
    '''update report.'''

    E.info("updating report")
    P.run_report(clean = False)

@follows(update_report)
def publish_report():
    '''publish report.'''

    E.info("publishing report")
    P.publish_report()

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
