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

This pipeline implements automated testing of CGAT pipelines. The
pipeline downloads test data from a specified URL, runs the associated
pipeline for each data set and compares the output with a reference.
The results are collected in a report.

Tests are setup in the pipeline configuration file.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

In order to run all tests, simply enter an empty directory and type::

   python <srcdir>/pipeline_testing.py make config

Edit the config files as required and then type::

   python <srcdir>/pipeline_testing.py make full
   python <srcdir>/pipeline_testing.py make build_report

The first command will download the data and run the pipelines while
the second will build a summary report.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.

Tests are described as section in the configuration file. A test
section starts with the prefix ``test_``. For following example is a
complete test setup::

   [test_mytest1]
   # pipeline to run
   pipeline=pipeline_mapping

   # pipeline target to run (default is 'full')
   # multiple targets can be specified as a comma separated list.
   target=full

   # filename suffixes to test
   suffixes=gtf.gz,bed.gz,tsv.gz,bam,nreads

   # regular expression of files to be excluded from
   # test for difference. Use | to separate multiple
   # regular expressions.
   regex_no_diff=rates.gff.gz

This configuration will run the test ``mytest1``. The associated
pipeline is :doc:`pipeline_mapping` and it will execute the target
``make full``. To check if the pipeline has completed successfully, it
will compare all files ending with any of the suffixes specified
(``gtf.gz``, ``bed.gz``, etc). The comparison will be done by building
a checksum of the whole file ignoring any comments (lines starting
with a ``#``).

Some files will be different at every run, for example if they use
some form of random initialization. Files matching the regular
expressions in ``regex_only_exist`` will thus only be tested if they
exist or not.

The test expects a file called :file:`test_mytest1.tgz` with the
test data at the download URL (parameter ``data_url``).

To define a default test for a pipeline, simply name the
test ``test_<pipeline name>``, for example::

   [test_mapping]
   suffixes=gtf.gz,bed.gz,tsv.gz,bam,nreads

Note that setting the ``target`` and ``pipeline`` options is
not necessary in this case as the default values suffice.

Input data
----------

The input data for each test resides a compressed tar-ball. The input
data should uncompress in a directory called :file:`<testname>.dir`
The tar-ball need also contain a file :file:`<testname>.ref`
containing the md5 checksums of files of a previous run of the test
that is being used as a reference.

The input data should contain all the data that is required for
running a test within a directory. It is best to minimize dependencies
between tests, though there is a mechanism for this (see below).

For example, the contents of a tar-ball will look light this::

   test_mytest1.dir/                     # test data root
   test_mytest1.dir/Brain-F2-R1.fastq.gz # test data
   test_mytest1.dir/Brain-F1-R1.fastq.gz
   test_mytest1.dir/hg19.fasta           # genomic data
   test_mytest1.dir/hg19.idx
   test_mytest1.dir/hg19.fa
   test_mytest1.dir/hg19.fa.fai
   test_mytest1.dir/pipeline.ini  # pipeline configuration file
   test_mytest1.dir/indices/      # configured to work in test dir
   test_mytest1.dir/indices/bwa/  # bwa indices
   test_mytest1.dir/indices/bwa/hg19.bwt
   test_mytest1.dir/indices/bwa/hg19.ann
   test_mytest1.dir/indices/bwa/hg19.pac
   test_mytest1.dir/indices/bwa/hg19.sa
   test_mytest1.dir/indices/bwa/hg19.amb
   test_mytest1.ref   # reference file

The reference file looks like this::

   test_mytest1.dir/bwa.dir/Brain-F2-R2.bwa.bam 0e1c4ee88f0249c21e16d93ac496eddf  
   test_mytest1.dir/bwa.dir/Brain-F1-R2.bwa.bam 01bee8af5bbb5b1d13ed82ef1bc3620d  
   test_mytest1.dir/bwa.dir/Brain-F2-R1.bwa.bam 80902c87519b6865a9ca982787280972  
   test_mytest1.dir/bwa.dir/Brain-F1-R1.bwa.bam 503c99ab7042a839e56147fb1a221f27  
   ...

This file is created by the test pipeline and called
:file:`test_mytest1.md5`.  When setting up a test, start with an empty
files and later add this file to the test data.

Pipeline dependencies
---------------------

Some pipelines depend on the output of other pipelines, most notable
is :doc:`pipeline_annotations`. To run a set of pipelines before other
pipelines name them in the option ``prerequisites``, for example::

   prerequisites=test_annnotations

Pipeline output
===============

The major output is in the database file :file:`csvdb`.

Code
====

"""
from ruffus import *

import sys
import pipes
import os
import re
import tarfile
import CGAT.Experiment as E
import CGAT.IOTools as IOTools

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


@files([(None, "%s.tgz" % x)
        for x in P.CONFIG.sections()
        if x.startswith("test")])
def setupTests(infile, outfile):
    '''setup tests.

    This method creates a directory in which a test will be run
    and downloads test data with configuration files.
    '''
    to_cluster = False

    track = P.snip(outfile, ".tgz")

    if os.path.exists(track + ".dir"):
        raise OSError('directory %s.dir already exists' % track)

    # create directory
    os.mkdir(track + ".dir")

    # run pipeline config
    pipeline_name = PARAMS.get(
        "%s_pipeline" % track,
        "pipeline_" + track[len("test_"):])

    statement = '''
    (cd %(track)s.dir;
    python %(pipelinedir)s/%(pipeline_name)s.py
    %(pipeline_options)s config) >& %(outfile)s.log
    '''
    P.run()

    # obtain data - should overwrite pipeline.ini file
    statement = '''
    wget --no-check-certificate -O %(track)s.tgz %(data_url)s/%(track)s.tgz'''
    P.run()

    tf = tarfile.open(outfile)

    tf.extractall()

    if not os.path.exists("%s.dir" % track):
        raise ValueError(
            "test package did not create directory '%s.dir'" % track)


def runTest(infile, outfile):
    '''run a test.'''

    track = P.snip(outfile, ".log")

    pipeline_name = PARAMS.get(
        "%s_pipeline" % track,
        "pipeline_" + track[len("test_"):])
    pipeline_targets = ' '.join(P.asList(
        PARAMS.get("%s_target" % track,
                   "full")))

    # do not run on cluster, mirror
    # that a pipeline is started from
    # the head node
    to_cluster = False

    statement = '''
    (cd %(track)s.dir;
    python %(pipelinedir)s/%(pipeline_name)s.py
    %(pipeline_options)s make %(pipeline_targets)s) >& %(outfile)s
    '''
    P.run()


@follows(setupTests)
@files([("%s.tgz" % x, "%s.log" % x)
        for x in P.asList(PARAMS.get("prerequisites", ""))])
def runPreparationTests(infile, outfile):
    '''run pre-requisite pipelines.'''
    runTest(infile, outfile)


@follows(runPreparationTests)
@files([("%s.tgz" % x, "%s.log" % x)
        for x in P.CONFIG.sections()
        if x.startswith("test") and
        x not in P.asList(PARAMS.get("prerequisites", ""))])
def runTests(infile, outfile):
    '''run a pipeline with test data.'''
    runTest(infile, outfile)


@transform((runPreparationTests, runTests),
           suffix(".log"),
           ".md5")
def buildCheckSums(infile, outfile):
    '''build checksums for files in the build directory.

    Files are uncompressed before computing the checksum
    as gzip stores meta information such as the time stamp.
    '''

    track = P.snip(infile, ".log")

    suffixes = P.asList(PARAMS.get(
        '%s_suffixes' % track,
        PARAMS["suffixes"]))

    if len(suffixes) == 0:
        raise ValueError('no file types defined for test')

    regex_pattern = ".*\(%s\)" % "\|".join(suffixes)

    regex_pattern = pipes.quote(regex_pattern)

    # ignore log files as time stamps will
    # be different
    statement = '''find %(track)s.dir
    -type f
    -not -regex ".*.log"
    -regex %(regex_pattern)s
    -exec %(scriptsdir)s/cgat_md5sum.sh {} \;
    >> %(outfile)s'''
    P.run()


@merge(buildCheckSums,
       "md5_compare.tsv")
def compareCheckSums(infiles, outfile):
    '''compare checksum files against
    existing reference data.
    '''

    outf = IOTools.openFile(outfile, "w")
    outf.write("\t".join((
        ("track", "status",
         "nfiles", "nref",
         "missing", "extra", "different",
         "files_missing", "files_extra", "files_different"))) + "\n")

    for infile in infiles:
        track = P.snip(infile, ".md5")
        reffile = track + ".ref"

        # regular expression of files to test only for existence
        regex_no_diff = PARAMS.get('%s_regex_no_diff' % track, None)
        if regex_no_diff:
            regex_no_diff = re.compile(regex_no_diff)

        if not os.path.exists(reffile):
            raise ValueError('no reference data defined for %s' % track)

        cmp_data = [x.split()[:2] for x in IOTools.openFile(infile)]
        ref_data = [x.split()[:2] for x in IOTools.openFile(reffile)]
        cmp_md5 = dict(cmp_data)
        ref_md5 = dict(ref_data)
        shared_files = set(cmp_md5).intersection(ref_md5)

        missing = set(ref_md5).difference(cmp_md5)
        extra = set(cmp_md5).difference(ref_md5)
        different = [x for x in shared_files if ref_md5[x] != cmp_md5[x]]
        # remove any files not to be checked for difference
        if regex_no_diff:
            different = [x for x in different
                         if not regex_no_diff.search(x)]

        if len(missing) + len(extra) + len(different) == 0:
            status = "OK"
        else:
            status = "FAIL"

        outf.write("\t".join(map(str, (
            track,
            status,
            len(cmp_md5),
            len(ref_md5),
            len(missing),
            len(extra),
            len(different),
            ",".join(missing),
            ",".join(extra),
            ",".join(different)))) + "\n")

    outf.close()


@transform(compareCheckSums,
           suffix(".tsv"),
           ".load")
def loadComparison(infile, outfile):
    '''load comparison data into database.'''
    P.load(infile, outfile)


@transform((runPreparationTests, runTests),
           suffix(".log"),
           ".report")
def runReports(infile, outfile):
    '''run a pipeline report.'''

    track = P.snip(outfile, ".report")

    pipeline_name = PARAMS.get(
        "%s_pipeline" % track,
        "pipeline_" + track[len("test_"):])

    statement = '''
    (cd %(track)s.dir; python %(pipelinedir)s/%(pipeline_name)s.py
    %(pipeline_options)s make build_report) >& %(outfile)s
    '''

    P.run()


@follows(runTests, runReports, loadComparison)
def full():
    pass


@files(None, 'reset.log')
def reset(infile, outfile):
    '''remove all data in pipeline.'''
    to_cluster = False

    statement = '''
    rm -rf test_* _cache _static _templates _tmp report;
    rm -f *.log csvdb *.load *.tsv'''
    P.run()

###################################################################
###################################################################
###################################################################
# primary targets
###################################################################


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
    sys.exit(P.main(sys.argv))
