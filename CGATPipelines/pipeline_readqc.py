##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Tildon Grant Belgard
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
##########################################################################

"""
====================
ReadQc pipeline
====================

:Author: David Sims
:Release: $Id$
:Date: |today|
:Tags: Python

The readqc pipeline imports unmapped reads from one or more
fastq and performs basic quality control steps:

   1. per position quality
   2. per read quality
   3. duplicates

For further details see http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/

Additionaly, optional bias analysis can be included through the configuration
file. This analysis is designed to help identify sequence contexts which bias
gene expression and asssess the consistency in the biases between samples.

Bias analysis utilises Sailfish to estimate transcript abundance. This
requires a multi-fasta transcripts file.

For further details see http://www.cs.cmu.edu/~ckingsf/software/sailfish/


Individual tasks are enabled in the configuration file.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning`
on general information how to use CGAT pipelines.

Configuration
-------------

No general configuration required.

Input
-----

Reads are imported by placing files or linking to files in the :term:
`working directory`.

The default file format assumes the following convention:

   <sample>-<condition>-<replicate>.<suffix>

``sample`` and ``condition`` make up an :term:`experiment`,
while ``replicate`` denotes the :term:`replicate` within an :term:`experiment`.
The ``suffix`` determines the file type.
The following suffixes/file types are possible:

sra
   Short-Read Archive format. Reads will be extracted using the :file:
   `fastq-dump` tool.

fastq.gz
   Single-end reads in fastq format.

fastq.1.gz, fastq2.2.gz
   Paired-end reads in fastq format.
   The two fastq files must be sorted by read-pair.

.. note::

   Quality scores need to be of the same scale for all input files.
   Thus it might be difficult to mix different formats.

Pipeline output
===============

The major output is a set of HTML pages and plots reporting on the quality of
the sequence archive

Example
=======

Example data is available at
http://www.cgat.org/~andreas/sample_data/pipeline_readqc.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_readqc.tgz
   tar -xvzf pipeline_readqc.tgz
   cd pipeline_readqc
   python <srcdir>/pipeline_readqc.py make full

Requirements:




Code
====

"""

###################################################
###################################################
###################################################
# load modules
###################################################

# import ruffus
from ruffus import *

# import useful standard python modules
import sys
import os
import re
import glob
import cStringIO
import numpy
import pandas
from pandas import DataFrame
from scipy.stats import linregress
import itertools as iter

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGATPipelines.PipelineMapping as PipelineMapping
import CGAT.Pipeline as P
import CGAT.CSV2DB as CSV2DB

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# load options from the config file
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])
PARAMS = P.PARAMS

# Include optional bias analysis for RNA-Seq
BIAS_ANALYSIS = P.isTrue("bias_analysis")
#########################################################################
#########################################################################
#########################################################################
# define input files
INPUT_FORMATS = ("*.fastq.1.gz", "*.fastq.gz", "*.sra", "*.csfasta.gz")
REGEX_FORMATS = regex(r"(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz)")

#########################################################################
#########################################################################
#########################################################################


@follows(mkdir(PARAMS["exportdir"]), mkdir(os.path.join(PARAMS["exportdir"],
                                                        "fastqc")))
@transform(INPUT_FORMATS,
           REGEX_FORMATS,
           r"\1.fastqc")
def runFastqc(infiles, outfile):
    '''convert sra files to fastq and check mapping qualities are in solexa format.
    Perform quality control checks on reads from .fastq files.'''
    m = PipelineMapping.FastQc(nogroup=PARAMS["readqc_no_group"],
                               outdir=PARAMS["exportdir"]+"/fastqc")
    statement = m.build((infiles,), outfile)
    P.run()

#########################################################################
#########################################################################
#########################################################################
##
#########################################################################


def FastqcSectionIterator(infile):
    data = []
    for line in infile:
        if line.startswith(">>END_MODULE"):
            yield name, status, header, data
        elif line.startswith(">>"):
            name, status = line[2:-1].split("\t")
            data = []
        elif line.startswith("#"):
            header = "\t".join([x for x in line[1:-1].split("\t") if x != ""])
        else:
            data.append(
                "\t".join([x for x in line[:-1].split("\t") if x != ""]))


@jobs_limit(1, "db")
@transform(runFastqc, suffix(".fastqc"), "_fastqc.load")
def loadFastqc(infile, outfile):
    '''load FASTQC stats.'''

    track = P.snip(infile, ".fastqc")

    filename = os.path.join(
        PARAMS["exportdir"], "fastqc", track + "*_fastqc", "fastqc_data.txt")

    for fn in glob.glob(filename):
        prefix = os.path.basename(os.path.dirname(fn))
        results = []

        for name, status, header, data in FastqcSectionIterator(
                IOTools.openFile(fn)):
            # do not collect basic stats, see loadFastQCSummary
            if name == "Basic Statistics":
                continue

            parser = CSV2DB.buildParser()
            (options, args) = parser.parse_args([])
            options.tablename = prefix + "_" + re.sub(" ", "_", name)
            options.allow_empty = True

            inf = cStringIO.StringIO("\n".join([header] + data) + "\n")
            CSV2DB.run(inf, options)
            results.append((name, status))

        # load status table
        parser = CSV2DB.buildParser()
        (options, args) = parser.parse_args([])
        options.tablename = prefix + "_status"
        options.allow_empty = True

        inf = cStringIO.StringIO(
            "\n".join(["name\tstatus"] +
                      ["\t".join(x) for x in results]) + "\n")
        CSV2DB.run(inf, options)

    P.touch(outfile)


def collectFastQCSections(infiles, section):
    '''iterate over all fastqc files and extract a particular section.'''
    results = []

    for infile in infiles:

        track = P.snip(infile, ".fastqc")

        filename = os.path.join(
            PARAMS["exportdir"], "fastqc", track + "*_fastqc",
            "fastqc_data.txt")

        for fn in glob.glob(filename):
            prefix = os.path.basename(os.path.dirname(fn))
            for name, status, header, data in FastqcSectionIterator(
                    IOTools.openFile(fn)):
                if name == section:
                    results.append((track, status, header, data))

    return results


@merge(runFastqc, "status_summary.tsv.gz")
def buildFastQCSummaryStatus(infiles, outfile):
    '''load fastqc status summaries into a single table.'''

    outf = IOTools.openFile(outfile, "w")
    first = True
    for infile in infiles:
        track = P.snip(infile, ".fastqc")
        filename = os.path.join(
            PARAMS["exportdir"], "fastqc", track + "*_fastqc",
            "fastqc_data.txt")

        for fn in glob.glob(filename):
            prefix = os.path.basename(os.path.dirname(fn))
            results = []

            names, stats = [], []
            for name, status, header, data in FastqcSectionIterator(
                    IOTools.openFile(fn)):
                stats.append(status)
                names.append(name)

            if first:
                outf.write("track\tfilename\t%s\n" % "\t".join(names))
                first = False

            outf.write("%s\t%s\t%s\n" %
                       (track, os.path.dirname(fn), "\t".join(stats)))
    outf.close()


@merge(runFastqc, "basic_statistics_summary.tsv.gz")
def buildFastQCSummaryBasicStatistics(infiles, outfile):
    '''load fastqc summaries into a single table.'''

    data = collectFastQCSections(infiles, "Basic Statistics")

    outf = IOTools.openFile(outfile, "w")
    first = True
    for track, status, header, rows in data:
        rows = [x.split("\t") for x in rows]
        if first:
            headers = [row[0] for row in rows]
            outf.write("track\t%s\n" % "\t".join(headers))
            first = False
        outf.write("%s\t%s\n" % (track, "\t".join([row[1] for row in rows])))
    outf.close()


@transform((buildFastQCSummaryStatus, buildFastQCSummaryBasicStatistics),
           suffix(".tsv.gz"), ".load")
def loadFastqcSummary(infile, outfile):
    P.load(infile, outfile, options="--add-index=track")


#########################################################################


@follows(loadFastqc, loadFastqcSummary)
def full():
    pass


#########################################################################


@follows()
def publish():
    '''publish files.'''
    P.publish_report()


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''

    E.info("starting documentation build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.'''

    E.info("updating documentation")
    P.run_report(clean=False)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
