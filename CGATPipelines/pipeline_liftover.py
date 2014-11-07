"""
====================
Liftover pipeline
====================

:Author: Andreas Heger
:Release: $Id: pipeline_liftover.py 2900 2010-04-13 14:38:00Z andreas $
:Date: |today|
:Tags: Python

The liftover pipeline maps a set of intervals from one or more genomes
to a target genome. It uses the :term:`liftover` tool from UCSC.

Overview
========

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineReporting`). To start with, use the files supplied with the
Example_ data.

Input
-----

Requirements
------------

Pipeline output
===============

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_rnaseq.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_rnaseq.tgz
   tar -xvzf pipeline_rnaseq.tgz
   cd pipeline_rnaseq
   python <srcdir>/pipeline_rnaseq.py make full

.. note:: 
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::

   liftover
      ucsc_ tool to convert coordinates between assemblies

.. _ucsc: http://genome.ucsc.edu/


.. todo::
   * make the merging step optional. Currently overlapping intervals are merged.


Code
====



"""
import sys
import tempfile
import optparse
import shutil
import itertools
import csv
import math
import random
import re
import glob
import os
import shutil
import collections
import csv
import gzip
from ruffus import *
import sqlite3

import CGAT.Experiment as E
import CGAT.Pipeline as P
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

###################################################################
###################################################################
###################################################################
##
###################################################################
if os.path.exists("pipeline_conf.py"):
    L.info("reading additional configuration from pipeline_conf.py")
    execfile("pipeline_conf.py")

PARAMS = P.getParameters()

###################################################################
###################################################################
# Helper functions mapping tracks to conditions, etc
###################################################################
import CGATPipelines.PipelineTracks as PipelineTracks

TRACKS = PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
    glob.glob("*.gtf.gz"), "(\S+).gtf.gz", exclude=(".mapped.gtf.gz", ))

#####################################################################
#####################################################################
#####################################################################


@transform(TRACKS.getTracks("%s.gtf.gz"),
           suffix(".gtf.gz"),
           '.psl.gz')
def convertGtf2Psl(infile, outfile):
    """convert a gtf to a psl file.

    This method only takes features of type 'exon' and
    skips all contigs that are not in the genome sequence
    (for example the variant human chromosomes).
    """

    track = outfile[:-len(".psl.gz")]
    genomefile = os.path.join(
        PARAMS["genome_dir"], PARAMS["%s_genome" % track])
    if not os.path.exists(genomefile + ".fasta"):
        raise IOError("genome %s does not exist" % genomefile)

    statement = """gunzip 
    < %(infile)s
    | awk '$3 == "exon"'
    | python %(scriptsdir)s/gff2gff.py
    --method=sanitize
    --sanitize-method=genome
    --skip-missing
    --genome=%(genomefile)s
    --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2psl.py
    --allow-duplicates
    --is-gtf
    --log=%(outfile)s.log
    | gzip > %(outfile)s
    """
    P.run()

###################################################################


@transform('*.bed.gz',
           suffix(".bed.gz"),
           '.psl.gz')
def convertBed2Psl(infile, outfile):
    """convert a bed to a psl file."""

    track = outfile[:-len(".bed.gz")]
    genomefile = os.path.join(
        PARAMS["genome_dir"], PARAMS["%s_genome" % track])
    if not os.path.exists(genomefile + ".fasta"):
        raise IOError("genome %s does not exist" % genomefile)

    statement = """gunzip < %(infile)s 
    | python %(scriptsdir)s/bed2psl.py 
         --genome=%(genomefile)s
         --log=%(outfile)s.log 
    | gzip > %(outfile)s
    """
    P.run()

###################################################################
###################################################################
###################################################################


@transform((convertGtf2Psl, convertBed2Psl), suffix(".psl.gz"), '.transcripts')
def mergeTranscripts(infile, outfile):
    """merge transcripts before mapping.

    Overlapping transcripts are combined in order to
    speed up the mapping process.
    """

    track = outfile[:-len(".transcripts")]
    genomefile = os.path.join(
        PARAMS["genome_dir"], PARAMS["%s_genome" % track])

    statement = """
        gunzip < %(infile)s 
        | awk '/psLayout/ { x = 4; next; } x > 0 { --x; next} { print; }' 
        | sort -k 14,14 -k 16,16n
	| %(cmd-farm)s 
		--split-at-column=14 
		--output-header 
		--renumber="%%06i" 
		--renumber-column=":id" 
		--log=%(outfile)s.log 
		--subdirs \
        "python %(scriptsdir)s/psl2assembly.py 
               --staggered=all 
               --method=region
               --method=transcript
                --threshold-merge-distance=0 
                --threshold-merge-overlap=3
		--genome=%(genomefile)s
		--mali-output-format=fasta 
		--log=%(outfile)s.log 
		--output-filename-pattern=%%DIR%%%(outfile)s.%%s" 
	> %(outfile)s"""

    P.run()


@transform(mergeTranscripts, suffix(".transcripts"), '.merged.mapped.psl')
def mapMergedTranscripts(infile, outfile):
    """map transcripts from PSL.

    Mapping from PSL is equivalent to first converting to genePred format 
    and using the option -gp.
    """

    track = outfile[:-len(".merged.mapped.psl")]
    chainfile = os.path.join(PARAMS["ucsc_dir"],
                             PARAMS["%s_genome" % track],
                             "liftOver",
                             "%sTo%s.over.chain.gz" %
                             (PARAMS["%s_genome" % track],
                              PARAMS["genome"].capitalize()))
    if not os.path.exists(chainfile):
        raise IOError("chain file %s does not exist" % chainfile)

    statement = """
        liftOver -minMatch=0.2 -minBlocks=0.01 -pslT 
                 %(infile)s.transcripts.psl 
                 <(gunzip < %(chainfile)s) 
                 %(outfile)s 
                 %(outfile)s.unmapped 
        >& %(outfile)s.log
        """
    P.run()


@transform((convertGtf2Psl, convertBed2Psl), suffix(".psl.gz"), '.mapped.psl')
def mapTranscripts(infile, outfile):
    """map transcripts from PSL.

    Mapping from PSL is equivalent to first converting to genePred format 
    and using the option -gp.
    """

    track = outfile[:-len(".mapped.psl")]
    chainfile = os.path.join(PARAMS["ucsc_dir"],
                             PARAMS["%s_genome" % track],
                             "liftOver",
                             "%sTo%s.over.chain.gz" %
                             (PARAMS["%s_genome" % track],
                              PARAMS["genome"].capitalize()))

    statement = """
        liftOver -minMatch=0.2 -minBlocks=0.01 -pslT 
                 <(gunzip < %(infile)s )
                 <(gunzip < %(chainfile)s) 
                 %(outfile)s 
                 %(outfile)s.unmapped 
        >& %(outfile)s.log
        """
    P.run()


@transform((mapMergedTranscripts, mapTranscripts), suffix(".psl"), '.gtf.gz')
def convertMappedPslToGtf(infile, outfile):
    '''convert to gtf for export.'''
    statement = """
    python %(scriptsdir)s/psl2gff.py --as-gtf 
    < %(infile)s 
    | gzip
    > %(outfile)s
    """
    P.run()


@transform(convertMappedPslToGtf, suffix(".gtf.gz"), '.summary')
def summary(infile, outfile):
    '''compute mapping stats.'''

    def _getfiles(filename):

        track = outfile[:-len(".mapped.summary")]
        if track.endswith(".merged"):
            xtrack = track[:-len(".merged")]
            finput = "%s.psl.gz" % xtrack
            fmerged = "%s.transcripts.transcripts.psl" % xtrack
            fmapped = "%s.mapped.psl" % track
        else:
            finput = "%s.psl.gz" % track
            fmerged = finput
            fmapped = "%s.mapped.psl" % track
        return track, finput, fmerged, fmapped

    outf = open(outfile, "w")
    outf.write("track\tinput\tmerged\tpmerged\tmapped\tpmapped\tpoutput\n")

    def countPSL(filename):
        if filename.endswith(".gz"):
            i = gzip.open(filename)
        else:
            i = open(filename)
        ll = [x[:10] for x in i.readlines() if not x.startswith("#")]
        if ll[0].startswith("psLayout"):
            return len(ll) - 5
        else:
            return len(ll)

    track, finput, fmerged, fmapped = _getfiles(outfile)
    ninput = countPSL(finput)
    # subtract header
    nmerged = countPSL(fmerged) - 5
    nmapped = countPSL(fmapped)

    outf.write("%s\t%i\t%i\t%s\t%i\t%s\t%s\n" %
               (track,
                ninput,
                nmerged,
                IOTools.prettyPercent(nmerged, ninput),
                nmapped,
                IOTools.prettyPercent(nmapped, nmerged),
                IOTools.prettyPercent(nmapped, ninput)))


@follows(convertMappedPslToGtf, summary)
def full():
    pass

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
