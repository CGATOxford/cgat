###############################################################################
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
#################################################################################
"""
================================
Optimizing cufflinks parameters
================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

The cufflinks optimization pipeline attempts to assess the quality of a variety of transcript
assemblies using various sets of parameters.


Overview
==========

Building transcripts is an important prerequisite for a number of projects that we undertake,
especially where the emphasis is on transcript structure or the identification of non-coding
transcripts (often lowly expressed).

In addition to the general issues  associated with transcript building, we are often working on
a variety of data from different library types and as such are not aware of the optimal set of
parameters for a transcript assembly from the outset of an analysis. This is especially important when 
we are faced with such things as high numbers of intronic reads from a non-polyA selected library
or low depth. Ideally we want to optimize the number of spliced reads that are incorporated into
transcript models, thus using the maximal amount of information that is present in any given datset.


Another issue in transcript assembly is the time that is required to run an assembly on an entire
transcriptome, creating a bottleneck in parameter optimization.

Given the above, the cufflinks optimization pipeline uses a reduced set of alignments (from a tophat run)
to optimize a variety of user specified parameters. It performs the following tasks

    * Reduction of input bam files into chr19 only bamfiles

    * Runs the rnaseq transcript building pipeline for all combinations
    of user specified parameters

    * Collects metrics in a single database csvdb from each run (metrics provided by the transcript
    building pipeline) and assesses how many reads contribute to transcripts (inc. spliced reads)

    * It also assesses the ratio of single vs. multi exon transfrags as a measure of overall transcriptome
    quality.


Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.

Input
-----

Mapped reads
++++++++++++

The principal input of this pipeline is a collection of reads mapped to a reference genome.
Mapped reads are imported by placing files are linking to files in the :term:`working directory`.

The default file format assumes the following convention:

   <sample>-<condition>-<replicate>.bam

``sample`` and ``condition`` make up an :term:`experiment`, while ``replicate`` denotes
the :term:`replicate` within an :term:`experiment`. 

  
Requirements
------------

The pipeline requires the results from :doc:`pipeline_annotations`. Set the configuration variable 
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|cufflinks_          |>=1.3.0            |transcription levels                            |
+--------------------+-------------------+------------------------------------------------+
|samtools            |>=0.1.16           |bam/sam files                                   |
+--------------------+-------------------+------------------------------------------------+
|bedtools            |                   |working with intervals                          |
+--------------------+-------------------+------------------------------------------------+


Code
====

"""

# load modules
from ruffus import *

import Experiment as E
import logging as L
import Database, CSV

import sys, os, re, shutil, itertools, math, glob, time, gzip, collections, random, operator

import numpy, sqlite3
import GFF, GTF, IOTools, IndexedFasta
import Tophat
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError
import pysam

import Expression

import PipelineGeneset
import PipelineMapping
import PipelineRnaseq
import PipelineMappingQC
import Stats

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import Pipeline as P
P.getParameters( 
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini" ],
    defaults = {
        'annotations_dir' : "",
        'paired_end' : False } )

PARAMS = P.PARAMS

PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
                                       "pipeline_annotations.py" )

###################################################################
###################################################################
###################################################################
# get options that are to be tested

options = P.asList(PARAMS["cufflinks_test_options"])

cufflinks_options = {}
for option in options:
    if option == "--min-isoform-fraction" \
            or option == "--pre-mrna-fraction" \
            or option == "--small-anchor-fraction" \
            or option == "--max-multiread-fraction": 
        cufflinks_options[option]=[0.1, 0.5, 0.75]
    elif option == "--junc-alpha":
        cufflinks_options[option]=[0.001, 0.01, 0.1]
    elif option == "--min-frags-per-transfrag":
        cufflinks_options[option]=[0,2,5,10]
    elif option == "--max-multiread-fraction":
        cufflinks_options[option]=[0.1, 0.5, 0.75]
    elif option == "--overhang-tolerance":
        cufflinks_options[option]=[0,2,5,8]
    elif option == "--overlap-radius":
        cufflinks_options[option]=[50, 100]
    else:
        raise ValueError("pipeline_cufflinks_optimization does not support parameter %s" % option)
c = 0
for x in itertools.product(*cufflinks_options.values()):
    c += 1
raw_input("pipeline_cufflinks_optimization will run %i trancript assemblies: hit enter to continue\n""" % c)

###################################################################
###################################################################
###################################################################
def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect( PARAMS["database"] )
    statement = '''ATTACH DATABASE '%s' as annotations''' % (PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute( statement )
    cc.close()

    return dbh

###################################################################
###################################################################
###################################################################
def updateFile(filename):
    '''
    create empty file for updating purposes
    '''
    
    outf = open(filename, "w")
    outf.write("file created for ruffus update")
    outf.close()

###################################################################
###################################################################
###################################################################
def options_generator(cufflinks_options):
    '''
    returns a generator for dealing with the cufflinks parameters
    for directory names and .ini file creation
    '''
    for option_values in itertools.product(*cufflinks_options.values()):
        yield " ".join(map(str,reduce(operator.add, zip(cufflinks_options.keys(), list(option_values)))))

###################################################################
###################################################################
###################################################################
def getDirectoryNames(options_generator):
    '''
    get the names for the directories
    '''
    for options in options_generator:
        yield "_".join(options.split(" ")).replace("--", "")

###################################################################
###################################################################
###################################################################
def getLogFileNames(options_generator):
    '''
    get filenames for log files
    '''
    for options in options_generator:
        yield "_".join(options.split(" ")).replace("--", "") + ".log"
        
###################################################################
###################################################################
###################################################################
@transform(glob.glob("*.accepted.bam"), suffix(".bam"), ".chr19.bam")
def reduceBamToChr19(infile, outfile):
    '''
    reduce the dataset for parameter testing.
    '''
    bam = pysam.Samfile(infile, "rb")
    outbam = pysam.Samfile(outfile, "wb", template = bam)
    for alignment in bam.fetch("chr19"):
        outbam.write(alignment)
    bam.close()
    outbam.close()

###################################################################
###################################################################
###################################################################
@transform(reduceBamToChr19, suffix(".bam"), ".bam.bai")
def indexBam(infile, outfile):
    '''
    index the chr19 bam file
    '''
    pysam.index(infile)
    
###################################################################
###################################################################
###################################################################
@follows(indexBam, mkdir([directory for directory in getDirectoryNames(options_generator(cufflinks_options))]))
@split(indexBam, [logfile for logfile in getLogFileNames(options_generator(cufflinks_options))])
def createLogFiles(infile, outfiles):
    '''
    creates sentinel files
    '''
    for outfile in outfiles:
        updateFile(outfile)

###################################################################
###################################################################
###################################################################
@follows(createLogFiles)
@transform(indexBam, suffix(".bai"), add_inputs([logfile for logfile in getLogFileNames(options_generator(cufflinks_options))]), ".log")
def linkBamToWorkingDirs(infiles, outfile):
    '''
    symlink the bam file and index to the working directories
    for execution of the transcript building pipeline
    '''

    bamfile = P.snip(infiles[0], ".bai")
    indexfile = infiles[0]
    directories = [P.snip(logfile, ".log") for logfile in infiles[1]]

    for directory in directories:
        os.symlink(os.path.abspath(bamfile), os.path.join(directory, bamfile))
        os.symlink(os.path.abspath(indexfile), os.path.join(directory, indexfile))
    updateFile(outfile)
    
###################################################################
###################################################################
###################################################################
@transform(createLogFiles, regex(r"(\S+).log")
           , r"\1/pipeline.ini")
def createConfigFiles(infile, outfile):
    '''
    create all of the relevant .ini files in each working
    directory in order to execute the transcript building
    '''
    # test options for cufflinks
    cuff_opts = P.snip(infile, ".log").split("_")
    cuff_options = []
    for opt in cuff_opts:
        if len(opt)>4: # not ideal to do my length but all I can think of at the moment
            cuff_options.append("--" + opt)
        else:
            cuff_options.append(opt)
    cuff_options = " ".join(cuff_options)
    
    # directory for output config
    outdir = P.snip(infile, ".log")

    outf = open(os.path.join(outdir, "pipeline.ini"), "w")
    config_headers = []
    for line in open(glob.glob("*.ini")[0]).readlines():
        if line.find("[cufflinks]") != -1:
            outf.write( "[cufflinks]\n\n# general cufflinks options\n\noptions=--upper-quartile-norm %s   \n" % cuff_options )
        else:
            outf.write(line)
    outf.close()
        
###################################################################
###################################################################
###################################################################
@follows(linkBamToWorkingDirs)
@transform(createConfigFiles, regex(r"(\S+).ini"), r"\1.log")
def executePipelineRnaseqTranscripts(infile, outfile):
    '''
    executes the transcripts building pipeline in each directory
    '''
    directory = os.path.dirname(infile)
    statement = '''cd ./%(directory)s; 
                   python %(scriptsdir)s/pipeline_rnaseqtranscripts.py -v5 -p10 make full'''
    P.run()

###################################################################
###################################################################
###################################################################
@files([os.path.join(PARAMS["general_pipeline_rnaseqtranscripts_dir"], x) for x in ["conf.py", "sphinxreport.ini"]]
           , [os.path.join(directoryname, y) for directoryname, y in itertools.product (getDirectoryNames(options_generator(cufflinks_options)), ["conf.py", "sphinxreport.ini"])])
def linkToPipelineRnaseqTranscriptsConfigFiles(infiles, outfiles):
    '''
    produces links to the configuration files for report building
    '''
    
    
    print infiles
    print outfiles

###################################################################
###################################################################
###################################################################
@follows(executePipelineRnaseqTranscripts)
def full():
    pass

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )

    

     
        


         
    
    
    
    

    

