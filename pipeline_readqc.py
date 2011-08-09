################################################################################
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

The pipeline can also be used to pre-process reads. Implemented tasks are:

   * :meth:`removeContaminants` - remove contaminants from read set

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

No general configuration required.

Removing contaminants
---------------------

Use the task :meth:`removeContaminants` to remove contaminants from read
sets.

Contaminant sequences are listed in the file :file:`contaminants.fasta`. 
If not given, a file with standard Illumina adapators will be created 
to remove adaptor contamination.

The task will create output files called :file:`nocontaminants-<infile>`.

The pipeline can then be re-run in order to add stats on the contaminant-removed
files.

.. note::

   Colour space filtering has not been implemented yet.

Input
-----

Reads are imported by placing files or linking to files in the :term:`working directory`.

The default file format assumes the following convention:

   <sample>-<condition>-<replicate>.<suffix>

``sample`` and ``condition`` make up an :term:`experiment`, while ``replicate`` denotes
the :term:`replicate` within an :term:`experiment`. The ``suffix`` determines the file type.
The following suffixes/file types are possible:

sra
   Short-Read Archive format. Reads will be extracted using the :file:`fastq-dump` tool.

fastq.gz
   Single-end reads in fastq format.

fastq.1.gz, fastq2.2.gz
   Paired-end reads in fastq format. The two fastq files must be sorted by read-pair.

.. note::

   Quality scores need to be of the same scale for all input files. Thus it might be
   difficult to mix different formats.

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|fastqc              |>=0.9.0            |read quality control                            |
+--------------------+-------------------+------------------------------------------------+
|sra-tools           |                   |extracting reads from .sra files                |
+--------------------+-------------------+------------------------------------------------+
|picard              |>=1.38             |bam/sam files. The .jar files need to be in your|
|                    |                   | CLASSPATH environment variable.                |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The major output is a set of HTML pages and plots reporting on the quality of the sequence archive

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_readqc.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_readqc.tgz
   tar -xvzf pipeline_readqc.tgz
   cd pipeline_readqc
   python <srcdir>/pipeline_readqc.py make full


Code
====

"""

# load modules
from ruffus import *
from rpy2.robjects import r as R

import Experiment as E
import logging as L
import Database
import sys, os, re, shutil, itertools, math, glob, time, gzip, collections, random
import numpy, sqlite3
import GTF, IOTools, IndexedFasta, FastaIterator
import Tophat
import rpy2.robjects as ro
import PipelineGeneset
import PipelineMapping
import Stats
import PipelineTracks
import Pipeline as P

USECLUSTER = True

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import Pipeline as P
P.getParameters( 
    ["%s.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini" ] )
PARAMS = P.PARAMS

#########################################################################
#########################################################################
#########################################################################
@follows(mkdir(PARAMS["exportdir"]), mkdir(os.path.join(PARAMS["exportdir"], "fastqc")) )
@transform( ("*.fastq.1.gz", 
              "*.fastq.gz",
              "*.sra",
	      "*.csfasta.gz"),
              regex( r"(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz)"),
              r"\1.fastqc")
def runFastqc(infiles, outfile):
    '''convert sra files to fastq and check mapping qualities are in solexa format. 
    Perform quality control checks on reads from .fastq files.'''
    to_cluster = USECLUSTER
    m = PipelineMapping.fastqc()
    statement = m.build((infiles,), outfile) 
    P.run()

#########################################################################
#########################################################################
#########################################################################
## adapter trimming
#########################################################################
# see http://intron.ccam.uchc.edu/groups/tgcore/wiki/013c0/Solexa_Library_Primer_Sequences.html
ILLUMINA_ADAPTORS = { "Genomic/ChIPSeq-Adapters1-1" : "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",
		      "Genomic/ChIPSeq-Adapters1-2" : "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		      "Genomic/ChIPSeq-PCR-1" : "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		      "Genomic/ChIPSeq-PCR-2" : "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT",
		      "Genomic/ChIPSeq-Adapters1-Genomic" : "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		      "Paired-End-Adapters-1" : "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",
		      "Paired-End-Adapters-2" : "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		      "Paired-End-PCR-1" : "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		      "Paired-End-PCR-2" : "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",
		      "Paired-End-Sequencing-1" : "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		      "Paired-End-Sequencing-2" : "CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT" }


@merge( None, "contaminants.fasta" )
def outputContaminants( infile, outfile ):
    '''output file with contaminants.'''
    outf = IOTools.openFile( outfile, "w")
    for key, value in ILLUMINA_ADAPTORS.iteritems():
        outf.write(">%s\n%s\n" % (key, value) )
    outf.close()

@transform( ("*.fastq.gz",
             "*.fastq.1.gz",
             "*.fastq.2.gz"),
	    regex( r"(?!nocontaminants)(\S+).(fastq.1.gz|fastq.gz|fastq.2.gz|csfasta.gz)"),
	    add_inputs(outputContaminants),
	    r"nocontaminants.\1.\2")
def removeContaminants( infiles, outfile ):
    '''remove adaptor contamination from fastq files.

    
    This method uses cutadapt.
    '''
    
    infile, contaminant_file = infiles

    adaptors = []
    for entry in FastaIterator.FastaIterator( IOTools.openFile( contaminant_file ) ):
        adaptors.append( "-a %s" % entry.sequence )
        
    adaptors= " ".join(adaptors)
    to_cluster = USECLUSTER

    statement = '''
    cutadapt 
    --discard
    %(adaptors)s
    --overlap=%(contamination_min_overlap_length)i
    --format=fastq
    %(contamination_options)s
    <( zcat < %(infile)s )
    2> %(outfile)s.log
    | gzip > %(outfile)s
    '''
    P.run()

@follows() 
def publish():
    '''publish files.'''
    P.publish_report( prefix = PARAMS["publish_prefix"] )

#########################################################################
#########################################################################
#########################################################################
@follows( runFastqc )
def full(): pass

@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    E.info( "starting documentation build process from scratch" )
    P.run_report( clean = True )

@follows( mkdir( "report" ) )
def update_report():
    '''update report.'''

    E.info( "updating documentation" )
    P.run_report( clean = False )

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )


