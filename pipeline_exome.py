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
readqc pipeline
====================

:Author: David Sims
:Release: $Id$
:Date: |today|
:Tags: Python

The exome pipeline imports unmapped reads from one or more
fastq and aligns them to the genome, then filters to the regions of interest and calls variants
in those regions:

   1. align to genome using gapped alignment (BWA, Stampy)
   2. filter to target regions (BEDTools)
   3. Call SNVs (SAMTools, SNVMix)
   4. Call indels (SAMTools)
   5. Filter variants using Ensembl variation / 1000 genomes data
   6. Annotate functional impact of remaining variants

For further details see http://www.bioinformatics.bbsrc.ac.uk/projects/exome/

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

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
|Stampy              |>=0.9.0            |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|BWA                 |                   |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|SAMtools            |                   |filtering, SNV / indel calling                  |
+--------------------+-------------------+------------------------------------------------+
|BEDTools            |                   |filtering, SNV / indel calling                  |
+--------------------+-------------------+------------------------------------------------+
|sra-tools           |                   |extracting reads from .sra files                |
+--------------------+-------------------+------------------------------------------------+
|picard              |>=1.38             |bam/sam files. The .jar files need to be in your|
|                    |                   | CLASSPATH environment variable.                |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The major output is a set of categorised, prioritised variant calls with functional annotations.

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/rnaseq.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_rnaseq.tgz
   tar -xvzf pipeline_rnaseq.tgz
   cd pipeline_rnaseq
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
import GTF, IOTools, IndexedFasta
import Tophat
import rpy2.robjects as ro
import PipelineGeneset
import PipelineMapping
import Stats


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

USECLUSTER = True

###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
import PipelineTracks

# collect sra nd fastq.gz tracks
#TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
#    glob.glob( "*.sra" ), "(\S+).sra" ) +\
#    PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
#    glob.glob( "*.fastq.gz" ), "(\S+).fastq.gz" ) +\
#    PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
#    glob.glob( "*.fastq.1.gz" ), "(\S+).fastq.1.gz" )

#ALL = PipelineTracks.Sample3()
#EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", "tissue" ) )
#CONDITIONS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", ) )
#TISSUES = PipelineTracks.Aggregate( TRACKS, labels = ("tissue", ) )


#########################################################################
#########################################################################
#########################################################################
#########################################################################
@follows(mkdir("bam"))
@transform( ("*.fastq.1.gz", 
              "*.fastq.gz",
              "*.sra"),
              regex( r"(\S+).(fastq.1.gz|fastq.gz|sra)"),
              r"bam/\1.bam")
def mapReads(infiles, outfile):
        '''Map reads to the genome using BWA '''
        to_cluster = USECLUSTER
        m = PipelineMapping.bwa()
        statement = m.build((infiles,), outfile) 
        print statement
        P.run()


#########################################################################
#########################################################################
#########################################################################
#########################################################################
@follows(mkdir("variants"))
@transform( "*.bam",
              regex( r"(\S+).bam"),
              r"variants/\1.bcf")
def runSAMtools(infiles, outfile):
        '''Perform SNV and indel called from gapped alignment using SAMtools '''
        to_cluster = USECLUSTER
        statement = []
        statement.append('''mkdir -p %(outfile)s;''' % locals())
        gd = PARAMS["genome_dir"]
        #print type(infiles)
        g = PARAMS["genome"]
        track = P.snip( os.path.basename( infiles ), ".bam" )
        statement.append('''samtools mpileup -ugf %(gd)s/%(g)s.fa %(infiles)s | bcftools view -bvcg - > %(outfile)s; ''' % locals() )
        statement = " ".join( statement )
        print statement
        P.run()

#########################################################################
#########################################################################
#########################################################################
@follows( runSAMtools )
def full(): pass


if __name__== "__main__":
    sys.exit( P.main(sys.argv) )


