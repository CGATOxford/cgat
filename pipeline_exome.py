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
     "../exome.ini",
     "exome.ini" ] )
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
@transform( ("*.fastq.1.gz", 
              "*.fastq.gz",
              "*.sra"),
              regex( r"(\S+).(fastq.1.gz|fastq.gz|sra)"),
              r"\1/bam/\1.bam")
def mapReads(infiles, outfile):
        '''Map reads to the genome using BWA '''
        to_cluster = USECLUSTER
        track = P.snip( os.path.basename(outfile), ".bam" )
        try: os.mkdir( track )
        except OSError: pass
        try: os.mkdir( '''%(track)s/bam''' % locals() )
        except OSError: pass
        m = PipelineMapping.bwa()
        statement = m.build((infiles,), outfile) 
        #print statement
        P.run()

#########################################################################
#########################################################################
#########################################################################
#########################################################################
@transform( mapReads,
              regex( r"(\S+)/bam/(\S+).bam"),
              r"\1/bam/\1.roi.bam")
def filterBamROI(infiles, outfile):
        '''Filter alignments in BAM format to regions of interest from a bed file.
           Todo: use multiple BAM files'''
        to_cluster = USECLUSTER
        statement = '''intersectBed -u -abam %(infiles)s -b %%(roi_bed)s > %(outfile)s; ''' % locals()
        statement += '''samtools index %(outfile)s; ''' % locals()
        #print statement
        P.run()


#########################################################################
#########################################################################
#########################################################################
#########################################################################
@transform( (mapReads, filterBamROI), 
            suffix(".bam"),
            ".readstats" )
def buildBAMStats( infile, outfile ):
    '''count number of reads mapped, duplicates, etc. '''

    to_cluster = USECLUSTER
    scriptsdir = PARAMS["general_scriptsdir"]

    statement = '''python %(scriptsdir)s/bam2stats.py --force 
                   --output-filename-pattern=%(outfile)s.%%s < %(infile)s > %(outfile)s'''

    P.run()


#########################################################################
#########################################################################
#########################################################################
#########################################################################
@merge( buildBAMStats, "bam_stats.load" )
def loadBAMStats( infiles, outfile ):
    '''import bam statisticis.'''

    scriptsdir = PARAMS["general_scriptsdir"]
    header = ",".join( [P.snip( os.path.basename(x), ".readstats") for x in infiles] )
    filenames = " ".join( [ "<( cut -f 1,2 < %s)" % x for x in infiles ] )
    tablename = P.toTable( outfile )
    E.info( "loading bam stats - summary" )
    statement = """python %(scriptsdir)s/combine_tables.py
                      --headers=%(header)s
                      --missing=0
                      --ignore-empty
                   %(filenames)s
                | perl -p -e "s/bin/track/"
                | perl -p -e "s/unique/unique_alignments/"
                | python %(scriptsdir)s/table2table.py --transpose
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s
            """
    P.run()

    for suffix in ("nm", "nh"):
        E.info( "loading bam stats - %s" % suffix )
        filenames = " ".join( [ "%s.%s" % (x, suffix) for x in infiles ] )
        tname = "%s_%s" % (tablename, suffix)
        
        statement = """python %(scriptsdir)s/combine_tables.py
                      --header=%(header)s
                      --skip-titles
                      --missing=0
                      --ignore-empty
                   %(filenames)s
                | perl -p -e "s/bin/%(suffix)s/"
                | python %(scriptsdir)s/csv2db.py
                      --table=%(tname)s 
                      --allow-empty
                >> %(outfile)s
                """
        P.run()


#########################################################################
#########################################################################
#########################################################################
#########################################################################
@transform((mapReads, filterBamROI),
              regex( r"(\S+)/bam/(\S+).bam"),
              r"\1/bam/\2.cov.bamstats")
def coverageStats(infiles, outfile):
        '''Generate coverage statistics for regions of interest from a bed file using BAMstats'''
        to_cluster = USECLUSTER
        statement = '''bamstats -i %(infiles)s -o %(outfile)s -f %%(roi_bed)s; ''' % locals()
        #print statement
        P.run()


#########################################################################
#########################################################################
#########################################################################
#########################################################################
@transform((mapReads, filterBamROI),
              regex( r"(\S+)/bam/(\S+).bam"),
              r"\1/bam/\2.cov.bedtools")
def coverageStatsBedtools(infiles, outfile):
        '''Generate coverage statistics for regions of interest from a bed file using bedtools'''
        to_cluster = USECLUSTER
        statement = '''coverageBed -abam %(infiles)s -b %%(roi_bed)s -hist > %(outfile)s;''' % locals()
        #print statement
        P.run()


#########################################################################
#########################################################################
#########################################################################
#########################################################################
@merge( coverageStats, "coverage_stats.load" )
def loadCoverageStats( infiles, outfile ):
    '''import coverage statistics.'''

    scriptsdir = PARAMS["general_scriptsdir"]
    header = ",".join( [P.snip( os.path.basename(x), ".cov.bamstats") for x in infiles] )
    #filenames = " ".join( [ "<( cut -f 1,2 < %s)" % x for x in infiles ] )
    tablename = P.toTable( outfile )
    E.info( "loading coverage stats - summary" )
    statement = """python %(scriptsdir)s/combine_tables.py
                      --headers=%(header)s
                      --missing=0
                      --ignore-empty
                   %(infiles)s
                | perl -p -e "s/bin/track/"
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s
            """
    P.run()


#########################################################################
#########################################################################
#########################################################################
#########################################################################
@transform(   (mapReads, filterBamROI),
              regex( r"(\S+)/bam/(\S+).bam"),
              r"\1/variants/\2.vcf")
def callVariantsSAMtools(infiles, outfile):
        '''Perform SNV and indel called from gapped alignment using SAMtools '''
        to_cluster = USECLUSTER
        statement = []
        outfolder = infiles[:infiles.find("/")]
        try: os.mkdir( '''%(outfolder)s/variants''' % locals() )
        except OSError: pass
        track = P.snip( os.path.basename(infiles), ".bam" )
        statement.append('''samtools mpileup -ugf %%(genome_dir)s/%%(genome)s.fa %(infiles)s | bcftools view -bvcg - > %(outfolder)s/variants/%(track)s.bcf; ''' % locals() )
        statement.append('''bcftools view %(outfolder)s/variants/%(track)s.bcf | vcfutils.pl varFilter -d 10 -D 100 > %(outfile)s; ''' % locals() )
        statement.append('''vcf-stats %(outfile)s > %(outfile)s.stats;''' % locals() )
        statement = " ".join( statement )
        #print statement
        P.run()

#########################################################################
#########################################################################
#########################################################################
#########################################################################
@transform( callVariantsSAMtools,
              regex( r"(\S+)/variants/(\S+).vcf"),
              r"\1/roi_variants/\2.vcf")
def filterVariantROI(infiles, outfile):
        '''Filter variant calls in vcf format to regions of interest from a bed file'''
        to_cluster = USECLUSTER
        outfolder = infiles[:infiles.find("/")]
        try: os.mkdir( '''%(outfolder)s/roi_variants''' % locals() )
        except OSError: pass
        statement = '''intersectBed -u -a %(infiles)s -b %%(roi_bed)s > %(outfile)s; ''' % locals()
        #print statement
        P.run()

#########################################################################
#########################################################################
#########################################################################
@follows( mapReads,
          filterBamROI,
          buildBAMStats,
          loadBAMStats,
          coverageStats,
          callVariantsSAMtools,
          filterVariantROI  )
def full(): pass


if __name__== "__main__":
    sys.exit( P.main(sys.argv) )


