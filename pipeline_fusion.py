################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_fusion.py 0001 2011-11-23 10:20:29Z Ian $
#
#   Copyright (C) 2009 Andreas Heger
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
===========================
Pipeline template
===========================

:Author: Ian Subery
:Release: $Id$
:Date: |today|
:Tags: Python

A pipeline to run various packages for detecting fusion transcripts in RNA-seq data.

Overview
========

The pipleline takes a set up reads in any legal format for tophat and performs the following tasks:

    * Runs tophat-fusion on each track seperately collecting the data in directiories of the form:
     track_tophat
    * Runs tophat-fusion-post on all track simulateously to build a single report of fusions found in the samples

Eventaully, this pipeline may include other algorithms for finding fusion transcripts and tasks for  QCing the
resulting calls. Currently the pipline uses the annotation provided by tophat  fusion. Eventaully we hope to be
able to use our own annontations.

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

Reads
+++++

As per the rnaseq pipeline.

Reads are imported by placing files are linking to files in the :term:`working directory`.

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
Optional inputs
+++++++++++++++

Requirements
------------

The pipeline requires the results from :doc:`pipeline_annotations`. Set the configuration variable 
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
| blast              | 2.2.25+           | Used by tophat-fusion-post                     |                       
+--------------------+-------------------+------------------------------------------------+
| blastall script    | 0.1               | Provides an alias to the now depcrecated       |
|                    |                   | blastall for those with Blast+                 |
+--------------------+-------------------+------------------------------------------------+
| tophat-fusion-post | 0.1.0             | post-postprocessing of tophat-fusion results   |
+--------------------+-------------------+------------------------------------------------+
| tophat-fusion      | 0.1.0             | read mapping                                   |
+--------------------+-------------------+------------------------------------------------+

In addition, the blast databases need to be availible and their location specified in the ini.

Pipeline output
===============

Currently the major output if in the tophatfusion_out directory. Integration with database coming soon.

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_template.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_template.tgz
   tar -xvzf pipeline_template.tgz
   cd pipeline_template
   python <srcdir>/pipeline_template.py make full

.. note:: 
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *

import sys, glob, gzip, os, itertools, re, math, types, collections, time
import optparse, shutil
import sqlite3
import Experiment as E
import IOTools
import Database
import PipelineMapping

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
PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
                                       "pipeline_annotations.py" )

USECLUSTER=True
###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
import PipelineTracks

###################################################################
###################################################################
###################################################################
def connect():
    '''connect to database.

    Use this method to connect to additional databases.

    Returns a database connection.
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
## worker tasks
###################################################################
@files(None,("ensGene.txt","ensGtp.txt","mcl","refGene_sorted.txt","blast_human"))
def prepare_directory(infiles,outfiles):
    ''' Prepares the directory neccesary for the pipeline. Links 
     the required annotation files and blast databases '''
    
    tophat_fusion_path = PARAMS['tophatfusion_installpath']
    annotations_path = os.path.join(tophat_fusion_path,'bin','annotation')

    for file in os.listdir(annotations_path):
        if not os.path.exists(file):
            E.info('Linking file %(annotations_path)s/%(file)s -> %(file)s' % locals())
            os.symlink('%(annotations_path)s/%(file)s' % locals(), file)
        else:
            E.info('File %s already exists, no need to link' % file)

    if not os.path.exists('blast_human'):
        E.info('Linking blast_db from %s -> blast_human' % PARAMS['blast_dbpath'])
        os.symlink(PARAMS['blast_dbpath'],'blast_human')
    else:
        E.info('blast_human directory already exists, no need to link')

##################################################################
##################################################################


@transform(("*.fastq.1.gz", 
             "*.fastq.gz",
             "*.sra",
             "*.csfasta.gz",
             "*.csfasta.F3.gz",
             ),
            regex( r"(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz|csfasta.F3.gz)"), 
            r"tophat_\1/accepted_hits.sam" )
def mapReadsWithTophatFusion( infiles, outfile ):
    '''map reads from .fastq or .sra files and find candidate fusions

    A list with known splice junctions expect from rnaseq pipeline
    '''
    
    job_options= "-pe dedicated %i -R y" % PARAMS["tophat_threads"]
    
    if "--butterfly-search" in PARAMS["tophat_options"]:
        # for butterfly search - require insane amount of
        # RAM.
        job_options += " -l mem_free=50G"

    to_cluster = USECLUSTER
    m = PipelineMapping.TopHat_fusion()
    infile = infiles

    #if a file of reference junctions, as generated by the rnaseq pipline,
    #has been specified in the ini, then pass this to tophat-fusion
    if not PARAMS['tophatfusion_reference_junctions']==None:
        reffile = PARAMS['tophatfusion_reference_junctions']
        tophat_options = PARAMS["tophat_options"] + " --raw-juncs %(reffile)s" % locals()

    tophatfusion_options = PARAMS["tophatfusion_options"]
    statement = m.build( (infile,), outfile ) 
    P.run()


###################################################################
###################################################################

@follows(prepare_directory)
@merge(mapReadsWithTophatFusion,'export/tophatfusion_out/result.html')
def postprocessTopHatFusion(infiles,outfile):
    ''' Uses tophat-fusion-post to postprocess and filter all of the
        tophat-fusion output into one report. Slow as it is not 
        cluster aware and spawns a large number of blast tasks'''

    to_cluster = USECLUSTER
    job_options = ' -pe dedicated %i -R y' % PARAMS["tophatfusion_postthreads"]
    statement = '''
                  module load tophatfusion;
                  tophat-fusion-post -p %(tophatfusion_postthreads)s
                                   %(tophatfusion_postoptions)s 
                                   %(bowtie_index_dir)s/%(genome)s
                  &> tophatfusion_out.log
                '''

    P.run()
    
    #put the results in the export directory.

    #if the export directory doesn't exist, create it
    if not os.path.exists('export'):
        os.mkdir('export')

    #otherwise if it does, then delete any out directory that is 
    #already there.
    elif os.path.exists('export/tophatfusion_out') and os.path.isdir('export/tophatfusion_out'):
        shutil.rmtree('export/tophatfusion.out')

    shutil.move('tophatfusion_out','export')

#################################################################
#################################################################

@follows(postprocessTopHatFusion)
@transform(mapReadsWithTophatFusion,
           suffix("accepted_hits.sam"),
           "junctions.bed.gz")
def doCleanUp(infile,outfile):
    ''' Does the clean-up that should have been done at the end of 
     the tophat run  but wasn't because the files were needed by 
     tophat-fusion-post '''

    to_cluster = USECLUSTER
    indir = P.snip(infile,'accepted_hits.sam')
    juncfile = indir + 'junctions.bed'
    statement = "gzip %(juncfile)s" % locals()
    P.run()


###################################################################
###################################################################
###################################################################
## primary targets
###################################################################

@follows( doCleanUp )
def full(): pass

@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    E.info( "starting report build process from scratch" )
    P.run_report( clean = True )

@follows( mkdir( "report" ) )
def update_report():
    '''update report.'''

    E.info( "updating report" )
    P.run_report( clean = False )

@follows( update_report )
def publish_report():
    '''publish report.'''

    E.info( "publishing report" )
    P.publish_report()

if __name__== "__main__":

    # P.checkFiles( ("genome.fasta", "genome.idx" ) )
    sys.exit( P.main(sys.argv) )
