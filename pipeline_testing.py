################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_snps.py 2870 2010-03-03 10:20:29Z andreas $
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
=================================================
pipeline_testing - automated testing of pipelines
=================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

This pipeline executes other pipelines for testing purposes.

Overview
========

This pipeline implements automated testing of CGAT pipelines. This pipeline is run
regularly on the code in the master repository in order to check if all pipelines
are functional.

This pipeline can also be run from within a cloned repository before checking in 
changes. This is highly recommended, as changes in shared scripts might break existing
workflows.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

In order to run all tests, simply enter an empty directory and type::

   python <srcdir>/pipeline_testing.py make full
   python <srcdir>/pipeline_testing.py make build_report

The first command will run all pipelines while the second will build a summary
report.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.

Input
-----

As input, the pipeline requires a directory with test cases. The directory is
specified in the ``data_dir`` configuration variable. Each data directory needs
to start with the name of the pipeline that is to be tested. Valid names for test
cases are thus ``pipeline_annotations``, ``pipeline_annotations_lizard``, 
``pipeline_rnaseq_illumina``, ``pipeline_rnaseq_solid``, etc.

A test case directory contains all the data that is required for running a test within 
a directory, such as the genomic sequence, various data files and configuration
files. For example, for testing :doc:`pipeline_rnaseq`, the following files are given:

   data files
      Brain-F1-R1.fastq.gz,Brain-F1-R2.fastq.gz,Brain-F2-R1.fastq.gz,Brain-F2-R2.fastq.gz,UHR-F1-R1.fastq.gz,UHR-F1-R2.fastq.gz
   
   genomic sequence
      hg19.fasta,hg19.idx

   bowtie indices
      hg19.1.ebwt,hg19.2.ebwt,hg19.3.ebwt,hg19.4.ebwt,hg19.fa,hg19.fa.fai,hg19.rev.1.ebwt,hg19.rev.2.ebwt
      hg19_cs.1.ebwt,hg19_cs.2.ebwt,hg19_cs.3.ebwt,hg19_cs.4.ebwt,hg19_cs.rev.1.ebwt,hg19_cs.rev.2.ebwt

   configuration files
      conf.py,pipeline.ini,sphinxreport.ini

Dependencies between pipelines are currently implemented within the script itself. The 
dependencies that are currently implemented are:

* :doc:`pipeline_annotations` is run before any other pipelines.

Optional inputs
+++++++++++++++

Requirements
------------

All the test cases need to be properly configured. In particular, all the required software
for all the pipelines need to be in place in order for all tests to pass.

Pipeline output
===============

The major output is in the database file :file:`csvdb`.

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

def splitTestName( outfile ):
    '''split the test name into pipeline and data part.

    The name of a test starts with the name of a pipeline followed by some description.
    '''
    
    pipeline_components = outfile.split("_")

    g = os.path.join( os.path.dirname( __file__ ),  "pipeline_*.py")

    pipelines = set([ os.path.basename(x) for x in glob.glob( g ) ] )
    
    p = []
    for x in range(len(pipeline_components)):
        p.append( pipeline_components[x] )
        if ("_".join(p) + ".py") in pipelines:
            pipeline_name = "_".join( p )
            test_description = "_".join( pipeline_components[x+1:] )
            break
    else:
        raise ValueError("pipeline for test %s not found" % outfile )

    return pipeline_name, test_description 

###################################################################
###################################################################
###################################################################
def runTest( infile, outfile, update = False ):
    '''run a test.'''

    test_name = P.snip(outfile, ".log")

    pipeline_name, test_description = splitTestName( test_name )
    
    if not update:
        E.info("building test directory" )
        try:
            shutil.rmtree( "%s.dir" % test_name )
        except OSError: pass
        os.mkdir( "%s.dir" % test_name )
        statement = '''ln -s %(infile)s/* %(test_name)s.dir'''
        P.run()
    
    statement = '''
	    (cd %(test_name)s.dir; python %(scriptsdir)s/%(pipeline_name)s.py 
                      %(pipeline_options)s make full) >& %(outfile)s
        ''' 
    
    P.run()

###################################################################
###################################################################
###################################################################
## general tests
###################################################################
@files( [ (os.path.join( PARAMS["data_dir"], x), x + ".log" ) for x in
             P.asList(PARAMS["prerequisites"]) ] )
def prepareTests( infile, outfile ):
    '''run pre-requisite pipelines.'''
    runTest( infile, outfile )

###################################################################
###################################################################
###################################################################
## run a test
###################################################################
@follows( prepareTests )
@files( [ (x, os.path.basename(x) + ".log" ) for x in \
              glob.glob( os.path.join( PARAMS["data_dir"], "pipeline_*")) ] )
def runTests( infile, outfile ):
    '''run a pipeline with test data.'''
    runTest( infile, outfile )

###################################################################
###################################################################
###################################################################
## update a test
###################################################################
@follows( prepareTests )
@files( [ (x, os.path.basename(x) + ".log" ) for x in \
              glob.glob( os.path.join( PARAMS["data_dir"], "pipeline_*")) ] )
def updateTests( infile, outfile ):
    '''run a pipeline with test data.'''
    runTest( infile, outfile, update = True )

###################################################################
###################################################################
###################################################################
## build reports
###################################################################
@transform( runTests, suffix(".log"), ".report" )
def runReports( infile, outfile ):
    '''run a pipeline report.'''
    
    test_name = P.snip(outfile, ".report") 

    pipeline_name, test_description = splitTestName( test_name )

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
@follows( runTests, runReports )
def full(): pass

###################################################################
###################################################################
###################################################################
## primary targets
###################################################################

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
