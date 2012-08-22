################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_cpg.py 2900 2011-05-24 14:38:00Z david $
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
======================
Read trimming pipeline
======================

:Author: David Sims 
:Release: $Id: trim_reads_fastq.py 2900 2011-05-24 14:38:00Z david $
:Date: |today|
:Tags: Python

This simple pipeline uses the fastx toolkit to trim reads to a specified length for a seq of fastq files.


Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline_cpg.ini` file. The pipeline looks for a configuration file in several places:

   1. The default configuration in the :term:`code directory`.
   2. A shared configuration file :file:`../pipeline.ini`.
   3. A local configuration :file:`pipeline.ini`.

The order is as above. Thus, a local configuration setting will
override a shared configuration setting and a default configuration
setting.

Configuration files follow the ini format (see the python
`ConfigParser <http://docs.python.org/library/configparser.html>` documentation).
The configuration file is organized by section and the variables are documented within 
the file. In order to get a local configuration file in the current directory, type::

    python <codedir>/pipeline_cpg.py config

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.


Input
-----

Reads
++++++

Input are :file:`.fastq.gz`-formatted files. 

Requirements
------------

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|FastX toolkit       |?                  |read trimming                                   |
+--------------------+-------------------+------------------------------------------------+


Pipline Output
==============

The results of the computation are stored in .trim.fastq files

Code
====

"""
import sys, tempfile, optparse, shutil, itertools, csv, math, random, re, glob, os, shutil, collections, gzip


import Experiment as E
import logging as L
from ruffus import *
import PipelineMapping

USECLUSTER = True

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################
import Pipeline as P
P.getParameters(  ["%s.ini" % __file__[:-len(".py")],  "../pipeline.ini", "pipeline.ini" ] )
PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],"pipeline_annotations.py" )


###################################################################
###################################################################
###################################################################
## TRIM READS
@follows(mkdir("trim"))
@transform( "*.gz", regex( r"(\S+).gz"), r"trim/\1.gz" )
def trimReads( infile, outfile ):
    '''trim reads with FastX'''
    to_cluster = True

    tmpdir_fastq = P.getTempDir()
    track = P.snip( os.path.basename( infile ), ".gz" )
    statement = """gunzip < %(infile)s | python %%(scriptsdir)s/fastq2fastq.py 
                       --change-format=sanger 
                       --guess-format=phred64 
                       --log=%(outfile)s.log
                   > %(tmpdir_fastq)s/%(track)s;""" % locals()
    statement += """zcat %(infile)s | fastx_trimmer -f %(first_base)s -l %(last_base)s -z -o %(outfile)s """

    P.run()


