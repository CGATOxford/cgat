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

The pipeline can also be used to pre-process reads (target ``process_reads``). 

Implemented tasks are:

   * :meth:`removeContaminants` - remove contaminants from read sets
   * :meth:`trim` - trim reads by a certain amount
   * :meth:`filter` - filter reads by quality score

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
import Fastq
import csv2db
import cStringIO

USECLUSTER = True

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import Pipeline as P
P.getParameters( 
    ["%s/pipeline.ini" % __file__[:-len(".py")],
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
    m = PipelineMapping.FastQc()
    statement = m.build((infiles,), outfile) 
    P.run()

#########################################################################
#########################################################################
#########################################################################
## 
#########################################################################
@transform( runFastqc, suffix(".fastqc"), "_fastqc.load" )
def loadFastqc( infile, outfile ):
    '''load FASTQC stats.'''
    
    track = P.snip( infile, ".fastqc" )

    def section_iterator( infile ):

        data = []
        for line in infile:
            if line.startswith( ">>END_MODULE" ): 
                yield name, status, header, data
            elif line.startswith(">>"):
                name, status = line[2:-1].split("\t")
                data = []
            elif line.startswith("#"):
                header = "\t".join([ x for x in line[1:-1].split("\t") if x != ""] )
            else:
                data.append( "\t".join([ x for x in line[:-1].split("\t") if x != ""] ) )

    filename = os.path.join( PARAMS["exportdir"], "fastqc", track + "*_fastqc", "fastqc_data.txt" )

    for fn in glob.glob( filename ):
        prefix = os.path.basename( os.path.dirname( fn ) )
        results = []
        
        for name, status, header, data in section_iterator(IOTools.openFile( fn )):

            parser = csv2db.buildParser()
            (options, args) = parser.parse_args([])
            options.tablename = prefix + "_" + re.sub(" ", "_", name ) 
            options.allow_empty= True

            inf = cStringIO.StringIO( "\n".join( [header] + data ) + "\n" )
            csv2db.run( inf, options )
            results.append( (name, status ) )

        # load status table
        parser = csv2db.buildParser()
        (options, args) = parser.parse_args([])
        options.tablename = prefix + "_status"
        options.allow_empty= True

        inf = cStringIO.StringIO( "\n".join( ["name\tstatus"] + ["\t".join( x ) for x in results ] ) + "\n" )
        csv2db.run( inf, options )

    P.touch( outfile )

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

#########################################################################
#########################################################################
#########################################################################
@merge( None, "contaminants.fasta" )
def outputContaminants( infile, outfile ):
    '''output file with contaminants.'''
    outf = IOTools.openFile( outfile, "w")
    for key, value in ILLUMINA_ADAPTORS.iteritems():
        outf.write(">%s\n%s\n" % (key, value) )
    outf.close()

#########################################################################
#########################################################################
#########################################################################
@transform( [ x for x in \
                  glob.glob("*.fastq.gz") + glob.glob("*.fastq.1.gz") + glob.glob("*.fastq.2.gz") \
                  if not x.startswith("nocontaminants")],
	    regex( r"(\S+).(fastq.1.gz|fastq.gz|fastq.2.gz|csfasta.gz)"),
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


#########################################################################
#########################################################################
#########################################################################
@transform( [ x for x in \
                  glob.glob("*.fastq.gz") + glob.glob("*.fastq.1.gz") + glob.glob("*.fastq.2.gz") \
                  if not x.startswith("processed.")],
	    regex( r"(\S+).(fastq.1.gz|fastq.gz|fastq.2.gz|csfasta.gz)"),
	    add_inputs(outputContaminants),
	    r"processed.\1.\2")
def processReads( infiles, outfile ):
    '''process reads.'''
    
    infile, contaminant_file = infiles

    do_sth = False
    to_cluster = True

    # fastx does not like quality scores below 64 (Illumina 1.3 format)
    # need to detect the scores and convert
    
    format = Fastq.guessFormat( IOTools.openFile(infile ) , raises = False)
    E.info( "%s: format guess: %s" % (infile, format))
    offset = Fastq.getOffset( format, raises = False )

    if PARAMS["process_remove_contaminants"]:
        statement.append( '''
        cutadapt 
              --discard
              %(adaptors)s
              --overlap=%(contamination_min_overlap_length)i
              --format=fastq
              %(contamination_options)s
              <( zcat < %(infile)s )
              2> %(outfile)s_contaminants.log
        ''' )
        do_sth = True
    else:
        statement = ['zcat %(infile)s' ]

    if PARAMS["process_artifacts"]:
        statement.append( 'fastx_artifacts_filter -Q %(offset)i -v %(artifacts_options)s 2> %(outfile)s_artifacts.log' )
        do_sth = True
        
    if PARAMS["process_trim"]:
        statement.append( 'fastx_trimmer -Q %(offset)i -v %(trim_options)s 2> %(outfile)s_trim.log' )
        do_sth = True

    if PARAMS["process_filter"]:
        statement.append( 'fastq_quality_filter -Q %(offset)i -v %(filter_options)s 2> %(outfile)s_filter.log')
        do_sth = True

    if do_sth:
        statement.append( "gzip" )
        statement = " | ".join( statement ) + " > %(outfile)s" 
        P.run()
    else:
        E.warn( "no filtering specified for %s - nothing done" % infile )

#########################################################################
#########################################################################
#########################################################################
@merge( removeContaminants, "filtering.summary.tsv.gz" )
def summarizeFiltering( infiles, outfile ):
    '''collect summary output from filtering stage.'''

    tracks = {}
    adapters = {}
    def _chunker( inf ):
        chunk = []
        for line in inf:
            if line.startswith("==="):
                if chunk: yield chunk
                chunk = []
            chunk.append( line )
            
    for f in infiles:
        track = f[len("nocontaminants."):]
        track = re.sub( "[.].*", "", track )
        results = {}
        lines = IOTools.openFile( f + ".log" ).readlines()
        del lines[0]
        for x, line in enumerate(lines):
            if not line.strip(): continue
            if ":" in line:
                if line.strip().startswith("Command"): continue
                param, value = line[:-1].split(":")
                param = re.sub( " ", "_", param.strip()).lower()
                value = re.sub( "[a-zA-Z ].*", "", value.strip() )
                results[param] = value
            else:
                break
            
        del lines[:x]
        results["unchanged_reads"] = int(results["processed_reads"]) - int(results["trimmed_reads"])
        headers = results.keys()
        tracks[track] = results
        
        results = {}
        for chunk in _chunker(lines):        
            adapter = re.search("=== (.*) ===", chunk[0]).groups()[0]
            length, removed = re.search( "Adapter '.*', length (\d+), was trimmed (\d+) times", chunk[2]).groups()
                
            results[adapter] = length, removed
        
        adapters[track] = results

    outf = IOTools.openFile( outfile, "w" )
    outf.write( "track\t%s\n" % "\t".join(headers))
    
    for track, results in tracks.iteritems():
        outf.write("%s\t%s\n" % (track, "\t".join( str(results[x]) for x in headers ) ) )
    outf.close()

@transform( summarizeFiltering,
            suffix(".summary.tsv.gz"),
            "_summary.load")
def loadFilteringSummary( infile, outfile ):
    '''load filtering summary.'''
    P.load(infile, outfile )


#########################################################################
#########################################################################
#########################################################################
@transform( [ x for x in glob.glob("*.fastq.gz") + glob.glob("*.fastq.1.gz") + glob.glob("*.fastq.2.gz")]
            , regex( r"(\S+).(fastq.1.gz|fastq.gz|fastq.2.gz|csfasta.gz)"), r"trim.\1.\2")
def trimReads (infile, outfile):
    '''trim reads to desired length using fastx

    '''

    E.warn( "deprecated - use processReads instead" )

    to_cluster = True
    statement = '''zcat %(infile)s | fastx_trimmer %(trim_options)s 2> %(outfile)s.log | gzip > %(outfile)s''' 
    P.run()

#########################################################################
#########################################################################
#########################################################################

@transform( [ x for x in glob.glob("*.fastq.gz") + glob.glob("*.fastq.1.gz") + glob.glob("*.fastq.2.gz")]
            , regex( r"(\S+).(fastq.1.gz|fastq.gz|fastq.2.gz|csfasta.gz)"), r"replaced.\1.\2")
def replaceBaseWithN(infile, outfile):
    '''replaces the specified base with N'''

    to_cluster = True
    statement = '''python %(scriptsdir)s/fastq2N.py -i %(infile)s %(replace_options)s'''
    P.run()
    
#########################################################################
#########################################################################
#########################################################################
@follows() 
def publish():
    '''publish files.'''
    P.publish_report()

#########################################################################
#########################################################################
#########################################################################
@follows( loadFastqc )
def full(): pass

@follows( loadFilteringSummary )
def cleanData(): pass

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


