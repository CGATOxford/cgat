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

The pipeline can also be used to pre-process reads (target ``process``). 

Implemented tasks are:

   * :meth:`removeContaminants` - remove contaminants from read sets
   * :meth:`trim` - trim reads by a certain amount
   * :meth:`filter` - filter reads by quality score
   * :meth:`sample` - sample a certain proportion of reads

Individual tasks are enabled in the configuration file.

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

###################################################
###################################################
###################################################
# load modules
###################################################

# import ruffus
from ruffus import *

# import useful standard python modules
import sys, os, re, shutil, itertools, math, glob, time, gzip, collections, random
import cStringIO

import CGAT.Experiment as E
import logging as L
import CGAT.Database as Database
import sys
import os
import re
import shutil
import string
import itertools
import math
import glob
import time
import gzip
import collections
import random
import numpy
import sqlite3
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
import CGAT.FastaIterator as FastaIterator
import CGAT.Tophat as Tophat
import rpy2.robjects as ro
import CGATPipelines.PipelineGeneset as PipelineGeneset
import CGATPipelines.PipelineMapping as PipelineMapping
import CGAT.Stats as Stats
import CGATPipelines.PipelineTracks as PipelineTracks
import CGAT.Pipeline as P
import CGAT.Fastq as Fastq
import CGAT.CSV2DB as CSV2DB

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
P.getParameters( 
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini" ] )
PARAMS = P.PARAMS

#########################################################################
#########################################################################
#########################################################################
# define input files
INPUT_FORMATS = ("*.fastq.1.gz", "*.fastq.gz", "*.sra", "*.csfasta.gz")
REGEX_FORMATS = regex( r"(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz)")

#########################################################################
#########################################################################
#########################################################################
@follows(mkdir(PARAMS["exportdir"]), mkdir(os.path.join(PARAMS["exportdir"], "fastqc")) )
@transform( INPUT_FORMATS,
            REGEX_FORMATS,
            r"\1.fastqc")
def runFastqc(infiles, outfile):
    '''convert sra files to fastq and check mapping qualities are in solexa format. 
    Perform quality control checks on reads from .fastq files.'''
    to_cluster = True
    m = PipelineMapping.FastQc(nogroup = PARAMS["readqc_no_group"] )
    statement = m.build((infiles,), outfile) 
    P.run()

#########################################################################
#########################################################################
#########################################################################
## 
#########################################################################
@jobs_limit( 1, "db" )
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

            parser = CSV2DB.buildParser()
            (options, args) = parser.parse_args([])
            options.tablename = prefix + "_" + re.sub(" ", "_", name ) 
            options.allow_empty= True

            inf = cStringIO.StringIO( "\n".join( [header] + data ) + "\n" )
            CSV2DB.run( inf, options )
            results.append( (name, status ) )

        # load status table
        parser = CSV2DB.buildParser()
        (options, args) = parser.parse_args([])
        options.tablename = prefix + "_status"
        options.allow_empty= True

        inf = cStringIO.StringIO( "\n".join( ["name\tstatus"] + ["\t".join( x ) for x in results ] ) + "\n" )
        CSV2DB.run( inf, options )

    P.touch( outfile )

#########################################################################
#########################################################################
#########################################################################
## adaptor trimming
#########################################################################
# these are the adaptors and PCR primers for used in various illumina libarary preps
# see https://cgatwiki.anat.ox.ac.uk/xwiki/bin/view/CGAT/Illumina+Sequencing#HIlluminaAdaptors.html
# currently included are primers/adaptors from:
# TruSeq DNA HT and RNA HT Sample Prep Kits; TruSeq DNA v1/v2/LT RNA v1/v2/LT and ChIP Sample Prep Kits;
# Oligonucleotide Sequences for TruSeq Small RNA Sample Prep Kits; Oligonucleotide Sequences for Genomic DNA;
# Oligonucleotide Sequences for the v1 and v1.5 Small RNA Sample Prep Kits; 
# Paired End DNA Oligonucleotide Sequences; Script Seq Adaptors;
# Oligonucleotide Sequences for the Multiplexing Sample Prep Oligo Only Kits.
ILLUMINA_ADAPTORS = { "Genomic-DNA-Adaptor" : "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",
                      "Genomic/Paired-End/Oligo-Only-Adaptor" : "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
                      "Genomic/TruSeq-Universal/PE/OO/ScriptSeq-Adaptor/PCR-Primer" : "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
                      "Genomic-PCR-Primer" : "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT",
                      "Paired-End-Adaptor" : "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",
                      "Paired-End-PCR-Primer" : "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",
                      "TruSeq-HT-Adaptor-I5" : "AATGATACGGCGACCACCGAGATCTACACNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT",
                      "TruSeq-HT-Adaptor-I7" : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
                      "TruSeq-LT-Adaptor-I6" : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
                      "TruSeq-Small-RNA-Adaptor" : "TGGAATTCTCGGGTGCCAAGG",
                      "TruSeq-Small-RNA-RT-Primer" : "GCCTTGGCACCCGAGAATTCCA",
                      "TruSeq-Small-RNA-PCR-Primer" : "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA",
                      "TruSeq-Small-RNA-PCR-Primer-I6" : "CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
                      "Oligo-Only-Adaptor" : "GATCGGAAGAGCACACGTCT",
                      "Oligo-Only-PCR-Primer" : "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
                      "Oligo-Only-PCR-Primer-I7" : "CAAGCAGAAGACGGCATACGAGATNNNNNNNTGACTGGAGTTC",
                      "Small-RNA-v1-RT-Primer" : "CAAGCAGAAGACGGCATACGA",
                      "Small-RNA-PCR-Primer" : "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA",
                      "ScriptSeq-Adaptor-I6" : "CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
                      "SmartIIA": "AAGCAGTGGTATCAACGCAGAGTAC",
                      }

#########################################################################
#########################################################################
#########################################################################
@merge( None, "contaminants.fasta" )
def outputContaminants( infile, outfile ):
    '''output file with contaminants.
    if contamination_reverse_complement is selected, then the reverse
    complement of each sequence is also written to the outfile. 
    '''
    outf = IOTools.openFile( outfile, "w")
    for key, value in ILLUMINA_ADAPTORS.iteritems(): 
        outf.write(">%s\n%s\n" % (key, value) )
        if PARAMS["contamination_reverse_complement"]:
            key_rc = key + "_rc"
            value_rc = value[::-1]
            value_rc = value_rc.translate(string.maketrans("ACTGN", "TGACN"))
            outf.write(">%s\n%s\n" % (key_rc, value_rc) )
        else: 
            continue
    outf.close()

def listAdaptors(infile): 
    adaptors = []
    for entry in FastaIterator.FastaIterator( IOTools.openFile( infile) ):
        adaptors.append("%s %s" % (PARAMS["contamination_trim_type"], entry.sequence) )
    adaptors = " ".join(adaptors)

    return adaptors

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

<<<<<<< local
    adaptors = []
    for entry in FastaIterator.FastaIterator( IOTools.openFile( contaminant_file ) ):
        adaptors.append( "-a %s" % entry.sequence )
        
    adaptors= " ".join(adaptors)
    to_cluster = True

=======
    adaptors = listAdaptors(contaminant_file)
    to_cluster = USECLUSTER
#    %(contamination_trim_type)s
>>>>>>> other
    statement = '''
    cutadapt 
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
def checkPairs( infile ):
    '''check for paired read files'''
    if infile.endswith( ".fastq.1.gz"):
        infile2 = P.snip( infile, ".fastq.1.gz") + ".fastq.2.gz"
        assert os.path.exists( infile2 ), "second part of read pair (%s) missing" % infile2
    else:
        infile2 = None
        
    return infile2

#########################################################################
#########################################################################
#########################################################################
@transform( [ x for x in \
                  glob.glob("*.fastq.gz") + glob.glob("*.fastq.1.gz")
                  if not x.startswith("processed.")],
	    regex( r"(\S+).(fastq.1.gz|fastq.gz|csfasta.gz)"),
	    add_inputs(outputContaminants),
	    r"processed.\1.\2")
def processReads( infiles, outfile ):
    '''process reads.'''

    infile, contaminant_file = infiles

    do_sth = False
    to_cluster = True

    infile2 = checkPairs( infile )

    if infile2:
        track = P.snip( outfile, ".fastq.1.gz" )        
        outfile2 = P.snip( outfile, ".fastq.1.gz" ) + ".fastq.2.gz"
    else:
        track = P.snip( outfile, ".fastq.gz" )


    if PARAMS["process_sample"] and infile2:
        E.warn( "sampling can not be combined with other processing for paired ended reads")
        statement = '''zcat %(infile)s
        | python %(scriptsdir)s/fastq2fastq.py 
                                   --sample=%(sample_proportion)f 
                                   --pair=%(infile2)s 
                                   --outfile-pair=%(outfile2)s 
                                   --log=%(outfile)s_sample.log
        | gzip 
        > %(outfile)s
        '''

        P.run()
        return

    # fastx does not like quality scores below 64 (Illumina 1.3 format)
    # need to detect the scores and convert
    format = Fastq.guessFormat( IOTools.openFile(infile ) , raises = False)
    E.info( "%s: format guess: %s" % (infile, format))
    offset = Fastq.getOffset( format, raises = False )

    if PARAMS["process_remove_contaminants"]:
        adaptors = listAdaptors(contaminant_file)
#              %(contamination_trim_type)s
        s = [ '''
        cutadapt 
              %(adaptors)s
              --overlap=%(contamination_min_overlap_length)i
              --format=fastq
              %(contamination_options)s
              <( zcat < %(infile)s )
              2>> %(outfile)s_contaminants.log
        ''' ]
        do_sth = True
    else:
        s = ['zcat %(infile)s' ]

    if PARAMS["process_artifacts"]:
        s.append( 'fastx_artifacts_filter -Q %(offset)i -v %(artifacts_options)s 2>> %(outfile)s_artifacts.log' )
        do_sth = True

    if PARAMS["process_trim"]:
        s.append( 'fastx_trimmer -Q %(offset)i -v %(trim_options)s 2>> %(outfile)s_trim.log' )
        do_sth = True

    if PARAMS["process_filter"]:
        s.append( 'fastq_quality_filter -Q %(offset)i -v %(filter_options)s 2>> %(outfile)s_filter.log')
        do_sth = True

    if PARAMS["process_sample"]:
        s.append( 'python %(scriptsdir)s/fastq2fastq.py --sample=%(sample_proportion)f --log=%(outfile)s_sample.log' )

    if not do_sth:
        E.warn( "no filtering specified for %s - nothing done" % infile )
        return

    s.append( "gzip" )
    if not infile2:
        statement = " | ".join( s ) + " > %(outfile)s" 
        P.run()
    else:
        tmpfile = P.getTempFilename(".")
        tmpfile1 = tmpfile + ".fastq.1.gz"
        tmpfile2 = tmpfile + ".fastq.2.gz"

        E.warn( "processing first of pair")
        # first read pair
        statement = " | ".join( s ) + " > %(tmpfile1)s" 
        P.run()

        # second read pair        
        E.warn( "processing second of pair")
        infile = infile2
        statement = " | ".join( s ) + " > %(tmpfile2)s" 
        P.run()

        # reconcile
        E.info("starting reconciliation" )
        statement = """python %(scriptsdir)s/fastqs2fastq.py
                           --method=reconcile
                           --output-pattern=%(track)s.fastq.%%i.gz
                           %(tmpfile1)s %(tmpfile2)s
                     > %(outfile)s_reconcile.log"""
        
        P.run()

        os.unlink( tmpfile1 )
        os.unlink( tmpfile2 )
        os.unlink( tmpfile )

<<<<<<< local
#########################################################################
#########################################################################
#########################################################################
=======
##################################################################
##################################################################
##################################################################
def parseCutadapt( lines ):
    '''parse cutadapt output.

    Multiple cutadapt outputs are joined.
    '''

    def _chunker( inf ):
        chunk = []
        for line in inf:
            if line.startswith("==="):
                if chunk: yield chunk
                chunk = []
            chunk.append( line )
        
    assert lines[0].startswith("cutadapt")
    results = {}

    del lines[0]
    for x, line in enumerate(lines):
        if not line.strip(): continue
        if ":" in line:
            if line.strip().startswith("Command"): continue
            param, value = line[:-1].split(":")
            param = re.sub( " ", "_", param.strip()).lower()
            value = re.sub( "[a-zA-Z ].*", "", value.strip() )
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except:
                    pass
            results[param] = value
        else:
            break

    del lines[:x]
    results["unchanged_reads"] = int(results["processed_reads"]) - int(results["trimmed_reads"])
    headers = results.keys()
    
    adapters = {}
    for chunk in _chunker(lines):        
        adapter = re.search("=== (.*) ===", chunk[0]).groups()[0]
        length, removed = re.search( "Adapter '.*', length (\d+), was trimmed (\d+) times", chunk[2]).groups()

        adapters[adapter] = length, removed

    return results, adapters

##################################################################
##################################################################
##################################################################
>>>>>>> other
@transform( processReads,
            suffix(""),
            ".tsv")
def summarizeProcessing( infile, outfile ):
    '''build processing summary.'''

    def _parseLog( inf, step ):

        inputs, outputs = [], []
        if step == "reconcile":
            for line in inf:
                x = re.search( "first pair: (\d+) reads, second pair: (\d+) reads, shared: (\d+) reads", line )
                if x:
                    i1, i2, o = map(int, x.groups())
                    inputs = [i1,i2]
                    outputs = [o,o]
                    break
        elif step == "contaminants":
            lines = inf.readlines()
            assert lines[0].startswith("cutadapt")
            lines = "@@@".join( lines )
            for part in lines.split("cutadapt")[1:]:
                results, adapters = parseCutadapt( ("cutadapt" + part).split("@@@") )
                inputs.append( results["processed_reads"] )
                outputs.append( results["unchanged_reads" ] )
        else:
            for line in inf:
                if line.startswith( "Input:"):
                    inputs.append( int( re.match( "Input: (\d+) reads.", line).groups()[0] ) )
                elif line.startswith( "Output:"):
                    outputs.append( int( re.match( "Output: (\d+) reads.", line).groups()[0] ) )

        return zip(inputs, outputs)
    

    infile2 = checkPairs( infile )
    if infile2: 
        track = P.snip( infile, ".fastq.1.gz")        
    else:
        track = P.snip( infile, ".fastq.gz" )

    outf = IOTools.openFile( outfile, "w")
    outf.write( "track\tstep\tpair\tinput\toutput\n")

    for step in "contaminants", "artifacts", "trim", "filter", "reconcile":
        fn = infile + "_%s.log" % step
        if not os.path.exists(fn): continue
        for x, v in enumerate( _parseLog( IOTools.openFile( fn ), step)):
            outf.write( "%s\t%s\t%i\t%i\t%i\n" % (track, step, x, v[0], v[1]) )

    outf.close()

#########################################################################
#########################################################################
#########################################################################
@jobs_limit( 1, "db" )
@transform( summarizeProcessing,
            regex(r"processed.(\S+).fastq.*.gz.tsv"),
            r"\1_processed.load")
def loadProcessingSummary( infile, outfile ):
    '''load filtering summary.'''
    P.load(infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@merge( summarizeProcessing, "processing_summary.tsv" )
def summarizeAllProcessing( infiles, outfile ):
    '''summarize processing information.'''

    outf = IOTools.openFile( outfile, "w" )
    data = []
    for infile in infiles:
        inf = IOTools.openFile( infile )
        for line in inf:
            track, step, pair, ninput, noutput = line[:-1].split("\t")
            if track == "track": continue
            data.append( (track, step, pair, ninput, noutput) )
            
    # sort by track, pair, input
    data.sort( key = lambda x: (x[0], x[2], -int(x[3])))
    first = True
    for key, v in itertools.groupby( data, lambda x: (x[0], x[2])):
        vals = list(v)
        track,pair = key
        ninput = int(vals[0][3])
        outputs = [int(x[4]) for x in vals]
        if first:
            outf.write( "track\tpair\tninput\t%s\t%s\t%s\t%s\n" % ("\t".join( [x[1] for x in vals] ),
                                                                   "noutput",
                                                                   "\t".join( ["percent_%s" % x[1] for x in vals] ),
                                                                   "percent_output" ))
            first = False
        outf.write( "%s\t%s\t%i\t%s\t%i\t%s\t%s\n" % ( track, pair, ninput, 
                                                       "\t".join( map(str,outputs)),
                                                       outputs[-1], 
                                                       "\t".join( [ "%5.2f" % (100.0 * x / ninput) for x in outputs ] ),
                                                       "%5.2f" % (100.0 * outputs[-1] / ninput)))
    outf.close()

#########################################################################
#########################################################################
#########################################################################
@jobs_limit( 1, "db" )
@transform( summarizeAllProcessing, suffix(".tsv"), ".load" )
def loadAllProcessingSummary( infile, outfile ):
    P.load( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@merge( removeContaminants, "filtering.summary.tsv.gz" )
def summarizeFiltering( infiles, outfile ):
    '''collect summary output from filtering stage.'''

    tracks = {}
    adapters = {}

    for f in infiles:
        track = f[len("nocontaminants."):]
        track = re.sub( "[.].*", "", track )
        result, adapter = parseCutadapt( IOTools.openFile( f + ".log" ) )
        tracks[track] = result
        adapters[track] = adapter
        header = result.keys()

    outf = IOTools.openFile( outfile, "w" )
    outf.write( "track\t%s\n" % "\t".join(headers))
    
    for track, results in tracks.iteritems():
        outf.write("%s\t%s\n" % (track, "\t".join( str(results[x]) for x in headers ) ) )
    outf.close()

<<<<<<< local
#########################################################################
#########################################################################
#########################################################################
=======
##################################################################
<<<<<<< local
>>>>>>> other
=======
@jobs_limit( 1, "db" )
>>>>>>> other
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
#########################################################################
@follows( loadProcessingSummary, loadAllProcessingSummary )
def process():
    '''process (filter,trim) reads.'''
    pass

#########################################################################
@follows( loadFastqc )
def full(): pass

#########################################################################
@follows( loadFilteringSummary )
def cleanData(): pass

#########################################################################
@follows() 
def publish():
    '''publish files.'''
    P.publish_report()

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


