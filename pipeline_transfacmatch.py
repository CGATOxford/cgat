################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_chipseq.py 2900 2010-04-13 14:38:00Z andreas $
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
#   along with this program; if not, write to the Free Software#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
"""

========================
Transfac match pipeline
========================

:Author: Nick Ilott
:Release: $Id: pipeline_transfacmatch 
:Date: |today|
:Tags: Python

The transfac match pipeline takes a several set of :term:`bed` or :term:`gtf` formatted files 
and scans intervals for transcription factor binding motifs.

It performs the following analyses:
   * transfac match motif analysis
      * If bed files are submitted then intervals must be specified with a unique name id
      * If gtf files are submitted then specification of interval locations should be provided 
   * enrichment of motifs in a foreground set of ids vs. background set of ids
      * foreground and background sets must be specified

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. The pipeline looks for a configuration file in several places:

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

    python <codedir>/pipeline_transfacmatch.py config

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.

Input
-----

Intervals / GTF entries
++++++++++++++++++++++++

Input are :term:`bed`-formatted files of intervals or :term:`gtf` formatted files. If intervals are 
submiited then they should be at least bed4 formatted, i.e., each interval should be labelled (uniquely).


Background and foreground sets
+++++++++++++++++++++++++++++++

Background and foreground sets are required to test for enrichments of TF motifs. These
are therefore subsets of the bed or gtf file. It consists of a :term:`tsv` file:

+---------+
|   id    |
+---------+
|  id1    |
|  id2    |
+---------+

Background and foreground files must take the form <name>.background.tsv <name>.foreground.tsv, repectively.

Motif library
++++++++++++++

The pipeline requires that a transfac matrix library be present. Its location is specified in the
pipeline.ini. 


Motif profiles
+++++++++++++++

The pipeline requires a file that contains a selection of matrices to search with defined cutoffs. 
Transfac match will search only matrices that are specified in profile. 
It is advised that all matrices in the transfac profile are used.

It will return only hits with scores that are
equal or higher than specified in profiles.


Requirements
------------

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|(transfac) match    |                   |sequence scanning for TF motifs                 |
+--------------------+-------------------+------------------------------------------------+

Pipline Output
==============

The results of the computation are all stored in an sqlite relational
database :file:`csvdb`.


"""


# load modules
from ruffus import *
import Experiment as E
import logging as L
import Database, CSV
import numpy as np
import fnmatch
import sqlite3
import list_overlap
import FastaIterator
import random
import GTF
import IOTools
import sys, os, re, shutil, itertools, math, glob, time, gzip, collections, random
from rpy2.robjects import r as R
from rpy2.robjects.packages import importr
from rpy2.robjects.numpy2ri import numpy2ri
import rpy2.robjects as robjects
import PipelineTracks
import PipelineTransfacMatch

#########################################################################
#########################################################################
#########################################################################
use_cluster = True

import Pipeline as P
P.getParameters( 
    ["%s.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini" ] )
PARAMS = P.PARAMS

################################################
# helper functions mapping input files
################################################

# input intervals / gtf entries
INPUT_FORMATS = ("*.gtf.gz", "*.bed.gz")
REGEX_FORMATS = regex(r"(\S+).(gtf.gz|bed.gz)")
for x in INPUT_FORMATS:
    input_file = glob.glob(x)
    if len( input_file ) != 0:
        INPUT_FILE = input_file[0]

# foreground and background sets
INPUT_FOREGROUND = glob.glob("*.foreground.tsv")
INPUT_BACKGROUND = glob.glob("*.background.tsv")

TRACK_FOREGROUND = [P.snip(x, ".foreground.tsv") for x in INPUT_FOREGROUND]
TRACK_BACKGROUND = [P.snip(x, ".background.tsv") for x in INPUT_BACKGROUND]

#########################################################################
#########################################################################
#########################################################################
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

# P.toTable()
#########################################################################
#########################################################################
#########################################################################
def filenameToTablename(filename):
    '''
    converts filename containing "." to tablename where "." converted to "_"
    '''
    return filename.replace(".", "_" )

# P.touch()
#########################################################################
#########################################################################
#########################################################################
def sentinelFile(filename):
    '''
    create empty file for updating purposes
    '''
    outf = open(filename, "w")
    outf.write("file created for ruffus update")
    outf.close()

#########################################################################
#########################################################################
#########################################################################
@follows(mkdir("fasta.dir"))
@transform( INPUT_FILE
           , REGEX_FORMATS
           , r"fasta.dir/\1.fasta")
def buildIntervalsFasta(infile, outfile):
    '''
    build fasta file from intervals. Alternatively
    if a gtf file is specified this function will
    use parameters specified in the .ini file to 
    use intervals upstream / downstream of tss 
    '''

    # define upstream and downstream extensions
    upstream = PARAMS["intervals_extension_upstream"]
    downstream = PARAMS["intervals_extension_downstream"]

    assert len(str(upstream)), """extension_upstream cannot be of %s type. 
                        If no extension is to be used specify 0""" % type(upstream)
    assert len(str(downstream)), """downstream extension cannot be of %s type. 
                        If no extension is to be used specify 0""" % type(downstream)

    # if input is gtf then convert to bed
    # with intervals defined by .ini file
    temp = P.getTempFile()
    if infile.endswith(".gtf.gz"):
        # the resulting temporary file will not be zipped
        concatenate = "cat"
        for gene in GTF.merged_gene_iterator(GTF.iterator(IOTools.openFile(infile))):
            if gene.strand == "+":
                temp.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (gene.contig, str(gene.start - upstream)
                                                             , str(gene.start + downstream), gene.gene_id
                                                             , ".", gene.strand))
            elif gene.strand == "-":
                temp.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (gene.contig, str(gene.end - downstream)
                                                             , str(gene.end + upstream), gene.gene_id
                                                             , ".", gene.strand))
        temp.close()
        inf = temp.name
    else:
        inf = infile
        concatenate = "zcat"

    # define statement
    # doesn't use strand information
    statement = '''%(concatenate)s %(inf)s | python %(scriptsdir)s/bed2fasta.py
                   --genome=%(genomedir)s/%(genome)s 
                   --log=%(outfile)s.log > %(outfile)s'''
    P.run()
    os.remove(inf)

###########    
# target
###########
@follows(buildIntervalsFasta)
def Fasta():
    pass

#########################################################################
#########################################################################
#########################################################################
# This section deals with comparisons in GC between the foreground and
# background sets
#########################################################################
#########################################################################
#########################################################################
@follows(mkdir("GC_content.dir"), buildIntervalsFasta)
@transform(INPUT_BACKGROUND + INPUT_FOREGROUND
           , regex(r"(\S+).tsv")
           , add_inputs(buildIntervalsFasta)
           , r"GC_content.dir/\1.gc.tsv")
def calculateGCContent(infiles, outfile):
    '''
    calculate the GC content across foreground and 
    background sets
    '''
    PipelineTransfacMatch.calculateCpGComposition(infiles[0], infiles[1], outfile)

#########################################################################
#########################################################################
#########################################################################
@transform(calculateGCContent, suffix(".tsv"), ".load")
def loadGCContent(infile, outfile):
    '''
    load the results the GC content for each background
    and foreground
    '''
    tablename = filenameToTablename(P.snip(os.path.basename(infile), ".tsv"))
    statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s
                   --log=%(outfile)s.log
                   --index=id
                   < %(infile)s > %(outfile)s'''
    P.run()

#########################################################################
#########################################################################
#########################################################################
if PARAMS["CpG_match_background"]:
    @collate(loadGCContent, regex("GC_content.dir/(.+)\.(?:background|foreground)\.gc\.load"),
             r"GC_content.dir/\1.cpg_matched.background.tsv")
    def matchBackgroundForCpGComposition(infiles, outfile):
        '''
        take the background set and subset it for intervals with
        a CpG distribution that is the same as the foreground set
        - this requires that the background set is sufficiently
        large
        '''
        track = re.match("GC_content.dir/(.+)\.(?:background|foreground)\.gc\.load", infiles[0]).groups()[0]
        input_background = "%s.background.tsv" % track
        input_foreground = "%s.foreground.tsv" % track
        
        PipelineTransfacMatch.matchBackgroundForCpGComposition( infiles
                                                                , input_background
                                                                , input_foreground
                                                                , PARAMS["database"]
                                                                , outfile )

    ###############################################
    ###############################################
    ###############################################
    @transform(matchBackgroundForCpGComposition
               , regex(r"(\S+).tsv")
               , add_inputs(buildIntervalsFasta)
               , r"\1.gc.tsv")
    def calculateCpGcomposition(infiles, outfile):
        '''
        calculate the GC content for the CpG matched data
        Should be the same as the CpG content of the foreground
        set
        '''
        PipelineTransfacMatch.calculateCpGComposition(infiles[0], infiles[1], outfile)

    ###############################################
    ###############################################
    ###############################################
    @transform(calculateCpGcomposition, suffix(".tsv"), ".load")
    def loadMatchedCpGComposition(infile, outfile):
        '''
        load the CpG compostion of matched background set
        '''
        tablename = filenameToTablename(P.snip(os.path.basename(outfile), ".load"))
        statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s 
                       --log=%(outfile)s.log
                       < %(infile)s > %(outfile)s'''
        P.run()

###############
# target
##############
if PARAMS["CpG_match_background"]:
    @follows(loadMatchedCpGComposition)
    def GC():
        pass
else:
    @follows(loadGCContent)
    def GC():
        pass

#########################################################################
#########################################################################
#########################################################################
@follows(mkdir("match.dir"))
@merge([buildIntervalsFasta
        , PARAMS["transfac_matrix"]
        , PARAMS["transfac_profile"]], "match.dir/match.result")
def runMatch(infiles, outfile):
    '''
    run transfac(R) match
    '''
    seq = infiles[0]
    mxlib = PARAMS["transfac_matrix"]
    mxprf = PARAMS["transfac_profile"]
    match_executable = "/ifs/data/biobase/transfac/match/bin/match_linux64"
    out = outfile
    
    statement = '''%(match_executable)s %(mxlib)s %(seq)s %(out)s %(mxprf)s'''
    P.run()
    
#########################################################################
#########################################################################
#########################################################################
@transform(runMatch, suffix(".result"), ".load")
def loadMatchResults(infile, outfile):
    '''
    load the results of the match analysis into 
    sqlite database
    '''
    temp = P.getTempFile()
    temp.write("seq_id\tmatrix_id\tposition\tstrand\tcore_score\tmatrix_score\tsequence\n")
    for details in PipelineTransfacMatch.match_iterator(infile):
        temp.write( "\t".join( map( str,[  details.seq_id
                                       , details.matrix_id 
                                       , details.position
                                       , details.strand
                                       , details.core_score 
                                       , details.matrix_score
                                       , details.sequence  ] ) ) + "\n" )
    inf = temp.name
    tablename = filenameToTablename(os.path.basename(infile))
    statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s
                   --log=%(outfile)s.log
                   --index=seq_id
                   < %(inf)s > %(outfile)s'''
    P.run()
    os.remove(inf)

#########################################################################
#########################################################################
#########################################################################
@transform(loadMatchResults, suffix(".load"), ".metrics")
def buildMatchMetrics(infile, outfile):
    '''
    match outputs transcription factors that are found in the supplied
    sequences. We are interested in the following metrics:
     * No. unique transcription factors found per sequence
     * Maximal number of TF motifs found per sequence
    '''
    tablename = filenameToTablename(os.path.basename(P.snip(infile, ".load"))) + "_result"
    PipelineTransfacMatch.frequencyMetrics(PARAMS["database"], tablename, outfile)

#########################################################################
#########################################################################
#########################################################################
@transform(buildMatchMetrics, suffix(".metrics"), ".load")
def loadMatchMetrics(infile, outfile):
    '''
    load match metrics
    '''
    tablename = filenameToTablename(os.path.basename(infile))
    statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s
                   --log=%(outfile)s.log
                   --index=seq_id
                   < %(infile)s > %(outfile)s'''
    P.run()

###########
# target
###########
@follows(loadMatchMetrics)
def Match():
    pass

#########################################################################
#########################################################################
#########################################################################
# The next section deals with testing for significant enrichment 
# between a foreground set of sequences and a background set
#########################################################################
#########################################################################
#########################################################################
if PARAMS["CpG_match_background"]:
    @follows(loadMatchResults, loadMatchedCpGComposition, mkdir("match_test.dir"))
    @collate([matchBackgroundForCpGComposition,calculateGCContent],
             regex(".+/(.+)\.(?:foreground.gc|cpg_matched.background)\.tsv"),
             r"match_test.dir/\1.cpg.matched.significance")
    def estimateEnrichmentOfTFBS(infiles, outfile):
        '''
        estimate the significance of trnascription factors that are associated with
        a foreground set of intervals vs a background set matched for CpG content
        '''
        # required files
        match_table = "match_result"

        #we don't know which order the foreground and backgorund will come in
        background = [infile for infile in infiles if re.search("background",infile)][0]
        foreground = ["%s.foreground.tsv" % re.match(".+/(.+)\.foreground\.gc\.tsv", infile).groups()[0] 
                      for infile in infiles if re.search("foreground", infile)][0]
        # run significance testing
        PipelineTransfacMatch.testSignificanceOfMatrices( background
                                                          , foreground
                                                          , PARAMS["database"]
                                                          , match_table
                                                          , outfile )

else:    
    @follows(loadMatchResults, mkdir("match_test.dir"))
    @collate(INPUT_BACKGROUND + INPUT_FOREGROUND, regex("(.+)\.(?:foreground|background)\.tsv")
           , r"match_test.dir/\1.significance")
    def estimateEnrichmentOfTFBS(infiles, outfile):
        '''
        estimate the significance of trnascription factors that are associated with
        a foreground set of intervals vs a background set
        '''
        # required files
        match_table = "match_result"

        #we don't know which order the foreground and backgorund will come in
        background = [infile for infile in infiles if re.search("background",infile)][0]
        foreground = [infile for infile in infiles if re.search("foreground",infile)][0]

        # run significance testing
        PipelineTransfacMatch.testSignificanceOfMatrices( background
                                                          , foreground
                                                          , PARAMS["database"]
                                                          , match_table
                                                          , outfile )

#########################################################################
#########################################################################
#########################################################################
@transform(estimateEnrichmentOfTFBS, suffix(".significance"), ".load")
def loadEnrichmentOfTFBS(infile, outfile):
    '''
    load the results of the enrichment
    '''
    tablename = filenameToTablename(os.path.basename(infile))
    statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s
                  --log=%(outfile)s.log 
                  --index=matrix_id
                  < %(infile)s > %(outfile)s'''
    P.run()
    
##############
# target
##############

@posttask(touch_file("complete.flag"))
@follows(loadEnrichmentOfTFBS)
def Significance():
    pass


@follows(Fasta
         , Match
         , Significance
         , GC)
def full():
    pass

#########################################################################
#########################################################################
#########################################################################
if __name__== "__main__":
    sys.exit( P.main(sys.argv) )
