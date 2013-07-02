################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_chains.py 2900 2010-04-13 14:38:00Z andreas $
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
===============
Chains pipeline
===============

:Author: Andreas Heger
:Release: $Id: pipeline_chains.py 2900 2010-04-13 14:38:00Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

build pairwise genomic alignments from a set of multiple alignments or combine several
pairwise alignments into new pairwise alignments.

Starting from a multiple genomic alignment in :term:`maf` format, this pipeline builds
a set of pairwise mappings between various query genomes and a single
reference or target genome. 

The output is a set of `chain <http://www.breyer.com/ucsc/htdocs/goldenPath/help/chain.html>`_
formatted files that can be used to map feature sets from one assembly to
another using :file:`liftOver`.

The script works by creating ``psl`` formatted files first. It will then filter and merge
these. ``chain`` formatted files are created in the last stage.

.. note::

   There is a limit on the divergence between assemblies beyond which features can not be mapped accurately.
   At high divergence, features should better mapped using tools like :file:`gmap`.

Liftover chain files are named <target>To<query>.over.chain.gz. Liftover files
downloaded from UCSC are unique for query, but can have overlapping targets.

.. note::
   This script removes all 1-to-many, many-to-1 and many-to-many mappings.
   It does so in a greedy way by removing all segments that overlap with any other segment.
   A single base overlap is sufficient for overlap and all overlapping segments are removed.

Nomenclature
++++++++++++

The nomenclature in the UCSC is :file:`TargetToQuery.chain` for mapping ``target``
to ``query`` (according to the UCSC documentation, ``target`` is the first entry
in ``chain`` files).

I have been using the nomenclature `QueryToTarget.psl`. 

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

To use the pipeline, run the following tasks

prepare
   extract pairwise alignments from a maf alignment
   and convert chain files to psl files

build
   build the required pairwise alignments


Code
----

"""
import sys
import tempfile
import optparse
import shutil
import itertools
import csv
import math
import random
import re
import glob
import os
import shutil
import collections

import CGAT.Experiment as E
import CGAT.Pipeline as P

import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
from ruffus import *

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import CGAT.Pipeline as P
P.getParameters( 
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini" ],
    defaults = { "maps" : "" } )
PARAMS = P.PARAMS

###################################################################
###################################################################
###################################################################
##
###################################################################
if os.path.exists("pipeline_conf.py"): 
    L.info( "reading additional configuration from pipeline_conf.py" )
    execfile("pipeline_conf.py")

PARAMS = P.getParameters()

###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
import CGATPipelines.PipelineTracks as PipelineTracks

###################################################################
def extractGenomes( filename ):
    '''extract the two genomes from a chain filename.
    It also decapitalizes them.
    '''
    return [ x[0].lower() + x[1:] for x in re.match("(\S+)To([^.]+)", filename ).groups() ]

def writeContigSizes( genome, outfile ):
    '''write contig sizes to outfile for UCSC tools.
    '''

    outf = IOTools.openFile( outfile, "w")
    fasta= IndexedFasta.IndexedFasta( os.path.join( PARAMS["genome_dir"], genome ) )
    for contig, size in fasta.getContigSizes( with_synonyms = False ).iteritems():
        outf.write( "%s\t%i\n" % (contig, size ) )
    outf.close()

@follows( mkdir( "export" ) )
def prepare():
    pass

####################################################################
####################################################################
####################################################################
## 
####################################################################
@follows(prepare)
@transform( "*.chain.gz", suffix(".chain.gz"), ".psl.gz" )
def convertChainToPsl( infile, outfile ):
    '''convert a chain file to a psl file.
    '''

    to_cluster = False
    
    target, query = extractGenomes( infile )
    
    E.debug( "query=%s, target=%s" % (query,target))

    statement = '''gunzip
    < %(infile)s 
    | %(cmd-farm)s --split-at-regex="^chain" --chunksize=1000 --max-lines=1000000 --log=%(outfile)s.log
    " python %(scriptsdir)s/chain2psl.py --log=%(outfile)s.log
      | pslSwap stdin stdout "
    | gzip
    >  %(outfile)s
    '''

    P.run()

##################################################################################
##################################################################################
##################################################################################
## extracting alignments from maf files
##################################################################################
if "maf_dir" in PARAMS and "maf_tracks" in PARAMS:
    @files( [ ( ("%s/*.maf.gz" % PARAMS["maf_dir"]) 
                , "%sTo%s.raw.psl.gz" % (PARAMS["%s_label" % track], PARAMS["maf_master"])
                , track) for track in P.asList(PARAMS["maf_tracks"]) ] )
    def extractPairwiseAlignmentSingleFile(infiles, outfile, track):
        '''build pairwise genomic aligment from maf files.'''

        try: os.remove( outfile )
        except OSError: pass

        genomefile = PARAMS["%s_genome" % track ]

        to_cluster = True

        for infile in infiles:

            E.info( "adding %s" % infile )

            statement = '''gunzip < %(infile)s 
                 | python %(scriptsdir)s/maf2psl.py 
                      --query=%(track)s
                      --target=%(maf_master)s
                      --log=%(outfile)s.log 
                 | python %(scriptsdir)s/psl2psl.py 
                      --method=filter-fasta 
                      --method=sanitize
                      --filename-queries=%(genomefile)s
                      --filename-target=%(genome)s
                      --log=%(outfile)s.log 
                 | gzip 
                 >> %(outfile)s
                 '''
            P.run()

    @transform( extractPairwiseAlignmentSingleFile
                , suffix(".raw.psl.gz")
                , ".psl.gz" )
    def buildGenomeAlignmentFromSingleFile( infile, outfile ):
        '''remove non-unique alignments in genomic infile.'''

        to_cluster = True

        statement = '''gunzip < %(infile)s 
             | sort -k10,10 -k12,12n
             | python %(scriptsdir)s/psl2psl.py 
                  --method=remove-overlapping-query
                  --log=%(outfile)s.log 
             | sort -k14,14 -k16,16n
             | python %(scriptsdir)s/psl2psl.py 
                  --method=remove-overlapping-target
                  --log=%(outfile)s.log 
             | gzip
             >> %(outfile)s
             '''
        P.run()

    @follows( mkdir( ("%s.dir" % track for track in P.asList(PARAMS["maf_tracks"]) ) ) )
    @files( [ (infile, "%s.dir/%s" % (track, os.path.basename(infile)), track) \
                  for infile, track in itertools.product( 
                glob.glob("%s/*.maf.gz" % PARAMS["maf_dir"]),  P.asList(PARAMS["maf_tracks"])) ] )
    def extractPairwiseAlignment(infile, outfile, track):
        '''build pairwise genomic aligment from maf files.'''

        genomefile = PARAMS["%s_genome" % track ]
        query=PARAMS["%s_label" % track]

        to_cluster = True

        statement = '''gunzip < %(infile)s 
                 | python %(scriptsdir)s/maf2psl.py 
                      --query=%(query)s
                      --target=%(maf_master)s
                      --log=%(outfile)s.log 
                 | python %(scriptsdir)s/psl2psl.py 
                      --method=filter-fasta 
                      --method=sanitize
                      --filename-queries=%(genomefile)s
                      --filename-target=%(genome)s
                      --log=%(outfile)s.log 
                 | gzip 
                 >> %(outfile)s
                 '''
        P.run()

    @follows(prepare)
    @collate( extractPairwiseAlignment,
              regex(r"^(\S+).dir.*" ),
              r"\1To%s.psl.gz" % PARAMS["maf_master"])
    def buildGenomeAlignment( infiles, outfile ):
        '''remove non-unique alignments in genomic infile.'''

        to_cluster = True

        infiles = " ".join(infiles)

        statement = '''zcat %(infiles)s 
             | sort -k10,10 -k12,12n
             | python %(scriptsdir)s/psl2psl.py 
                  --method=remove-overlapping-query
                  --log=%(outfile)s.log 
             | sort -k14,14 -k16,16n
             | python %(scriptsdir)s/psl2psl.py 
                  --method=remove-overlapping-target
                  --log=%(outfile)s.log 
             | gzip
             >> %(outfile)s
             '''
        P.run()

else:
    @follows(prepare)
    def buildGenomeAlignment():
        pass

@follows( buildGenomeAlignment, convertChainToPsl )
@files( [ (None, x + ".over.psl.gz", x) for x in P.asList(PARAMS["maps"]) ] )
def buildIndirectMaps( infile, outfile, track ):
    '''build a map between query and target, linking
    via intermediate targets.'''

    to_cluster = True
    
    path = P.asList(PARAMS["%s_path" % track])

    E.info( "path=%s" % str(path))

    statement = []

    for stage, part in enumerate(path):
        filename = part + ".over.psl.gz"
        if not os.path.exists( filename ):
            raise ValueError( "required file %s for %s (stage %i) not exist." % (filename, outfile, stage ))

        if stage == 0:
            statement.append( '''gunzip < %(filename)s''' % locals() )
        else:
            statement.append( '''
               pslMap stdin <(gunzip < %(filename)s) stdout
            ''' % locals() )

    statement.append( "gzip" )
    
    statement = " | ".join( statement ) + " > %(outfile)s " % locals()

    P.run() 

##################################################################################
##################################################################################
##################################################################################
## 
##################################################################################
@transform( buildIndirectMaps, 
            regex( r"(.*).psl.gz") ,
            r"export/\1.chain.gz" )
def convertPslToChain( infile, outfile ):
    '''convert a psl to a chain file.

    see http://genomewiki.ucsc.edu/index.php/Minimal_Steps_For_LiftOver
    '''

    to_cluster = True

    target, query = extractGenomes( infile )
    
    tmpfilename1 = P.getTempFilename( ".")
    tmpfilename2 = P.getTempFilename( ".")

    writeContigSizes( target, tmpfilename1 )
    writeContigSizes( query, tmpfilename2 )

    statement = '''gunzip
    < %(infile)s
    | pslSwap stdin stdout
    | python %(scriptsdir)s/psl2chain.py --log=%(outfile)s.log
    | chainSort stdin stdout
    | gzip
    > %(outfile)s.sorted.chain.gz;
    checkpoint; 
    gunzip < %(outfile)s.sorted.chain.gz 
    | chainNet stdin %(tmpfilename1)s %(tmpfilename2)s stdout /dev/null
    | netChainSubset stdin <( zcat %(outfile)s.sorted.chain ) stdout
    | gzip
    > %(outfile)s'''
    P.run()

    os.unlink( tmpfilename1 )
    os.unlink( tmpfilename2 )

##################################################################################
##################################################################################
##################################################################################
## 
##################################################################################
@transform( (buildGenomeAlignment, convertChainToPsl, buildIndirectMaps), 
            suffix(".psl.gz"),
            ".psl.stats")
def buildPslStats( infile, outfile ):
    '''compute alignment coverage statistics in chain files
    '''
    
    to_cluster = True

    statement = '''
    gunzip < %(infile)s 
    | python %(scriptsdir)s/psl2stats.py --log=%(outfile)s.log
    > %(outfile)s'''
    
    P.run()

##################################################################################
##################################################################################
##################################################################################
## 
##################################################################################
@transform( convertPslToChain,
            suffix(".chain.gz"),
            ".chain.stats")
def buildChainStats( infile, outfile ):
    '''compute alignment coverage statistics in chain files
    '''
    
    to_cluster = True

    statement = '''
    gunzip < %(infile)s 
    | python %(scriptsdir)s/chain2psl.py --log=%(outfile)s.log
    | python %(scriptsdir)s/psl2stats.py --log=%(outfile)s.log
    > %(outfile)s'''
    
    P.run()

##################################################################################
##################################################################################
##################################################################################
## 
##################################################################################
@follows( buildIndirectMaps, 
          convertPslToChain, 
          buildPslStats,
          buildChainStats )
def full():
    pass

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
    sys.exit( P.main(sys.argv) )

