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

:Author: Andreas Heger
:Release: $Id: pipeline_chains.py 2900 2010-04-13 14:38:00Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

build pairwise genomic alignments from a set of multiple alignments.

Starting from a multiple genomic alignment in maf, this pipeline builds
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
import sys, tempfile, optparse, shutil, itertools, csv, math, random, re, glob, os, shutil, collections

import Experiment as E
import Pipeline as P
from ruffus import *

PARAMS = P.getParameters()

def getGenomes( filename ):
    '''get the two genomes.
    '''
    return re.match("(\S+)To([^.]+)", filename ).groups()

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
             | python %(scriptsdir)s/blat2blat.py 
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
         | python %(scriptsdir)s/blat2blat.py 
              --method=remove-overlapping-query
              --log=%(outfile)s.log 
         | sort -k14,14 -k16,16n
         | python %(scriptsdir)s/blat2blat.py 
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
             | python %(scriptsdir)s/blat2blat.py 
                  --method=filter-fasta 
                  --method=sanitize
                  --filename-queries=%(genomefile)s
                  --filename-target=%(genome)s
                  --log=%(outfile)s.log 
             | gzip 
             >> %(outfile)s
             '''
    P.run()

@collate( extractPairwiseAlignment,
          regex(r"^(\S+).dir.*" ),
          r"\1To%s.psl.gz" % PARAMS["maf_master"])
def buildGenomeAlignment( infiles, outfile ):
    '''remove non-unique alignments in genomic infile.'''

    to_cluster = True

    infiles = " ".join(infiles)

    statement = '''zcat %(infiles)s 
         | sort -k10,10 -k12,12n
         | python %(scriptsdir)s/blat2blat.py 
              --method=remove-overlapping-query
              --log=%(outfile)s.log 
         | sort -k14,14 -k16,16n
         | python %(scriptsdir)s/blat2blat.py 
              --method=remove-overlapping-target
              --log=%(outfile)s.log 
         | gzip
         >> %(outfile)s
         '''
    P.run()

@transform( "*.chain.gz", suffix(".chain.gz"), ".psl.gz" )
def convertChainToPsl( infile, outfile ):
    '''convert a chain file to a psl file.'''

    to_cluster = False
    
    target, query = getGenomes( infile )
    
    E.debug( "query=%s, target=%s" % (query,target))

    statement = '''gunzip
    < %(infile)s 
    | %(cmd-farm)s --split-at-regex="^chain" --chunksize=1000 --max-lines=1000000
    " python %(scriptsdir)s/chain2psl.py --log=%(outfile)s.log
      | pslSwap stdin stdout "
    | gzip
    >  %(outfile)s
    '''
    
    P.run()

@files( [ (None, x + ".psl.gz", x) for x in P.asList(PARAMS["maps"]) ] )
def buildIndirectMaps( infile, outfile, track ):
    '''build a map between query and target, linking
    via intermediate targets.'''

    to_cluster = True
    
    path = PARAMS["%s_path" % track]

    E.info( "path=%s" % str(path))

    statement = []

    for stage,part in enumerate(path):
        filename = part + ".psl.gz"
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

@transform( (buildGenomeAlignment, convertChainToPsl, buildIndirectMaps), 
            suffix(".psl.gz"),
            ".stats")
def buildAlignmentStats( infile, outfile ):
    '''compute alignment coverage statistics'''
    
    to_cluster = True

    statement = '''
    gunzip < %(infile)s |
    python %(scriptsdir)s/psl2stats.py 
        --log=%(outfile)s.log 
    > %(outfile)s'''
    
    P.run()

@transform( buildIndirectMaps, 
            suffix(".psl.gz"),
            ".chain.gz" )
def convertPslToChain( infile, outfile ):
    '''convert a psl to a chain file.'''

    to_cluster = True

    statement = '''gunzip
    < %(infile)s
    | pslSwap stdin stdout
    | python %(scriptsdir)s/psl2chain.py
    | gzip
    > %(outfile)s'''
    P.run()

@follows( convertChainToPsl, buildGenomeAlignment )
def prepare(): 
    pass

@follows( buildIndirectMaps, buildAlignmentStats, convertPslToChain )
def build():
    pass

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )

