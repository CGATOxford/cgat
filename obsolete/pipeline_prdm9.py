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
pipeline_prdm9.py
=================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Analysis of prdm9 in 17 mouse strains.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

"""
from ruffus import *
import sys
import glob
import gzip
import os
import itertools
import CGAT.CSV as CSV
import re
import math
import types
import optparse
import shutil
import sqlite3
import CGAT.GFF as GFF
import CGAT.GTF as GTF
import CGAT.Experiment as E
import CGAT.Pipeline as P
import CGAT.IOTools as IOTools
import CGAT.Genomics as Genomics
import CGAT.Database as Database
import CGAT.FastaIterator as FastaIterator
import PipelineGeneset as PGeneset
import PipelineGO as PGO
import scipy.stats
import CGAT.Stats as Stats
import alignlib
import CGAT.Mali as Mali

PARAMS = P.getParameters()

@files( (("../data/znf.data", "profile.fasta" ) , ))
def createProfile( infile, outfile ):
    '''convert mali to profile
    '''
    outfile = open(outfile, "w")
    for line in open(infile):
        if line.startswith("#"): continue
        data = re.split("\s+", line[:-1],1)
        print data
        pid, sequence = data
        outfile.write(">%s\n%s\n" % (pid, sequence ))
    outfile.close()

def getParts( src ):
    '''split a wrap-around alignment'''

    result = None
    r = []
    last_s = src.getColTo()
    for p in range( src.getRowFrom(), 
                    src.getRowTo() ):
        s = src.mapRowToCol(p)
        if s < 0: continue
        if last_s >= s:
            if result:
                r.append( result )
            result = alignlib.makeAlignmentVector()
        last_s = s
        result.addPair( s, p, 0 )

    if result:
        r.append( result )
    return r

def alignIndels( all_alleles, colcounts, extend_by = 0 ):
    '''align all indel-regions.'''

    aa = alignlib.makeAlignatorDPFull( alignlib.ALIGNMENT_LOCAL, 0, 0 )     
    alignator = alignlib.makeMultipleAlignatorSimple( aa)

    ids = all_alleles.keys()

    for x,c in enumerate(colcounts):
        if c <= 1: continue
        sequences = alignlib.StringVector()
        for sid in ids:
            for allele in all_alleles[sid]:
                sequences.append( allele[x] )

        mali = alignlib.makeMultAlignment()
        alignator.align( mali, sequences )
        realigned = []
        for line in str(alignlib.MultAlignmentFormatPlain( mali, sequences )).split("\n")[:-1]:
            data = line[:-1].split("\t")
            realigned.append( data[1] )
        assert len(realigned) == len(sequences)

        l = max( [len(r) for r in realigned] )
        i = 0
        for sid in ids:
            for allele in all_alleles[sid]:
                if realigned[i]: allele[x] = realigned[i] 
                else: allele[x] = "-" * l 
                i += 1
                
        colcounts[x] = l

def _alignToProfile( infile, outfile, 
                     min_score = 0 ):
    '''align sequences in *infile* against mali

    Only alignments with a score higher than *min_score* are accepted.

    Output multiple alignment in fasta format to *outfile* and a table
    in :file:`outfile.log`.
    '''

    mali = Mali.Mali()
    mali.readFromFile( open("../data/mouse.fasta") )
    src_mali = Mali.convertMali2Alignlib( mali )
    
    E.debug( "read mali: %i sequences x %i columns" % (mali.getNumSequences(), mali.getNumColumns() ))

    # add pseudocounts
    profile_mali = mali.getClone()
    n = profile_mali.getNumColumns() 
    for x in "ACGT": 
        for y in range(0,2):
            profile_mali.addSequence( "%s%i" % (x,y), 0, n, x * n )


    profile_mali = Mali.convertMali2Alignlib( profile_mali )
    alignlib.setDefaultEncoder( alignlib.getEncoder( alignlib.DNA4 ) )
    alignlib.setDefaultLogOddor( alignlib.makeLogOddorUniform() )

    # bg = alignlib.FrequencyVector()
    # bg.extend( ( 0.3, 0.1, 0.2, 0.2, 0.2) )
    # alignlib.setDefaultRegularizor( alignlib.makeRegularizorTatusov(
    #         alignlib.makeSubstitutionMatrixDNA4(),
    #         bg,
    #         "ACGTN",
    #         10.0, 1.0) )

    profile = alignlib.makeProfile( profile_mali )
    
    alignment_mode = alignlib.ALIGNMENT_WRAP

    alignator = alignlib.makeAlignatorDPFull( alignment_mode,
                                              -5.0,
                                              -0.5 )
    
    map_seq2profile = alignlib.makeAlignmentVector()
    map_rseq2profile = alignlib.makeAlignmentVector()
    profile.prepare()

    # print profile

    build_mali = alignlib.makeMultAlignment()
    m = alignlib.makeAlignmentVector()
    m.addDiagonal( 0, n, 0 )
    build_mali.add( src_mali, m )

    outf = open( outfile, "w" )
    outf_log = open( outfile + ".info", "w" )
    outf_log.write( "read_id\tlength\tstart\tend\tparts\tcovered\tpcovered\tscore\tmali_start\tmali_end\tmali_covered\tmali_pcovered\n" )

    sequences, aa = alignlib.StringVector(), alignlib.AlignandumVector()
    ids = []

    for pid in mali.getIdentifiers():
        sequences.append( re.sub( "-", "", mali[pid] ) )
        ids.append( pid )

    # print str(alignlib.MultAlignmentFormatPlain( build_mali, sequences ))

    c = E.Counter()

    for s in FastaIterator.FastaIterator( open(infile)):

        E.debug("adding %s" % s.title )
        c.input += 1
        rsequence = Genomics.complement(s.sequence)
        seq = alignlib.makeSequence( s.sequence )
        rseq = alignlib.makeSequence( rsequence )

        alignator.align( map_seq2profile, seq, profile )
        alignator.align( map_rseq2profile, rseq, profile )

        if map_seq2profile.getScore() > map_rseq2profile.getScore():
            m, seq, sequence = map_seq2profile, seq, s.sequence
        else:
            m, seq, sequence = map_rseq2profile, rseq, rsequence

        if m.getLength() == 0:
            c.skipped += 1
            continue

        if m.getScore() < min_score: 
            c.skipped += 1
            continue

        r = getParts( m )

        covered = 0
        for mm in r:
            build_mali.add( mm )
            sequences.append( sequence )
            ids.append( s.title )
            covered += mm.getLength() - mm.getNumGaps()

        mali_covered = m.getColTo() - m.getColFrom()

        outf_log.write( "\t".join( map(str, (
                        s.title,
                        len(s.sequence),
                        m.getRowFrom(),
                        m.getRowTo(),
                        len(r),
                        covered,
                        "%5.2f" % (100.0 * covered / len(s.sequence) ),
                        m.getScore(),
                        m.getColFrom(),
                        m.getColTo(),
                        mali_covered,
                        "%5.2f" % ((100.0 * mali_covered) / mali.getNumColumns())
                        ) ) ) + "\n" )

        c.output += 1

    #build_mali.expand( aa )
    result = str(alignlib.MultAlignmentFormatPlain( build_mali, 
                                                    sequences, 
                                                    alignlib.UnalignedStacked ))

    for pid, data in zip(ids, result.split("\n") ):
        start, sequence, end = data.split("\t")
        outf.write(">%s/%i-%i\n%s\n" % (pid, int(start)+1, int(end), sequence) )


    outf.close()
    outf_log.close()

    E.info( "%s\n" % str(c) )

@follows(createProfile)
@files( [ (x, "%s_%03i_na.mali" % (x[:-3],f), f) 
          for x,f in itertools.product( glob.glob("*.fa"), (0, 80 ) )  ] )
def alignToProfile( infile, outfile, threshold ):
    _alignToProfile( infile, outfile, min_score = threshold )

@transform( alignToProfile
            , suffix( ".mali")
            , ".import" )
def importMaliStats( infile, outfile ):
    '''import stats.'''
    
    table = P.toTable( outfile ) + "_info" 
    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --index=read_id 
              --table=%(table)s 
    < %(infile)s.info
    > %(outfile)s
    '''
    P.run()

@transform( alignToProfile
            , suffix( "_na.mali")
            , "_aa.mali" )
def buildMalisCodon( infile, outfile ):
    '''build codon alignments

    The frame is given by the master multiple alignment.

    Sequences with indels that are not multiples of three
    are removed.

    Sequences with stop codons are removed.
    TODO
    '''
    
    statement = '''
    cat %(infile)s
    | python %(scriptsdir)s/fasta2fasta.py --method=translate 
    | perl -p -e "s/[a-z]/-/g unless (/^>/)"
    | python %(scriptsdir)s/mali2mali.py --method=remove-all-gaps --log=%(outfile)s.log --allow-duplicates

    > %(outfile)s
    '''
    
    P.run()


@transform( alignToProfile
            , suffix( "_na.mali")
            , "_aa.mali" )
def buildMalisAA( infile, outfile ):
    '''translate multiple alignments.

    The frame is given by the master multiple alignment.

    Indels in other sequences are ignored.
    '''
    
    statement = '''
    cat %(infile)s
    | perl -p -e "s/[a-z]/-/g unless (/^>/)"
    | python %(scriptsdir)s/mali2mali.py --method=remove-all-gaps --log=%(outfile)s.log --allow-duplicates
    | python %(scriptsdir)s/fasta2fasta.py --method=translate 
    > %(outfile)s
    '''
    
    P.run()

@transform( alignToProfile
            , suffix( "_na.mali")
            , "_na.columns" )
def computeColumnStatsNA( infile, outfile ):
    '''compute stats per column.

    Only columns are counted that are non-insert with respect to
    the master alignment. The master alignment columns are
    later removed in order to count only sequences within 
    a strain.
    '''
    
    statement = '''
    cat %(infile)s
    | perl -p -e "s/[a-z]/-/g unless (/^>/)"
    | python %(scriptsdir)s/mali2mali.py --method=remove-all-gaps --log=%(outfile)s.log --allow-duplicates
    | python %(scriptsdir)s/fasta2fasta.py --exclude="(^mouse)" --log=%(outfile)s.log
    | python %(scriptsdir)s/mali2table.py --section=all --log=%(outfile)s.log --allow-duplicates
    > %(outfile)s
    '''
    
    P.run()

@transform( buildMalisAA
            , suffix( "_aa.mali")
            , "_aa.columns" )
def computeColumnStatsAA( infile, outfile ):
    '''compute stats per column.

    Only columns are counted that are non-insert with respect to
    the master alignment. The master alignment columns are
    later removed.
    '''
    
    statement = '''
    cat %(infile)s
    | perl -p -e "s/[X?]/-/g unless (/^>/)"
    | python %(scriptsdir)s/fasta2fasta.py --exclude="(^mouse)" --log=%(outfile)s.log
    | python %(scriptsdir)s/mali2table.py --section=all --log=%(outfile)s.log --allow-duplicates --alphabet=aa
    > %(outfile)s
    '''
    
    P.run()

@transform( (computeColumnStatsNA, computeColumnStatsAA)
            , suffix( ".columns")
            , "_columns.import" )
def importColumnStats( infile, outfile ):
    '''import stats.'''
    
    table = P.toTable( outfile ) 
    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --index=column 
              --table=%(table)s 
    < %(infile)s
    > %(outfile)s
    '''
    P.run()


@follows( alignToProfile
          , importMaliStats
          , computeColumnStatsNA
          , buildMalisAA
          , computeColumnStatsAA
          , importColumnStats
          )
def full(): pass

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )


