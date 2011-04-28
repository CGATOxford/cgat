################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_vitaminD.py 2870 2010-03-03 10:20:29Z andreas $
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
pipeline_expression.py - analysis of expression data
====================================================

:Author: Andreas Heger
:Release: $Id: pipeline_vitaminD.py 2870 2010-03-03 10:20:29Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

This pipeline predicts gene that are differentially expressed.

This pipeline uses the bioconductor package.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

"""
import sys, tempfile, optparse, shutil, itertools, csv, math, random, re, glob, os, shutil, collections, ConfigParser

import Experiment as E
import Pipeline as P
from ruffus import *
import csv
import sqlite3
import IndexedFasta, IndexedGenome, FastaIterator, Genomics
import IOTools
import MAST, GTF, GFF, Bed, Stats
import cStringIO
import pysam
import numpy
import gzip
import Expression
import fileinput

if not os.path.exists("conf.py"):
    raise IOError( "could not find configuration file conf.py" )

execfile("conf.py")

TARGET_ANNOTATION= 'ensembl_regions.gff'
TARGET_GENESET= 'ensembl.gtf'
TARGET_PROMOTORS = 'promotors.gtf'
TARGET_TSS = 'tss.gtf'
TARGET_REPEATS = 'repeats.gff'
TARGET_TRANSCRIPTS = 'transcripts.gtf.gz'
TARGET_PROBESET = 'probeset.gtf'
TARGET_TRANSCRIPTS_TSS = 'transcripts_tss.gtf'
TARGET_TRANSCRIPTS_PROMOTORS = 'transcripts_promotors.gtf'
TARGET_ANNOTATOR_GENETERRITORIES='annotator_geneterritories.gff'
TARGET_MAPPABILITY='mappability.bed'

PARAMS = P.getParameters()

@files( ((None, "probeset2gene.table" ), ) )
def buildProbeset2Gene( infile, outfile ):
    '''build map relating a probeset to an ENSEMBL gene_id'''
    Expression.buildProbeset2Gene( infile, outfile )

@follows( buildProbeset2Gene )
def prepare(): pass    

@files( [ ( (x, "%s.map" % x), "%s_levels.import" % x[:-len("_series_matrix.txt.gz")] ) 
          for x in glob.glob("*_series_matrix.txt.gz") ] )
def importFromSeries( infiles, outfile ):
    '''import expression levels from a GEO series.'''
    tablename = P.toTable( outfile )
    
    tmpf = P.getTempFile()

    infile_data, infile_map = infiles

    map_header = IOTools.readMap( open(infile_map, "r") )
    if "ID_REF" not in map_header:
        map_header["ID_REF"] = "probeset"

    inf = gzip.open( infile_data, "r" )

    for line in inf:
        if line.startswith("!"): continue
        if not line.strip(): continue
        line = re.sub('"',"", line)
        if line.startswith( "ID_REF"):
            line = "\t".join( [map_header[x] for x in line[:-1].split("\t")] ) + "\n" 
                              
        tmpf.write( line )
        
    tmpf.close()
    tmpname = tmpf.name
    
    header = map_header["ID_REF"]
    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=%(header)s \
              --table=%(tablename)s \
    < %(tmpname)s > %(outfile)s
    '''

    P.run()
    os.unlink(tmpname)

@split( "import.map", "*.CEL" )
def importCEL( infile, outfiles ):
    '''import CEL files.'''
    
    map_cel = IOTools.readMap( open(infile, "r"), has_header = True )
    
    indir = PARAMS["datadir"]
    
    for old, new in map_cel.iteritems():
        oldname = os.path.join( indir, old + ".CEL" )
        newname = os.path.join( ".", new + ".CEL" )

        if not os.path.exists( oldname ):
            raise IOError( "input file %s does not exist" % oldname )
        if os.path.exists(newname): continue

        os.symlink( os.path.abspath( oldname), os.path.abspath( newname ) )
    
import rpy
from rpy import r as R

@merge( "*.CEL.gz", "expression.out" )
def estimateExpression( infiles, outfile ):
    '''estimate expression levels.'''
    
    R.library( "affy" )

    E.info( "reading data" )

    raw_data = R.ReadAffy( infiles )

    E.info( "RMA normalization" )

    eset = R.rma( raw_data )

    R.boxplot( raw_data )
    R.boxplot( eset )
    
    print R.as_list(R.assayData(eset))

def getTreatmentsAndControls( infile ):
    
    layout = ConfigParser.RawConfigParser()
    
    layout.read( infile )
    
    tablename = "levels"
    controls = layout.get("sam", "control").split(",")
    treatments = layout.get("sam", "treatment").split(",")
    
    dbhandle = sqlite3.connect( PARAMS["database"] )
    cc = dbhandle.cursor()
 
    treatment_columns = ",".join( ["%s" % x for x in treatments ] )
    control_columns = ",".join( ["%s" % x for x in controls ] )

    presence = "rma_presence"

    cc = dbhandle.cursor()
    statement = """SELECT l.probeset, 
              %(treatment_columns)s,
              %(control_columns)s
                        FROM %(tablename)s AS l, %(presence)s AS p
              WHERE p.probeset = l.probeset""" % locals()

    cc.execute( statement )

    r = zip(*cc.fetchall())

    probesets = r[0]
    treatments = r[1:len(treatments)+1]
    controls = r[len(treatments)+1:]

    return probesets, treatments, controls
        
@merge( "*.CEL", ("rma.Rdata", "levels.tsv" ) )
def computeExpressionLevels( infiles, outfiles ):
    '''normalize data using gcrma libary.
    
    output a file with the R object and
    another as human readable table.
    '''
    
    outfile_r, outfile_table = outfiles

    R.library( "simpleaffy" )
    R.library( "gcrma" )

    E.info( "reading data" )

    raw_data = R('''raw.data = ReadAffy()''')

    E.info( "normalization" )

    R('''gcrma.eset = call.exprs( raw.data, "%(normalization_method)s" )''' % PARAMS )

    E.info( "saving data" )
    R('''save( gcrma.eset, raw.data, file = "%s") ''' % outfile_r)

    data = R('''as.list(assayData(gcrma.eset))''')['exprs']
    probesets, headers = R('''dimnames( assayData(gcrma.eset)$exprs )''')
    headers = [ re.sub(".CEL", "", x) for x in headers ]

    outf = open( outfile_table, "w" )
    outf.write( "probeset\t%s\n" % "\t".join( headers) )

    for probeset, data in zip( probesets, data ):
        outf.write( "%s\t%s\n" % (probeset,
        "\t".join( map(str, data)) ) )
    outf.close()

@merge( computeExpressionLevels, "levels.import")
def importExpressionLevels( infiles, outfile ):
    '''import presence/absence data.'''
    
    infile_data, infile_table = infiles

    tablename = P.toTable( outfile )

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=probeset \
              --table=%(tablename)s \
    < %(infile_table)s > %(outfile)s
    '''

    P.run()

@transform( computeExpressionLevels, suffix(".Rdata"), ".presence")
def computePresenceCalls( infiles, outfile ):
    '''
    modified workflow from https://stat.ethz.ch/pipermail/bioconductor/2005-July/009758.html

    Calculate MAS5 detection pvalues and Present/Marginal/Absent                                                                                                                                                                            
    calls.

    Note that this analysis only applies to 3' arrays, but not for
    other arrays like HuGene. The simpleaffy package was designed
    to correct for 3' bias due to the 7' oligo dT primers used
    for building cDNAs. It requires PM and MM probes, which might
    not exist on other arrays.
    '''

    R.library( "simpleaffy" )
    R.library( "gcrma" )

    infile_r, infile_table = infiles

    # load the normalized data
    # (raw.data and gcrma.eset)
    R('''load( '%s' )''' % infile_r)

    # calculate detection calls
    call_eset = R('''call.eset = detection.p.val( raw.data, 
                                 alpha1=%(mas5_alpha1)f,
                                 alpha2=%(mas5_alpha2)f,
                  )''' % PARAMS )

    # output to files
    R('''write.table( call.eset$call, '%s.calls', 
                      sep='\t', 
                      quote=FALSE,
                      col.names = NA, row.names = TRUE)''' % outfile)
    R('''write.table( call.eset$pval, '%s.pval', 
                      sep='\t',
                      quote=FALSE,
                      col.names = NA, row.names = TRUE)''' % outfile)

    # build list of probesets to keep 
    # keep all probesets where at least 3 out of 4 replicates 
    # are indicated as present (P/M) in at least one treatment group
    R('''absent = rowSums( call.eset$call == 'A')''' ) 
    keep = R('''remove = names(absent[absent <= length(colnames(assayData( gcrma.eset )$exprs))-4] )''')
    outf = open( outfile, "w" )
    outf.write("probeset\n")
    outf.write("\n".join(map(str, keep)) + "\n" )
    outf.close()


@transform( computePresenceCalls, suffix(".presence"), "_presence.import")
def importPresence( infile, outfile ):
    '''import presence/absence data.'''
    
    tablename = P.toTable( outfile )

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=probeset \
              --table=%(tablename)s \
    < %(infile)s > %(outfile)s
    '''

    P.run()

@transform( "*.layout", suffix(".layout"), ".sam")
def computeDifferentialExpressionSAM( infile, outfile ):
    '''check for differential expression.
    
    Probesets that are called absent in either control or treatment
    are removed.
    '''

    probesets, treatments, controls = getTreatmentsAndControls( infile )
    
    E.info( "computing SAM for %i probesets, %i samples, %i control" % (
            len(probesets),
            len(treatments),
            len(controls)))

    target_path = os.path.join( os.path.abspath( PARAMS["exportdir"] ), "SAM" )
    if not os.path.exists( target_path): 
        try:
            os.makedirs( target_path )
        except OSError: 
            pass

    genes, summary, cutoffs = Expression.SAM()( probesets, treatments, controls,
                                                pattern = os.path.join(target_path, outfile + "%s"),
                                                fdr = float(PARAMS["sam_fdr"]),
                                                ngenes = float(PARAMS["sam_ngenes"]),
                                                ndelta = float(PARAMS["sam_ndelta"]),
                                                npermutations = PARAMS["sam_permutations"],
                                                method = PARAMS["sam_method"] )

    if summary == None:
        E.warn( "no cutoff when running sam for %s" % infile )

    if cutoffs:
        logs = open( outfile + ".cutoffs", "w" )
        logs.write("\t".join( cutoffs[0]._fields) + "\n" )
        for x in cutoffs:
            logs.write( "\t".join(map(str,x) ) + "\n" )
        logs.close()

    outs = open( outfile,"w")
    outs.write( "cluster_id\tmean1\tstd1\tmean2\tstd2\tpvalue\tqvalue\tdiff\tfold\tcalled\n" )
        
    for s in genes:
        outs.write ("%s\t%f\t%f\t%f\t%f\t%e\t%f\t%f\t%f\t%i\n" % \
                        (s.probeset, 
                         s.mean1,
                         s.stddev1,
                         s.mean2,
                         s.stddev2,
                         s.pvalue, 
                         s.qvalue, 
                         s.difference,
                         s.fold,
                         s.called,
                         ) )
    outs.close()




@follows( importExpressionLevels,
          importPresence )
def full():
    pass

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )
