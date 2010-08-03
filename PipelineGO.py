################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: PipelineGO.py 2877 2010-03-27 17:42:26Z andreas $
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
:Release: $Id: PipelineGO.py 2877 2010-03-27 17:42:26Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

Pipeline components - GO analysis

Tasks related to gene set GO analysis.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----


"""
import sys, tempfile, optparse, shutil, itertools, csv, math, random, re, glob, os, shutil, collections
import sqlite3

import Experiment as E
import Pipeline as P
import Stats

try:
    PARAMS = P.getParameters()
except IOError:
    pass

############################################################
############################################################
############################################################
## get GO assignments
def createGO( infile, outfile ):
    '''get GO assignments from ENSEMBL'''
    
    statement = '''
        python %(scriptsdir)s/GO.py 
                     --filename-dump=%(outfile)s 
                     --host=%(go_host)s 
                     --user=anonymous 
                     --database=%(go_database)s 
                     --port=5306 > %(outfile)s.log
        '''

    P.run()

############################################################
############################################################
############################################################
## get GO Slim assignments
############################################################
def createGOSlim( infile, outfile ):
    '''get GO assignments from ENSEMBL'''
    
    statement = '''wget %(go_url_goslim)s'''
    P.run()

    statement = '''wget %(go_url_ontology)s'''
    P.run()

    
    statement = '''
        /net/cpp-group/server/lib/perl5/5.10.0/bin/map2slim -outmap go2goslim.map goslim_goa.obo gene_ontology.obo
    '''
    P.run()

    statement = '''
        python %(scriptsdir)s/GO.py \
                --go2goslim \
                --filename-ontology=gene_ontology.obo \
                --slims=go2goslim.map \
                --log=%(outfile)s.log \
        < %(infile)s > %(outfile)s
        '''
    P.run()

############################################################
def runGOFromFiles( outfile,
                    outdir,
                    fg_file,
                    bg_file,
                    go_file,
                    ontology_file = None,
                    samples = 1000 ):
    '''check for GO enrichment within a gene list.

    The gene list is given in ``fg_file``. It is compared
    against ``bg_file`` using the GO assignments from
    ``go_file``. Results are saved in ``outfile`` and
    ``outdir``.
    '''

    to_cluster = True
    
    if ontology_file == None: ontology_file = PARAMS["go_ontology"]

    statement = '''
    python %(scriptsdir)s/GO.py 
        --filename-input=%(go_file)s 
        --genes=%(fg_file)s 
        --background=%(bg_file)s 
        --sample=%(samples)s
        --fdr 
        --filename-ontology=%(ontology_file)s 
        --output-filename-pattern='%(outdir)s/%%(go)s.%%(section)s' 
    > %(outfile)s'''

    P.run()    

    dbhandle = sqlite3.connect( PARAMS["database"] )

############################################################
def runGOFromDatabase( outfile, outdir, 
                       statement_fg, 
                       statement_bg, 
                       go_file,
                       ontology_file = None,
                       samples = 1000 ):
    '''Take gene lists from the SQL database using
    ``statement_foreground`` and ``statement_background``
    '''

    dbhandle = sqlite3.connect( PARAMS["database"] )
    
    cc = dbhandle.cursor()
    fg = set( [x[0] for x in cc.execute( statement_fg).fetchall() ] )
    bg = set( [x[0] for x in cc.execute( statement_bg).fetchall() ] )

    if len(fg) == 0:
        P.touch( outfile )
        return

    fg_file = os.path.join( outdir, "foreground" )
    bg_file = os.path.join( outdir, "background" )
    outf = open( fg_file, "w")
    outf.write("\n".join( map(str, fg ) ) + "\n" )
    outf.close()
    outf = open( bg_file, "w")
    outf.write("\n".join( map(str, bg ) ) + "\n" )
    outf.close()
    
    runGOFromFiles( outfile, outdir, 
                    fg_file, bg_file, 
                    go_file,
                    ontology_file = ontology_file,
                    samples = samples )

############################################################
############################################################
############################################################
##
############################################################
def loadGO( infile, outfile, tablename ):
    '''import GO results into individual tables.'''

    indir = infile + ".dir"

    if not os.path.exists( indir ):
        P.touch( outfile )
        return

    statement = '''
    python %(toolsdir)s/cat_tables.py %(indir)s/*.overall |\
    csv2db.py %(csv2db_options)s \
              --allow-empty \
              --index=category \
              --index=goid \
              --table=%(tablename)s \
    > %(outfile)s
    '''
    P.run()


############################################################
############################################################
############################################################
##
############################################################
def loadGOs( infiles, outfile, tablename ):
    '''import GO results into a single table.

    This method also computes a global QValue over all
    tracks, genesets and annotation sets.
    '''

    header = False

    tempf1 = P.getTempFile()

    pvalues = []

    for infile in infiles:
        indir = infile + ".dir"

        if not os.path.exists( indir ):
            continue

        track, geneset, annotationset = re.search("^(\S+)_vs_(\S+)\.(\S+)", infile ).groups()

        for filename in glob.glob( os.path.join(indir, "*.overall") ):
            for line in open(filename, "r" ):
                if line.startswith("#"): continue
                data = line[:-1].split("\t")
                if line.startswith("code"):
                    if header: continue
                    tempf1.write( "track\tgeneset\tannotationset\t%s" % line )
                    header = True
                    assert data[10] == "pover" and data[11] == "punder", "format error, expected pover-punder, got %s-%s" % (data[10], data[11])
                    continue
                tempf1.write( "%s\t%s\t%s\t%s" % (track, geneset, annotationset, line) )
                pvalues.append( min( float(data[10]), float(data[11]) ) )

    tempf1.close()

    E.info( "analysing %i pvalues" % len(pvalues ))
    fdr = Stats.doFDR( pvalues )
    E.info( "got %i qvalues" % len(fdr.mQValues ))
    qvalues = ["global_qvalue" ] + fdr.mQValues

    tempf2 = P.getTempFile()

    for line, qvalue in zip( open(tempf1.name,"r"), qvalues ):
        tempf2.write( "%s\t%s\n" % (line[:-1], str(qvalue)) )

    tempf2.close()
    tempfilename = tempf2.name
    print tempf1.name
    print tempf2.name

    statement = '''
    csv2db.py %(csv2db_options)s 
              --allow-empty 
              --index=category 
              --index=track,geneset,annotationset
              --index=geneset
              --index=annotationset
              --index=goid 
              --table=%(tablename)s 
    < %(tempfilename)s
    > %(outfile)s
    '''
    P.run()

    #os.unlink( tempf1.name )
    #os.unlink( tempf2.name )
    
