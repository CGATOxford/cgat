'''PipelineMedis.py

'''

import sys, tempfile, optparse, shutil, itertools, csv, math, random, re, glob, os, shutil, collections
import sqlite3
import cStringIO

import Experiment as E
import Pipeline as P

import csv
import IndexedFasta, IndexedGenome, FastaIterator, Genomics, Database
import IOTools
import GTF, GFF, Bed
# import Stats

import pysam
import numpy
import gzip
import fileinput

from rpy2.robjects import r as R
import rpy2.robjects as ro

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################
P.getParameters( 
    ["%s.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini" ] )

PARAMS = P.PARAMS

if os.path.exists("pipeline_conf.py"):
    E.info( "reading additional configuration from pipeline_conf.py" )
    execfile("pipeline_conf.py")

#########################################################################
#########################################################################
#########################################################################
def buildDMRStats( tables, method, outfile ):
    '''build dmr summary statistics.
    
    Creates some diagnostic plots in

    <exportdir>/<method> directory.

    Tables should be labeled <tileset>_<design>_<method>.

    '''

    dbhandle = sqlite3.connect( PARAMS["database"] )

    def togeneset( tablename ):
        return re.match("([^_]+)_", tablename ).groups()[0]

    keys_status = "OK", "NOTEST", "FAIL", "NOCALL"

    outf = IOTools.openFile( outfile, "w" )
    outf.write( "\t".join( ("tileset", "design", "track1", "track2", "tested",
                            "\t".join( [ "status_%s" % x for x in keys_status ] ),
                            "significant",
                            "twofold" ) ) + "\n" )

    all_tables = set(Database.getTables( dbhandle ))
    outdir = os.path.join( PARAMS["exportdir"], "diff_methylation" )

    for tablename in tables:

        prefix = P.snip( tablename, "_%s" % method )
        tileset, design = prefix.split("_")

        def toDict( vals, l = 2 ):
            return collections.defaultdict( int, [ (tuple( x[:l]), x[l]) for x in vals ] )

        E.info( "collecting data from %s" % tablename )
        tested = toDict( Database.executewait( dbhandle,
                                               """SELECT group1, group2, COUNT(*) FROM %(tablename)s 
                                GROUP BY group1,group2""" % locals() ).fetchall() )
        status = toDict( Database.executewait( dbhandle,
                                               """SELECT group1, group2, status, COUNT(*) FROM %(tablename)s 
                                GROUP BY group1,group2,status""" % locals() ).fetchall(), 3 )
        signif = toDict( Database.executewait( dbhandle,
                                               """SELECT group1, group2, COUNT(*) FROM %(tablename)s 
                                WHERE significant
                                GROUP BY group1,group2""" % locals() ).fetchall() )
        fold2 = toDict( Database.executewait( dbhandle,
                """SELECT group1, group2, COUNT(*) FROM %(tablename)s 
                                WHERE (lfold >= 1 or lfold <= -1) AND significant
                                GROUP BY group1,group2,significant""" % locals() ).fetchall() )
        
        groups = tested.keys()

        for group1, group2 in groups:
            outf.write( "\t".join(map(str, (
                            tileset,
                            design,
                            group1,
                            group2,
                            tested[(group1,group2)],
                            "\t".join( [ str(status[(group1,group2,x)]) for x in keys_status]),
                            signif[(group1,group2)],
                            fold2[(group1,group2)] ) ) ) + "\n" )

        ###########################################
        ###########################################
        ###########################################
        # plot length versus P-Value
        data = Database.executewait( dbhandle, 
                                     '''SELECT end - start, pvalue 
                             FROM %(tablename)s
                             WHERE significant'''% locals() ).fetchall()

        # require at least 10 datapoints - otherwise smooth scatter fails
        if len(data) > 10:
            data = zip(*data)

            pngfile = "%(outdir)s/%(tileset)s_%(design)s_%(method)s_pvalue_vs_length.png" % locals()
            R.png( pngfile )
            R.smoothScatter( R.log10( ro.FloatVector(data[0]) ),
                             R.log10( ro.FloatVector(data[1]) ),
                             xlab = 'log10( length )',
                             ylab = 'log10( pvalue )',
                             log="x", pch=20, cex=.1 )

            R['dev.off']()

    outf.close()

