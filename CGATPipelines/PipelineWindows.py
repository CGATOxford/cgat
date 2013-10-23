import os
import re
import sys

from ruffus import *

import CGAT.Experiment as E
import CGAT.Pipeline as P
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import CGAT.Stats as Stats

try:
    PARAMS = P.getParameters()
except IOError:
    pass

#########################################################################
#########################################################################
#########################################################################
def countReadsWithinWindows( bedfile, windowfile, outfile, counting_method = "midpoint" ):
    '''count reads in *bedfile* within *windowfile*.

    '''
    to_cluster = True

    job_options = "-l mem_free=4G"

    if counting_method == "midpoint":
        f = '''| awk '{a = $2+($3-$2)/2; printf("%s\\t%i\\t%i\\n", $1, a, a+1)}' '''
    elif countig_method == "nucleotide":
        f = ""
    else: 
        raise ValueError("unknown counting method: %s" % counting_method )

    statement = '''
    zcat %(bedfile)s
    %(f)s
    | coverageBed -a stdin -b %(windowfile)s -split
    | sort -k1,1 -k2,2n 
    | gzip
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
def aggregateWindowsReadCounts( infiles, outfile ):
    '''aggregate tag counts for each window.

    coverageBed outputs the following columns:
    1) Contig
    2) Start
    3) Stop
    4) Name
    5) The number of features in A that overlapped (by at least one base pair) the B interval.
    6) The number of bases in B that had non-zero coverage from features in A.
    7) The length of the entry in B.
    8) The fraction of bases in B that had non-zero coverage from features in A.

    For bed: use column 5
    For bed6: use column 7
    For bed12: use column 13

    This method uses the maximum number of reads found in any interval as the tag count.

    Tiles with no counts will not be output.
    '''
    
    to_cluster = True

    src = " ".join( [ '''<( zcat %s | awk '{printf("%%s:%%i-%%i\\t%%i\\n", $1,$2,$3,$4 );}' ) ''' % x for x in infiles] )
    tmpfile = P.getTempFilename( "." )
    statement = '''paste %(src)s > %(tmpfile)s'''
    P.run()
    
    tracks = [ re.sub( "\..*", '', os.path.basename(x) ) for x in infiles ]

    outf = IOTools.openFile( outfile, "w")
    outf.write( "interval_id\t%s\n" % "\t".join( tracks ) )
    
    for line in open( tmpfile, "r" ):
        data = line[:-1].split("\t")
        genes = list(set([ data[x] for x in range(0,len(data), 2 ) ]))
        values = [ int(data[x]) for x in range(1,len(data), 2 ) ]
        if sum(values) == 0: continue
        assert len(genes) == 1, "paste command failed, wrong number of genes per line: '%s'" % line
        outf.write( "%s\t%s\n" % (genes[0], "\t".join(map(str, values) ) ) )
    
    outf.close()

    os.unlink(tmpfile)


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
                            "up", "down",
                            "twofold",
                            "twofold_up", "twofold_down",
                            ) ) + "\n" )

    all_tables = set(Database.getTables( dbhandle ))
    outdir = os.path.join( PARAMS["exportdir"], "diff_methylation" )

    for tablename in tables:

        prefix = P.snip( tablename, "_%s" % method )
        tileset, design = prefix.split("_")

        def toDict( vals, l = 2 ):
            return collections.defaultdict( int, [ (tuple( x[:l]), x[l]) for x in vals ] )

        E.info( "collecting data from %s" % tablename )
        
        tested = toDict( Database.executewait( dbhandle,
                                               """SELECT treatment_name, control_name, COUNT(*) FROM %(tablename)s 
                                GROUP BY treatment_name,control_name""" % locals() ).fetchall() )
        status = toDict( Database.executewait( dbhandle,
                                               """SELECT treatment_name, control_name, status, COUNT(*) FROM %(tablename)s 
                                GROUP BY treatment_name,control_name,status""" % locals() ).fetchall(), 3 )
        signif = toDict( Database.executewait( dbhandle,
                                               """SELECT treatment_name, control_name, COUNT(*) FROM %(tablename)s 
                                WHERE significant
                                GROUP BY treatment_name,control_name""" % locals() ).fetchall() )
        fold2 = toDict( Database.executewait( dbhandle,
                """SELECT treatment_name, control_name, COUNT(*) FROM %(tablename)s 
                                WHERE (l2fold >= 1 or l2fold <= -1) AND significant
                                GROUP BY treatment_name,control_name,significant""" % locals() ).fetchall() )

        up = toDict( Database.executewait( dbhandle,
                                                """SELECT treatment_name, control_name, COUNT(*) FROM %(tablename)s 
                                WHERE l2fold > 0 AND significant
                                GROUP BY treatment_name,control_name,significant""" % locals() ).fetchall() )

        down = toDict( Database.executewait( dbhandle,
                                             """SELECT treatment_name, control_name, COUNT(*) FROM %(tablename)s 
                                WHERE l2fold < 0 AND significant
                                GROUP BY treatment_name,control_name,significant""" % locals() ).fetchall() )

        fold2up = toDict( Database.executewait( dbhandle,
                                           """SELECT treatment_name, control_name, COUNT(*) FROM %(tablename)s 
                                WHERE l2fold > 1 AND significant
                                GROUP BY treatment_name,control_name,significant""" % locals() ).fetchall() )

        fold2down = toDict( Database.executewait( dbhandle,
                                             """SELECT treatment_name, control_name, COUNT(*) FROM %(tablename)s 
                                WHERE l2fold < -1 AND significant
                                GROUP BY treatment_name,control_name,significant""" % locals() ).fetchall() )
        
        groups = tested.keys()

        for treatment_name, control_name in groups:
            k = (treatment_name,control_name)
            outf.write( "\t".join(map(str, (
                            tileset,
                            design,
                            treatment_name,
                            control_name,
                            tested[k],
                            "\t".join( [ str(status[(treatment_name,control_name,x)]) for x in keys_status]),
                            signif[(k)],
                            up[k], down[k],
                            fold2[k],
                            fold2up[k], fold2down[k] ) ) ) + "\n" )
                            

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

