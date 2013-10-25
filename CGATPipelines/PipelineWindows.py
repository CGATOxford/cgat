import os
import re
import sys
import collections

import CGAT.Experiment as E
import CGAT.Pipeline as P
import CGAT.IOTools as IOTools

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
def buildDMRStats( infile, outfile, method ):
    '''build dmr summary statistics.
    '''
    results = collections.defaultdict( lambda : collections.defaultdict(int) )
    status =  collections.defaultdict( lambda : collections.defaultdict(int) )
    x = 0
    for line in IOTools.iterate( IOTools.openFile( infile ) ):
        key = (line.treatment_name, line.control_name )
        r,s = results[key], status[key]
        r["tested"] += 1
        s[line.status] += 1

        is_significant = line.significant == "1"
        up = float(line.l2fold) > 0
        down = float(line.l2fold) < 0
        fold2up = float(line.l2fold) > 1
        fold2down = float(line.l2fold) < -1
        fold2 = fold2up or fold2down

        if up: r["up"] += 1
        if down: r["down"] += 1

        if is_significant:
            r["significant"] += 1
            if up: r["significant_up"] += 1
            if down: r["significant_down"] += 1
            if fold2: r["fold2"] += 1
            if fold2up: r["l2fold_up"] += 1
            if fold2down: r["l2fold_down"] += 1
            
    header1, header2 = set(), set()
    for r in results.values(): header1.update( r.keys() )
    for s in status.values(): header2.update( s.keys() )
    
    header = ["method", "treatment", "control" ]
    header1 = list(header1)
    header2 = list(header2)

    outf = IOTools.openFile( outfile, "w" )
    outf.write( "\t".join(header + header1 + header2) + "\n" )

    for treatment,control in results.keys():
        key = (treatment,control)
        r = results[key]
        s = status[key]
        outf.write( "%s\t%s\t%s\t" % (method,treatment, control))
        outf.write( "\t".join( [str(r[x]) for x in header1 ] ) + "\t" )
        outf.write( "\t".join( [str(s[x]) for x in header2 ] ) + "\n" )


        # ###########################################
        # ###########################################
        # ###########################################
        # # plot length versus P-Value
        # data = Database.executewait( dbhandle, 
        #                              '''SELECT end - start, pvalue 
        #                      FROM %(tablename)s
        #                      WHERE significant'''% locals() ).fetchall()

        # # require at least 10 datapoints - otherwise smooth scatter fails
        # if len(data) > 10:
        #     data = zip(*data)

        #     pngfile = "%(outdir)s/%(tileset)s_%(design)s_%(method)s_pvalue_vs_length.png" % locals()
        #     R.png( pngfile )
        #     R.smoothScatter( R.log10( ro.FloatVector(data[0]) ),
        #                      R.log10( ro.FloatVector(data[1]) ),
        #                      xlab = 'log10( length )',
        #                      ylab = 'log10( pvalue )',
        #                      log="x", pch=20, cex=.1 )

        #     R['dev.off']()

#    outf.close()

