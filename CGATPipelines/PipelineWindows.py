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

