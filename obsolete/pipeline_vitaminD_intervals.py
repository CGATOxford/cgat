'''ChIP-Seq tasks associated with intervals.
'''

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
import sqlite3
import cStringIO

import CGAT.Experiment as E
import CGAT.Pipeline as P

import csv
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IndexedGenome as IndexedGenome
import CGAT.FastaIterator as FastaIterator
import CGAT.Genomics as Genomics
import CGAT.IOTools as IOTools
import CGAT.MAST as MAST
import CGAT.GTF as GTF
import CGAT.GFF as GFF
import CGAT.Bed as Bed
import CGAT.Stats as Stats

import pysam
import numpy
import gzip
import fileinput

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################
P.getParameters( 
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini" ] )

PARAMS = P.PARAMS

if os.path.exists("conf.py"):
    E.info( "reading additional configuration from conf.py" )
    execfile("conf.py")

############################################################
############################################################
############################################################
## 
############################################################
def getPeakShift( infile ):
    '''get peak shift for filename infile (.macs output file).

    returns None if no shift found'''

    shift = None
    with open(infile, "r") as ins:
        rx = re.compile("#2   d: (\d+)")
        for line in ins:
            x = rx.search(line)
            if x: 
                shift = int(x.groups()[0])
                break
    return shift

############################################################
############################################################
############################################################
def getMappedReads( infile ):
    '''return number of reads mapped.
    '''
    for lines in open(infile,"r"):
        data = lines[:-1].split("\t")
        if data[1].startswith( "without duplicates"):
            return int(data[0])
    return

############################################################
############################################################
############################################################
def getMinimumMappedReads( infiles ):
    '''find the minimum number of mapped reads in infiles.'''
    v = []
    for infile in infiles:
        x = getMappedReads( infile )
        if x: v.append( x )
    if len(v) == 0:
        raise P.PipelineError( "could not find mapped reads in files %s" % (str(infiles)))
    return min(v)
    
############################################################
############################################################
############################################################
##
############################################################
def buildNormalizedBAM( infiles, outfile, normalize = True ):
    '''build a normalized BAM file.

    Infiles are merged and duplicated reads are removed. 
    If *normalize* is set, reads are removed such that all 
    files will have approximately the same number of reads.
    '''

    min_reads = getMinimumMappedReads( glob.glob("*.readstats") )
    
    samfiles = []
    num_reads = 0
    for infile, statsfile in infiles:
        samfiles.append( pysam.Samfile( infile, "rb" ) )
        num_reads += getMappedReads( statsfile )

    threshold = float(min_reads) / num_reads 

    pysam_out = pysam.Samfile( outfile, "wb", template = samfiles[0] )

    ninput, noutput, nduplicates = 0, 0, 0

    # iterate over mapped reads
    last_contig, last_pos = None, None
    for pysam_in in samfiles:
        for read in pysam_in.fetch():

            ninput += 1
            if read.rname == last_contig and read.pos == last_pos:
                nduplicates += 1
                continue

            if normalize and random.random() <= threshold:
                pysam_out.write( read )
                noutput += 1

            last_contig, last_pos = read.rname, read.pos

        pysam_in.close()

    pysam_out.close()

    logs = open( outfile + ".log", "w")
    logs.write("# min_reads=%i, threshold= %5.2f\n" % \
                   (min_reads, threshold))
    logs.write("set\tcounts\tpercent\n")
    logs.write("ninput\t%i\t%5.2f%%\n" % (ninput, 100.0) )
    nwithout_dups = ninput - nduplicates
    logs.write("duplicates\t%i\t%5.2f%%\n" % (nduplicates,100.0*nduplicates/ninput))
    logs.write("without duplicates\t%i\t%5.2f%%\n" % (nwithout_dups,100.0*nwithout_dups/ninput))
    logs.write("target\t%i\t%5.2f%%\n" %   (min_reads,100.0*min_reads/nwithout_dups))
    logs.write("noutput\t%i\t%5.2f%%\n" % (noutput,100.0*noutput/nwithout_dups))
    
    logs.close()
    
    # if more than one samfile: sort
    if len(samfiles) > 1:
        tmpfilename = P.getTempFilename()
        pysam.sort( outfile, tmpfilename )
        shutil.move( tmpfilename + ".bam", outfile )
        os.unlink( tmpfilename )

    pysam.index( outfile )

    P.info( "buildNormalizedBam: %i input, %i output (%5.2f%%), should be %i" % (ninput, noutput, 100.0*noutput/ninput, min_reads ))

############################################################
############################################################
############################################################
def buildBAMStats( infile, outfile ):
    '''calculate bamfile statistics
    '''

    # no bedToBigBed
    # to_cluster = True
    outs = open(outfile, "w" )
    outs.write( "reads\tcategory\n" )
    for line in pysam.flagstat( infile ):
        data = line[:-1].split( " ")
        outs.write( "%s\t%s\n" % (data[0], " ".join(data[1:]) ) )
    pysam_in = pysam.Samfile( infile, "rb" )


    outs_dupl = open( outfile + ".duplicates", "w" )
    outs_dupl.write( "contig\tpos\tcounts\n" )

    outs_hist = open( outfile + ".histogram", "w" )
    outs_hist.write( "duplicates\tcounts\tcumul\tfreq\tcumul_freq\n" )

    last_contig, last_pos = None, None
    ninput, nduplicates = 0, 0

    duplicates = collections.defaultdict( int )
    counts = collections.defaultdict( int )
    count = 0

    for read in pysam_in.fetch():

        ninput += 1

        if read.rname == last_contig and read.pos == last_pos:
            count += 1
            nduplicates += 1
            continue

        if count > 1:
            outs_dupl.write("%s\t%i\t%i\n" % (last_contig, last_pos, count) )
            counts[count] += 1

        count = 1
        last_contig, last_pos = read.rname, read.pos

    outs.write("%i\tduplicates (%5.2f%%)\n" % (nduplicates, 100.0* nduplicates / ninput))
    outs.write("%i\twithout duplicates (%5.2f%%)\n" % (ninput - nduplicates,
                                                       100.0*(ninput - nduplicates)/ninput))
    pysam_in.close()
    outs.close()
    outs_dupl.close()

    keys = counts.keys()
    # count per position (not the same as nduplicates, which is # of reads)
    c = 0
    total = sum( counts.values() )
    for k in sorted(keys):
        c += counts[k]
        outs_hist.write("%i\t%i\t%i\t%f\t%f\n" % (k, counts[k], c, 
                                                  100.0 * counts[k] / total,
                                                  100.0 * c / total) )
    outs_hist.close()
    

############################################################
############################################################
############################################################
def exportIntervalsAsBed( infile, outfile ):
    '''export sequences for all intervals.'''

    dbhandle = sqlite3.connect( PARAMS["database"] )

    if infile.endswith("_macs.import"):
        track = infile[:-len("_macs.import")]
    else:
        track = infile[:-len("_intervals.import")]
        
    if track.startswith("control"): return

    cc = dbhandle.cursor()
    statement = "SELECT contig, start, end, interval_id, peakval FROM %s_intervals ORDER by contig, start" % track
    cc.execute( statement )

    outs = open( outfile, "w")

    for result in cc:
        contig, start, end, interval_id,peakval = result
        peakval = int(min(peakval,1000))
        outs.write( "%s\t%i\t%i\t%s\t%i\n" % (contig, start, end, str(interval_id), peakval) )

    cc.close()
    outs.close()

############################################################
############################################################
############################################################
def exportPeaksAsBed( infile, outfile ):
    '''export sequences for all peaks.'''

    dbhandle = sqlite3.connect( PARAMS["database"] )

    if infile.endswith("_macs.import"):
        track = infile[:-len("_macs.import")]
    else:
        track = infile[:-len("_intervals.import")]
        
    if track.startswith("control"): return
    
    peakwidth = PARAMS["peakwidth"]
    
    cc = dbhandle.cursor()
    statement = '''SELECT contig, peakcenter - %(peakwidth)i, peakcenter + %(peakwidth)i,
                          interval_id, peakval FROM %(track)s_intervals ORDER by contig, start''' % locals()
    cc.execute( statement )

    outs = open( outfile, "w")

    for result in cc:
        contig, start, end, interval_id,peakval = result
        peakval = int(min(peakval,1000))
        outs.write( "%s\t%i\t%i\t%s\t%i\n" % (contig, start, end, str(interval_id), peakval) )

    cc.close()
    outs.close()

############################################################
############################################################
############################################################
def importCombinedIntervals( infile, outfile ):
    '''import combined intervals.

    Also, re-evaluate the intervals by counting reads within
    the interval. In contrast to the initial pipeline, the
    genome is not binned. In particular, the meaning of the
    columns in the table changes to:

    nProbes: number of reads in interval
    PeakCenter: position with maximum number of reads in interval
    AvgVal: average coverage within interval

    If *replicates* is true, only replicates will be considered
    for the counting. Otherwise the counts aggregate both replicates
    and conditions.
    '''

    tmpfile = tempfile.NamedTemporaryFile(delete=False)

    headers = ("AvgVal","DisttoStart","GeneList","Length","PeakCenter","PeakVal","Position","interval_id","nCpGs","nGenes","nPeaks","nProbes","nPromoters", "contig","start","end" )

    tmpfile.write( "\t".join(headers) + "\n" )

    avgval,contig,disttostart,end,genelist,length,peakcenter,peakval,position,start,interval_id,ncpgs,ngenes,npeaks,nprobes,npromoters = \
        0,"",0,0,"",0,0,0,0,0,0,0,0,0,0,0,

    samfiles, offsets = [], []

    track = infile[:-len(".bed")]
    cellline, condition, replicate = splitTrack( track )
    
    if cellline and condition:
        for replicate in REPLICATES:
            fn = FILEPATTERN % locals()
            samfiles.append( pysam.Samfile( fn + ".norm.bam",  "rb" ) )
            if os.path.exists( fn + ".macs" ):
                offsets.append( getPeakShift( fn + ".macs" ) )

    elif condition:
        for cellline in CELLLINES:
            for replicate in REPLICATES:
                fn = FILEPATTERN % locals()
                samfiles.append( pysam.Samfile( fn + ".norm.bam",  "rb" ) )
                if os.path.exists( fn + ".macs" ):
                    offsets.append( getPeakShift( fn + ".macs" ) )

    mlength = int(PARAMS["merge_min_interval_length"])
    for line in open(infile, "r"):
        contig, start, end, interval_id = line[:-1].split()[:4]
        start, end = int(start), int(end)
        
        # remove very short intervals
        if end-start < mlength: continue

        npeaks, peakcenter, length, avgval, peakval, nprobes = countPeaks( contig, start, end, samfiles, offsets )

        tmpfile.write( "\t".join( map( str, (avgval,disttostart,genelist,length,peakcenter,peakval,position,interval_id,ncpgs,ngenes,npeaks,nprobes,npromoters, contig,start,end) )) + "\n" )

    tmpfile.close()

    tmpfilename = tmpfile.name
    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=interval_id \
              --table=%(track)s_intervals \
    < %(tmpfilename)s > %(outfile)s
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )
    os.unlink( tmpfile.name )

############################################################
############################################################
############################################################
def mergeBedFiles( infiles, outfile ):
    '''generic method for merging bed files.
    '''

    if len(infiles) < 2:
        raise ValueError( "expected at least two files to merge into %s" % outfile )

    infile = " ".join( infiles )
    statement = '''
        cat %(infile)s |\
        mergeBed -i stdin |\
        cut -f 1-3 |\
        awk '{printf("%%s\\t%%i\\n",$0, ++a); }'
        > %(outfile)s 
        ''' 

    P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
def intersectBedFiles( infiles, outfile ):
    '''generic method for merging bed files.
    '''

    if len(infiles) != 2:
        raise ValueError( "expected two files to merge into %s" % outfile )

    statement = '''
        intersectBed -a %s -b %s |\
        cut -f 1,2,3,4,5 |\
        awk 'BEGIN { OFS="\\t"; } {$4=++a; print;}'
        > %%(outfile)s 
        ''' % (infiles[0], infiles[1])

    P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
def subtractBedFiles( infile, subtractfile, outfile ):
    '''subtract intervals in *subtractfile* from *infile*
    and store in *outfile*.
    '''
    statement = '''
        intersectBed -v -a %(infile)s -b %(subtractfile)s |\
        cut -f 1,2,3,4,5 |\
        awk 'BEGIN { OFS="\\t"; } {$4=++a; print;}'
        > %%(outfile)s 
        ''' % locals()

    P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
def summarizeMACS( infiles, outfile ):
    '''run MACS for peak detection.'''
    def __get( line, stmt ):
        x = line.search(stmt )
        if x: return x.groups() 

    map_targets = [
        ("unique tags in treatment: (\d+)", "tag_treatment_unique",()),
        ("total tags in treatment: (\d+)", "tag_treatment_total",()),
        ("unique tags in control: (\d+)", "tag_control_unique",()),
        ("total tags in control: (\d+)", "tag_control_total",()),
        ("#2 number of paired peaks: (\d+)", "paired_peaks",()),
        ("#2   min_tags: (\d+)","min_tags", ()),
        ("#2   d: (\d+)", "shift", ()),
        ("#2   scan_window: (\d+)", "scan_window", ()),
        ("#3 Total number of candidates: (\d+)", "ncandidates",("positive", "negative") ),
        ("#3 Finally, (\d+) peaks are called!",  "called", ("positive", "negative") ) ]

    mapper, mapper_header = {}, {}
    for x,y,z in map_targets: 
        mapper[y] = re.compile( x )
        mapper_header[y] = z

    keys = [ x[1] for x in map_targets ]

    outs = open(outfile,"w")

    headers = []
    for k in keys:
        if mapper_header[k]:
            headers.extend( ["%s_%s" % (k,x) for x in mapper_header[k] ])
        else:
            headers.append( k )
    outs.write("track\t%s" % "\t".join(headers) + "\n" )

    for infile in infiles:
        results = collections.defaultdict(list)
        with open( infile ) as f:
            for line in f:
                if "diag:" in line: break
                for x,y in mapper.items():
                    s = y.search( line )
                    if s: 
                        results[x].append( s.groups()[0] )
                        break
                
        row = [ infile[:-len(".macs")] ]
        for key in keys:
            val = results[key]
            if len(val) == 0: v = "na"
            else: 
                c = len(mapper_header[key])
                if c >= 1: assert len(val) == c, "key=%s, expected=%i, got=%i, val=%s, c=%s" %\
                   (key,
                    len(val),
                    c,
                    str(val), mapper_header[key])
                v = "\t".join( val )
            row.append(v)
        outs.write("\t".join(row) + "\n" )

    outs.close()

###################################################################
###################################################################
###################################################################
##
###################################################################
def importMACSSummary( infile, outfile ):
    '''import regions of interest.'''
    
    table = outfile[:-len(".import")]

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=track \
              --table=%(table)s \
    < %(infile)s > %(outfile)s
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
def importMACS( infile, outfile, suffix = ".norm.bam" ):
    '''import MACS results.

    Imports only positive peaks. It filters peaks by p-value,
    q-value and fold change and imports the diagnostic data.

    Does re-counting of peakcenter, peakval, ... using
    the normalized tag counts.
    '''

    track = infile[:-len(".macs")]    
    infilename = infile + "_peaks.xls"
    filename_diag = infile + "_diag.xls"
    filename_r = infile + "_model.r"
    
    if not os.path.exists(infilename):
        E.warn("could not find %s" % infilename )
        outs = open(outfile,"w")
        outs.close()
        return

    shift = getPeakShift( infile )
    assert shift != None, "could not determine peak shift from MACS file %s" % infile

    samfiles = [ pysam.Samfile( track + suffix, "rb" ) ]
    offsets = [ shift / 2 ]
    outtemp = P.getTempFile()

    outtemp.write( "\t".join( ( \
                "interval_id", 
                "contig", "start", "end",
                "npeaks", "peakcenter", 
                "length", 
                "avgval", "peakval",
                "nprobes",
                "pvalue", "fold", "qvalue",
                "macs_summit", "macs_nprobes",
                )) + "\n" )
    id = 0

    ## get thresholds
    max_qvalue = float(PARAMS["macs_max_qvalue"])
    # min, as it is -10log10
    min_pvalue = float(PARAMS["macs_min_pvalue"])
    min_fold = float(PARAMS["macs_min_fold"])
    
    counter = E.Counter()
    with open( infilename, "r" ) as ins:
        for line in ins:
            if line.startswith("#"): continue
            if line.startswith( "chr\tstart"): continue
            counter.input += 1
            data = line[:-1].split("\t")
            if len(data) == 9:
                contig,start,end,length,summit,ntags,pvalue,fold,qvalue = data
            elif len(data) == 8:
                contig,start,end,length,summit,ntags,pvalue,fold = data
                qvalue = 1.0
            else:
                raise ValueError( "could not parse line %s" % line )
            
            pvalue, qvalue, summit, fold = float(pvalue), float(qvalue), int(summit), float(fold)

            if qvalue > max_qvalue or pvalue < min_pvalue or fold < min_fold: 
                counter.skipped += 1
                continue

            # these are 1-based coordinates
            start, end = int(start)-1, int(end)
            assert start < end
            # macs can have negative start coordinates
            start = max(start, 0)
            npeaks, peakcenter, length, avgval, peakval, nreads = countPeaks( contig, start, end, samfiles, offsets )

            outtemp.write ( "\t".join( map(str, ( \
                            id, contig, start, end, npeaks, peakcenter, length, avgval, peakval, nreads,
                            pvalue, fold, qvalue,
                            start + summit - 1, 
                            ntags) ) ) + "\n" )
            id += 1                        
            counter.output += 1
    outtemp.close()

    tablename = "%s_intervals" % track
    tmpfilename = outtemp.name

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=interval_id \
              --table=%(tablename)s \
    < %(tmpfilename)s > %(outfile)s
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )

    # import diagnotic data
    if os.path.exists( filename_diag ):

        tablename = "%s_macsdiag" % track

        statement = '''
        sed "s/FC range.*/fc\\tnpeaks\\tp90\\tp80\\tp70\\tp60\\tp50\\tp40\\tp30\\tp20/" < %(filename_diag)s |\
       python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
                  --map=fc:str \
                  --table=%(tablename)s \
        > %(outfile)s
        '''

        P.run( **dict( locals().items() + PARAMS.items() ) )


    # create plot
    if os.path.exists( filename_r ):

        target_path = os.path.join( os.getcwd(), "doc", "_static", "MACS" )
        try:
            os.makedirs( target_path )
        except OSError: 
            # ignore "file exists" exception
            pass

        statement = '''
        R --vanilla < %(track)s.macs_model.r > %(outfile)s
        '''
        
        P.run( **dict( locals().items() + PARAMS.items() ) )

        shutil.copyfile(
            "%s.macs_model.pdf" % track,
            os.path.join( target_path, "%s_model.pdf" % track) )
        
    os.unlink( tmpfilename )

    E.info("%s: %s" % (track, str(counter))) 

############################################################
############################################################
############################################################
def runMACS( infile, outfile ):
    '''run MACS for peak detection.

    The output bed files contain the P-value as their score field.
    '''
    to_cluster = True

    if infile.endswith( ".norm.bam"):

        track = infile[:-len(".norm.bam")]
        if track.startswith("control"):
            P.touch( outfile )
            return

        format = "bam"
        suffix = ".norm.bam"

    elif infile.endswith( ".bam"):

        track = infile[:-len(".bam")]
        if track.startswith("control"):
            P.touch( outfile )
            return

        format = "bam"
        suffix = ".norm.bam"
        
    elif infile.endswith(".bed.gz"):

        track = infile[:-len(".bed.gz")]
        if track.startswith("control"):
            outs = open( outfile, "w")
            outs.close()
            return

        format = "bed"
        suffix = ".bed.gz"
        
    control = getControl( track )
    
    if control != None:
        control += suffix
    else:
        E.info("%s: no control for track %s" % (outfile, track ) )
        control = None
            
    if control: control = "-c %s" % control
    else: control = ""


        
    statement = '''
    macs -t %(infile)s %(control)s \
    --diag \
    --name=%(outfile)s \
    --format=%(format)s \
    %(macs_options)s >& %(outfile)s''' 
    
    P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
def getCounts( contig, start, end, samfiles, offsets = [] ):
    '''count reads per position.'''
    assert len(offsets) == 0 or len(samfiles) == len(offsets)

    length = end - start
    counts = numpy.zeros( length )

    nreads = 0

    if offsets:
        # if offsets are given, shift tags. 
        for samfile, offset in zip(samfiles,offsets):

            # for peak counting I follow the MACS protocoll,
            # see the function def __tags_call_peak in PeakDetect.py
            # In words
            # Only take the start of reads (taking into account the strand)
            # add d/2=offset to each side of peak and start accumulate counts.

            # for counting, extend reads by offset
            # on + strand shift tags upstream
            # i.e. look at the downstream window
            xstart, xend = max(0, start - offset), max(0, end - offset)
            for read in samfile.fetch( contig, xstart, xend ):
                if read.is_reverse: continue
                nreads += 1
                rstart = max( 0, read.pos - xstart - offset)
                rend = min( length, read.pos - xstart + offset) 
                counts[ rstart:rend ] += 1

            # on the - strand, shift tags downstream
            xstart, xend = max(0, start + offset), max(0, end + offset)
            for read in samfile.fetch( contig, xstart, xend ):
                if not read.is_reverse: continue
                nreads += 1
                rstart = max( 0, read.pos + read.rlen - xstart - offset)
                rend = min( length, read.pos + read.rlen - xstart + offset) 
                counts[ rstart:rend ] += 1
    else:
        for samfile in samfiles:
            for read in samfile.fetch( contig, start, end ):
                nreads += 1
                rstart = max( 0, read.pos - start )
                rend = min( length, read.pos - start + read.rlen ) 
                counts[ rstart:rend ] += 1
    return nreads, counts

############################################################
############################################################
############################################################
def countPeaks( contig, start, end, samfiles, offsets = None):
    '''update peak values within interval contig:start-end.

    If offsets is given, tags are moved by the offset
    before summarizing.
    '''

    nreads, counts = getCounts( contig, start, end, samfiles, offsets )

    # nreads can be 0 if the intervals overlap only slightly
    # and due to the binning, no reads are actually in the overlap region.
    # However, these intervals should be small and have already be deleted via 
    # the merge_min_interval_length cutoff, so I keep the assertion.
    assert nreads > 0, "no reads in interval %s:%i-%i" % (contig, start, end)

    length = end - start            
    nprobes = nreads
    avgval = numpy.mean( counts )
    peakval = max(counts)

    # set other peak parameters
    peaks = numpy.array( range(0,length) )[ counts >= peakval ]
    npeaks = len( peaks )
    # peakcenter is median coordinate between peaks
    # such that it is a valid peak in the middle
    peakcenter = start + peaks[npeaks//2] 

    return npeaks, peakcenter, length, avgval, peakval, nreads

############################################################
############################################################
############################################################
##
############################################################
def makeIntervalCorrelation( infiles, outfile, field ):
    '''compute correlation of interval properties between sets
    '''

    dbhandle = sqlite3.connect( PARAMS["database"] )

    tracks, idx = [], []
    for infile in infiles:
        track = infile[:-len(".bed")]
        cc = dbhandle.cursor()
        statement = "SELECT contig, start, end, %(field)s FROM %(track)s_intervals" % locals()
        cc.execute( statement )
        ix = IndexedGenome.IndexedGenome()
        for contig, start, end, peakval in cc:
            ix.add( contig, start, end, peakval )        
        idx.append( ix )
        tracks.append( track )
    outs = open( outfile, "w" )
    outs.write( "id\tcontig\tstart\tend\t" + "\t".join( tracks ) + "\n" )

    for bed in Bed.iterator( infile = open( "merged.bed", "r") ):
        
        row = []
        for ix in idx:
            try:
                intervals = list(ix.get( bed.contig, bed.start, bed.end ))
            except KeyError:
                row.append( "" )
                continue
        
            if len(intervals) == 0:
                peakval = ""
            else:
                peakval = str( (max( [ x[2] for x in intervals ] )) )
            row.append( peakval )

        outs.write( str(bed) + "\t" + "\t".join( row ) + "\n" )

    outs.close()

