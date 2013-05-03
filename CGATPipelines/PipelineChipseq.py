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
import CGAT.GTF as GTF
import CGAT.GFF as GFF
import CGAT.Bed as Bed
import MACS
import CGAT.WrapperZinba as WrapperZinba
# import Stats

import CGATPipelines.PipelineMapping as PipelineMapping
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
    ["%s.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini" ] )

PARAMS = P.PARAMS

if os.path.exists("pipeline_conf.py"):
    E.info( "reading additional configuration from pipeline_conf.py" )
    execfile("pipeline_conf.py")

############################################################
############################################################
############################################################
## 
############################################################
def getPeakShiftFromMacs( infile ):
    '''get peak shift for filename infile (.macs output file).

    returns None if no shift found'''

    shift = None
    with IOTools.openFile(infile, "r") as ins:
        rx = re.compile("#2 predicted fragment length is (\d+) bps")
        r2 = re.compile("#2 Use (\d+) as shiftsize, \d+ as fragment length" )
        for line in ins:
            x = rx.search(line)
            if x: 
                shift = int(x.groups()[0])
                break
            x = r2.search(line)
            if x: 
                shift = int(x.groups()[0])
                E.warn( "shift size was set automatically - see MACS logfiles" )
                break
            
    return shift

############################################################
############################################################
############################################################
## 
############################################################
def getPeakShiftFromZinba( infile ):
    '''get peak shift for filename infile (.zinba output file).

    returns None if no shift found
    '''

    shift = None

    # search for
    # $offset
    # [1] 125

    with IOTools.openFile(infile, "r") as ins:
        lines = ins.readlines()
        for i, line in enumerate(lines):
            if line.startswith("$offset"):
                shift = int(lines[i+1].split()[1])
                break
            
    return shift

def getPeakShift( track ):
    '''get peak shift for a track.'''

    if os.path.exists( "%s.macs" % track ):
        return getPeakShiftFromMacs( "%s.macs" % track )
    elif os.path.exists( "%s.zinba" % track ):
        return getPeakShiftFromZinba( "%s.zinba" % track )
    
############################################################
############################################################
############################################################
def getMappedReads( infile ):
    '''return number of reads mapped. '''
    for lines in IOTools.openFile(infile,"r"):
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
def getExonLocations(filename):
    '''return a list of exon locations as Bed entries
    from a file contain a one ensembl gene ID per line
    '''
    fh = IOTools.openFile(filename,"r")
    ensembl_ids = []
    for line in fh:
        ensembl_ids.append(line.strip())
    fh.close()

    dbhandle = sqlite3.connect(PARAMS["annotations_database"])
    cc = dbhandle.cursor()

    gene_ids = []
    n_ids = 0
    for ID in ensembl_ids:
        gene_ids.append('gene_id="%s"' % ID)
        n_ids += 1

    statement = "select contig,start,end from geneset_cds_gtf where "+" OR ".join(gene_ids)
    
    cc.execute( statement)

    region_list = []
    n_regions = 0
    for result in cc:
        b = Bed.Bed()
        b.contig, b.start, b.end = result
        region_list.append( b )
        n_regions +=1

    cc.close()

    E.info("Retrieved exon locations for %i genes. Got %i regions" % (n_ids,n_regions) )

    return(region_list)

############################################################
############################################################
############################################################
def getBedLocations(filename):
    '''return a list of regions as (contig,start,end) tuples
    from a bed file'''
    fh = open(filename,"r")
    region_list = []
    n_regions = 0

    for line in fh:
        if line.strip() != "" and line[0] !="#":
            fields = line.split("\t")
            contig, start, end = fields[0], int(fields[1]), int(fields[2])
            region_list.append((contig,start,end))
            n_regions +=1

    fh.close()

    #E.info("Read in %i regions from %s" % ( n_regions, filename) )
    return (region_list)

############################################################
############################################################
############################################################
def buildQuicksectMask(bed_file):
    '''return Quicksect object containing the regions specified
       takes a bed file listing the regions to mask 
    '''
    mask = IndexedGenome.Quicksect()

    n_regions = 0
    for bed in Bed.iterator( IOTools.openFile( bed_file ) ):
        #it is neccessary to extend the region to make an accurate mask
        mask.add(bed.contig,(bed.start-1),(bed.end+1),1)
        n_regions += 1

    E.info("Built Quicksect mask for %i regions" % n_regions)

    return(mask)

############################################################
############################################################
############################################################
def buildBAMforPeakCalling( infiles, outfile, dedup, mask):
    ''' Make a BAM file suitable for peak calling.

        Infiles are merged and unmapped reads removed. 

        If specificied duplicate reads are removed. 
        This method use Picard.

        If a mask is specified, reads falling within
        the mask are filtered out. 

        This uses bedtools.

        The mask is a quicksect object containing
        the regions from which reads are to be excluded.
    '''

    #open the infiles, if more than one merge and sort first using samtools.

    samfiles = []
    num_reads = 0
    nfiles = 0

    statement = []
    
    tmpfile = P.getTempFilename()

    if len(infiles) > 1 and isinstance(infiles,str)==0:
        # assume: samtools merge output is sorted
        # assume: sam files are sorted already
        statement.append( '''samtools merge @OUT@ %s''' % (infiles.join(" ")) )
        statement.append( '''samtools sort @IN@ @OUT@''')

    if dedup:
        statement.append( '''MarkDuplicates
                                       INPUT=@IN@
                                       ASSUME_SORTED=true 
                                       REMOVE_DUPLICATES=true
                                       QUIET=true
                                       OUTPUT=@OUT@
                                       METRICS_FILE=%(outfile)s.picardmetrics
                                       VALIDATION_STRINGENCY=SILENT 
                   > %(outfile)s.picardlog ''' )

    if mask:
        statement.append( '''intersectBed -abam @IN@ -b %(mask)s -wa -v > @OUT@''' )

    statement.append('''mv @IN@ %(outfile)s''' )
    statement.append('''samtools index %(outfile)s''' )

    statement = P.joinStatements( statement, infiles )
    P.run()

############################################################
############################################################
############################################################
def buildSimpleNormalizedBAM( infiles, outfile, nreads ):
    '''normalize a bam file to given number of counts
       by random sampling
    '''
    infile,countfile = infiles

    pysam_in = pysam.Samfile (infile,"rb")
    
    fh = IOTools.openFile(countfile,"r")
    readcount = int(fh.read())
    fh.close()
    
    threshold = float(nreads) / float(readcount)

    pysam_out = pysam.Samfile( outfile, "wb", template = pysam_in )

    # iterate over mapped reads thinning by the threshold
    ninput, noutput = 0,0
    for read in pysam_in.fetch():
         ninput += 1
         if random.random() <= threshold:
             pysam_out.write( read )
             noutput += 1

    pysam_in.close()
    pysam_out.close()
    pysam.index( outfile )

    E.info( "buildNormalizedBam: %i input, %i output (%5.2f%%), should be %i" % (ninput, noutput, 100.0*noutput/ninput, nreads ))

############################################################                                                        
############################################################                                                        
############################################################                                                        
####### Depreciate this function? ##########################
############################################################                                                        
def buildNormalizedBAM( infiles, outfile, normalize = True ):
    '''build a normalized BAM file.

    Infiles are merged and duplicated reads are removed. 
    If *normalize* is set, reads are removed such that all 
    files will have approximately the same number of reads.

    Note that the duplication here is wrong as there
    is no sense of strandedness preserved.
    '''

    min_reads = getMinimumMappedReads( glob.glob("*.readstats") )
    
    samfiles = []
    num_reads = 0
    for infile, statsfile in infiles:
        samfiles.append( pysam.Samfile( infile, "rb" ) )
        num_reads += getMappedReads( statsfile )

    threshold = float(min_reads) / num_reads 

    E.info( "%s: min reads: %i, total reads=%i, threshold=%f" % (infiles, min_reads, num_reads, threshold) )

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

    logs = IOTools.openFile( outfile + ".log", "w")
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

    E.info( "buildNormalizedBam: %i input, %i output (%5.2f%%), should be %i" % (ninput, noutput, 100.0*noutput/ninput, min_reads ))


############################################################
############################################################
############################################################
def buildBAMStats( infile, outfile ):
    '''calculate bamfile statistics - currently only single-ended
    duplicates.
    '''

    # no bedToBigBed
    # to_cluster = True
    outs = IOTools.openFile(outfile, "w" )
    outs.write( "reads\tcategory\n" )
    for line in pysam.flagstat( infile ):
        data = line[:-1].split( " ")
        outs.write( "%s\t%s\n" % (data[0], " ".join(data[1:]) ) )

    pysam_in = pysam.Samfile( infile, "rb" )

    outs_dupl = IOTools.openFile( outfile + ".duplicates", "w" )
    outs_dupl.write( "contig\tpos\tcounts\n" )

    outs_hist = IOTools.openFile( outfile + ".histogram", "w" )
    outs_hist.write( "duplicates\tcounts\tcumul\tfreq\tcumul_freq\n" )

    last_contig, last_pos = None, None
    ninput, nduplicates = 0, 0

    duplicates = collections.defaultdict( int )
    counts = collections.defaultdict( int )
    count = 0

    # count nh, nm tags
    nh, nm = [], []

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
    '''export macs peaks as bed files.
    '''

    dbhandle = sqlite3.connect( PARAMS["database"] )
    
    if outfile.endswith( ".gz" ):
        compress = True
        track = P.snip( outfile, ".bed.gz" )
    else:
        compress = False
        track = P.snip( outfile, ".bed" )

    tablename = "%s_intervals" % P.quote(track)

    cc = dbhandle.cursor()
    statement = "SELECT contig, start, end, interval_id, peakval FROM %s ORDER by contig, start" % tablename
    cc.execute( statement )

    outs = IOTools.openFile( "%s.bed" % track, "w")

    for result in cc:
        contig, start, end, interval_id,peakval = result
        # peakval is truncated at a 1000 as this is the maximum permitted
        # score in a bed file.
        peakval = int(min(peakval,1000))
        outs.write( "%s\t%i\t%i\t%s\t%i\n" % (contig, start, end, str(interval_id), peakval) )

    cc.close()
    outs.close()

    if compress:
        E.info( "compressing and indexing %s" % outfile )
        use_cluster = True
        statement = 'bgzip -f %(track)s.bed; tabix -f -p bed %(outfile)s'
        P.run()

############################################################
############################################################
############################################################
def exportMacsIntervalsAsBed( infile, outfile, foldchange ):
    '''export sequences for all intervals.'''

    dbhandle = sqlite3.connect( PARAMS["database"] )

    track = P.toTable(os.path.basename( infile ) )
    assert track.endswith("_macs")
    track = track[:-len("_macs")]

    cc = dbhandle.cursor()
    statement = "SELECT contig, start, end, interval_id, fold FROM %(track)s_macs_intervals where fold >= %(foldchange)s ORDER by contig, start" % locals()
    cc.execute( statement )

    outs = open( outfile, "w")

    for result in cc:
        contig, start, end, interval_id,fold = result
        outs.write( "%s\t%i\t%i\t%s\t%d\n" % (contig, start, end, str(interval_id), fold) )
    cc.close()
    outs.close()

############################################################
############################################################
############################################################
def exportPeaksAsBed( infile, outfile ):
    '''export peaks as bed files.'''

    dbhandle = sqlite3.connect( PARAMS["database"] )

    if infile.endswith("_macs.load"):
        track = infile[:-len("_macs.load")]
    else:
        track = infile[:-len("_intervals.load")]
        
    if track.startswith("control"): return
    
    peakwidth = PARAMS["peakwidth"]
    
    cc = dbhandle.cursor()
    statement = '''SELECT contig, peakcenter - %(peakwidth)i, peakcenter + %(peakwidth)i,
                          interval_id, peakval FROM %(track)s_intervals ORDER by contig, start''' % locals()
    cc.execute( statement )

    outs = IOTools.openFile( outfile, "w")

    for result in cc:
        contig, start, end, interval_id,peakval = result
        # peakval is truncated at a 1000 as this is the maximum permitted
        # score in a bed file.
        peakval = int(min(peakval,1000))
        outs.write( "%s\t%i\t%i\t%s\t%i\n" % (contig, start, end, str(interval_id), peakval) )

    cc.close()
    outs.close()

############################################################
############################################################
############################################################
def mergeBedFiles( infiles, outfile ):
    '''generic method for merging bed files. '''

    if len(infiles) < 2:
        raise ValueError( "expected at least two files to merge into %s" % outfile )

    infile = " ".join( infiles )
    statement = '''
        zcat %(infile)s 
        | mergeBed -i stdin 
        | cut -f 1-3 
        | awk '{printf("%%s\\t%%i\\n",$0, ++a); }'
        | bgzip
        > %(outfile)s 
        ''' 

    P.run()

############################################################
############################################################
############################################################
def mergeIntervalsWithScores( infile, outfile, dist, method ):
    '''merge adjacent intervals (within dist) and integrate scores (using method) from merged peaks.
       Assume bed file sorted by position and score in column 5.
       Methods: mean, max, length weighted mean'''

    intervals = open( infile, "r")
    merged = open( outfile, "w")

    topline = intervals.readline()
    last_contig, last_start, last_end, last_id, last_score = topline[:-1].split("\t")[:5]
    last_start = int(last_start)
    last_end = int(last_end)
    last_score = int(last_score)
    for line in intervals:
        data = line[:-1].split("\t")
        contig, start, end, interval_id, score = data[:5]
        start = int(start)
        end = int(end)
        score = int(score)
        if (contig == last_contig) and ((last_end+dist) >= start):
            if method == "mean":
                newscore = (score + last_score) / 2
            elif method == "length_weighted_mean":
                length1 = end - start
                length2 = last_end - last_start
                newscore = ((score*length1)+(last_score*length2))/(length1+length2)
            elif method == "max":
                newscore = max(last_score, score)
            last_end=end
            last_score=newscore
        else:
            merged.write("%(last_contig)s\t%(last_start)i\t%(last_end)i\t%(last_id)s\t%(last_score)s\n" % locals())
            data = line[:-1].split("\t")
            last_contig, last_start, last_end, last_id, last_score = data[:5]
            last_start = int(last_start)
            last_end = int(last_end)
            last_score = int(last_score)
    intervals.close()
    merged.close()

############################################################
############################################################
############################################################
def intersectBedFiles( infiles, outfile ):
    '''merge :term:`bed` formatted *infiles* by intersection
    and write to *outfile*.

    Only intervals that overlap in all files are retained.
    Interval coordinates are given by the first file in *infiles*.

    Bed files are normalized (overlapping intervals within 
    a file are merged) before intersection. 

    Intervals are renumbered starting from 1.
    '''

    if len(infiles) == 1:

        shutil.copyfile( infiles[0], outfile )

    elif len(infiles) == 2:
        
        if P.isEmpty( infiles[0] ) or P.isEmpty( infiles[1] ):
            P.touch( outfile )
        else:
            statement = '''
        intersectBed -u -a %s -b %s 
        | cut -f 1,2,3,4,5 
        | awk 'BEGIN { OFS="\\t"; } {$4=++a; print;}'
        | bgzip > %%(outfile)s 
        ''' % (infiles[0], infiles[1])
            P.run()
        
    else:

        tmpfile = P.getTempFilename(".")

        # need to merge incrementally
        fn = infiles[0]
        if P.isEmpty( infiles[0] ): 
            P.touch( outfile )
            return
            
        statement = '''mergeBed -i %(fn)s > %(tmpfile)s'''
        P.run()
        
        for fn in infiles[1:]:
            if P.isEmpty( infiles[0] ): 
                P.touch( outfile)
                os.unlink( tmpfile )
                return

            statement = '''mergeBed -i %(fn)s | intersectBed -u -a %(tmpfile)s -b stdin > %(tmpfile)s.tmp; mv %(tmpfile)s.tmp %(tmpfile)s'''
            P.run()

        statement = '''cat %(tmpfile)s
        | cut -f 1,2,3,4,5 
        | awk 'BEGIN { OFS="\\t"; } {$4=++a; print;}'
        | bgzip
        > %(outfile)s '''
        P.run()

        os.unlink( tmpfile )

############################################################
############################################################
############################################################
def subtractBedFiles( infile, subtractfile, outfile ):
    '''subtract intervals in *subtractfile* from *infile*
    and store in *outfile*.
    '''

    if P.isEmpty( subtractfile ):
        shutil.copyfile( infile, outfile )
        return
    elif P.isEmpty( infile ):
        P.touch( outfile )
        return

    statement = '''
        intersectBed -v -a %(infile)s -b %(subtractfile)s 
        | cut -f 1,2,3,4,5 
        | awk 'BEGIN { OFS="\\t"; } {$4=++a; print;}'
        | bgzip > %(outfile)s ; tabix -p bed %(outfile)s
        ''' 

    P.run()

############################################################
############################################################
############################################################
def summarizeMACS( infiles, outfile ):
    '''run MACS for peak detection.

    This script parses the MACS logfile to extract 
    peak calling parameters and results.
    '''

    def __get( line, stmt ):
        x = line.search(stmt )
        if x: return x.groups() 

    # mapping patternts to values.
    # tuples of pattern, label, subgroups
    map_targets = [
        ("tags after filtering in treatment: (\d+)", "tag_treatment_filtered",()),
        ("total tags in treatment: (\d+)", "tag_treatment_total",()),
        ("tags after filtering in control: (\d+)", "tag_control_filtered",()),
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

    outs = IOTools.openFile(outfile,"w")

    headers = []
    for k in keys:
        if mapper_header[k]:
            headers.extend( ["%s_%s" % (k,x) for x in mapper_header[k] ])
        else:
            headers.append( k )
    outs.write("track\t%s" % "\t".join(headers) + "\n" )

    for infile in infiles:
        results = collections.defaultdict(list)
        with IOTools.openFile( infile ) as f:
            for line in f:
                if "diag:" in line: break
                for x,y in mapper.items():
                    s = y.search( line )
                    if s: 
                        results[x].append( s.groups()[0] )
                        break
                
        row = [ P.snip( os.path.basename(infile), ".macs" ) ]
        for key in keys:
            val = results[key]
            if len(val) == 0: v = "na"
            else: 
                c = len(mapper_header[key])
                # append missing data (no negative peaks without control files)
                v = "\t".join( map(str, val + ["na"] * (c - len(val)) ))
            row.append(v)
            # assert len(row) -1 == len( headers )
        outs.write("\t".join(row) + "\n" )

    outs.close()

############################################################
############################################################
############################################################
def summarizeMACSsolo( infiles, outfile ):
    '''run MACS for peak detection.'''
    def __get( line, stmt ):
        x = line.search(stmt )
        if x: return x.groups() 

    map_targets = [
        ("total tags in treatment: (\d+)", "tag_treatment_total",()),
        ("#2 number of paired peaks: (\d+)", "paired_peaks",()),
        ("#2   min_tags: (\d+)","min_tags", ()),
        ("#2   d: (\d+)", "shift", ()),
        ("#2   scan_window: (\d+)", "scan_window", ()),
        ("#3 Total number of candidates: (\d+)", "ncandidates",("positive",) ),
        ("#3 Finally, (\d+) peaks are called!",  "called", ("positive",) ) ]

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

        row = [ P.snip( os.path.basename(infile), ".macs" ) ]
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
    
############################################################
############################################################
############################################################
def summarizeMACSFDR( infiles, outfile ):
    '''compile table with peaks that would remain after filtering
    by fdr.
    '''
    
    fdr_thresholds = numpy.arange( 0, 1.05, 0.05 )

    outf = IOTools.openFile( outfile, "w")
    outf.write( "track\t%s\n" % "\t".join( map(str, fdr_thresholds) ) )

    for infile in infiles:
        called = []
        track = P.snip( os.path.basename(infile), ".macs" )
        infilename = infile + "_peaks.xls.gz"
        inf = IOTools.openFile( infilename )
        peaks = list( MACS.iteratePeaks(inf) )
        
        for threshold in fdr_thresholds:
            called.append( len( [ x for x in peaks if x.fdr <= threshold ] ) )
            
        outf.write( "%s\t%s\n" % (track, "\t".join( map(str, called ) ) ) )

    outf.close()

############################################################
############################################################
############################################################
def loadMACS( infile, outfile, bamfile, tablename = None ):
    '''load MACS results in *tablename*

    This method loads only positive peaks. It filters peaks by p-value,
    q-value and fold change and loads the diagnostic data and
    re-calculates peakcenter, peakval, ... using the supplied bamfile.

    If *tablename* is not given, it will be :file:`<track>_intervals`
    where track is derived from ``infile`` and assumed to end
    in :file:`.macs`.

    This method creates two optional additional files:

    * if the file :file:`<track>_diag.xls` is present, load MACS 
    diagnostic data into the table :file:`<track>_macsdiag`.
    
    * if the file :file:`<track>_model.r` is present, call R to
    create a MACS peak-shift plot and save it as :file:`<track>_model.pdf`
    in the :file:`export/MACS` directory.

    This method creates :file:`<outfile>.tsv.gz` with the results
    of the filtering.
    '''

    track = P.snip( os.path.basename(infile), ".macs" )
    folder = os.path.dirname(infile)
    if len(folder) > 0:
        infilename =  folder + "/" + track + "_peaks.xls"
        filename_diag = folder + "/" + track + "_diag.xls"
        filename_r = folder + "/" + track + "_model.r"
        filename_rlog = folder + "/" + track + ".r.log"
        filename_pdf = track + "_model.pdf"
    else: 
        infilename =  track + "_peaks.xls"
        filename_diag = track + "_diag.xls"
        filename_r = track + "_model.r"
        filename_rlog = track + ".r.log"
        filename_pdf = track + "_model.pdf"

    if not os.path.exists(infilename):
        E.warn("could not find %s" % infilename )
        P.touch( outfile )
        return

    # create plot by calling R
    if os.path.exists( filename_r ):
        if len(folder) > 0:
            statement = '''R --vanilla < %(filename_r)s > %(filename_rlog)s; mv %(filename_pdf)s %(folder)s/%(filename_pdf)s; '''
        else:
            statement = '''R --vanilla < %(filename_r)s > %(filename_rlog)s; '''
        P.run()

    # filter peaks
    shift = getPeakShiftFromMacs( infile )
    assert shift != None, "could not determine peak shift from MACS file %s" % infile

    E.info( "%s: found peak shift of %i" % (track, shift ))

    samfiles = [ pysam.Samfile( bamfile, "rb" ) ]
    offsets = [ shift / 2 ]

    outtemp = P.getTempFile()
    tmpfilename = outtemp.name

    outtemp.write( "\t".join( ( \
                "interval_id", 
                "contig", "start", "end",
                "npeaks", "peakcenter", 
                "length", 
                "avgval", 
                "peakval",
                "nprobes",
                "pvalue", "fold", "qvalue",
                "macs_summit", "macs_nprobes",
                )) + "\n" )
    id = 0

    ## get thresholds
    max_qvalue = float(PARAMS["macs_max_qvalue"])
    # min, as it is -10log10
    min_pvalue = float(PARAMS["macs_min_pvalue"])
    
    counter = E.Counter()
    with IOTools.openFile( infilename, "r" ) as ins:
        for peak in MACS.iteratePeaks( ins ):

            if peak.fdr > max_qvalue:
                counter.removed_qvalue += 1
                continue
            elif peak.pvalue < min_pvalue:
                counter.removed_pvalue += 1
                continue

            assert peak.start < peak.end

            npeaks, peakcenter, length, avgval, peakval, nreads = countPeaks( peak.contig, peak.start, peak.end, samfiles, offsets )

            outtemp.write ( "\t".join( map(str, ( \
                            id, peak.contig, peak.start, peak.end, 
                            npeaks, peakcenter, length, avgval, peakval, nreads,
                            peak.pvalue, peak.fold, peak.fdr,
                            peak.start + peak.summit - 1, 
                            peak.tags) ) ) + "\n" )
            id += 1                        
            counter.output += 1

    outtemp.close()

    # output filtering summary
    outf = IOTools.openFile( "%s.tsv.gz" % outfile, "w" )
    outf.write( "category\tcounts\n" )
    outf.write( "%s\n" % counter.asTable() )
    outf.close()

    E.info( "%s filtering: %s" % (track, str(counter)))
    if counter.output == 0:
        E.warn( "%s: no peaks found" % track )

    # load data into table
    if tablename == None:
        tablename = "%s_macs_intervals" % track
    statement = '''python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                       --allow-empty
                       --index=interval_id 
                       --index=contig,start
                       --table=%(tablename)s 
                   < %(tmpfilename)s > %(outfile)s '''
    P.run()
    os.unlink( tmpfilename )

    # load diagnostic data
    if os.path.exists( filename_diag ):

        tablename = "%s_macsdiag" % track
        statement = '''
        cat %(filename_diag)s 
        | sed "s/FC range.*/fc\\tnpeaks\\tp90\\tp80\\tp70\\tp60\\tp50\\tp40\\tp30\\tp20/" 
        | python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                  --map=fc:str 
                  --table=%(tablename)s 
        >> %(outfile)s
        '''
        P.run()        

############################################################
############################################################
############################################################
def runMACS( infile, outfile, controlfile = None ):
    '''run MACS for peak detection from BAM files.

    The output bed files contain the P-value as their score field.
    '''
    to_cluster = True

    if controlfile: control = "--control=%s" % controlfile
    else: control = ""
        
    statement = '''
    macs14 
    -t %(infile)s 
    %(control)s 
    --diag 
    --name=%(outfile)s 
    --format=BAM
    %(macs_options)s 
    >& %(outfile)s
    ''' 
    
    P.run() 
    
    # compress macs bed files and index with tabix
    for suffix in ('peaks', 'summits'):
        statement = '''
        bgzip -f %(outfile)s_%(suffix)s.bed; 
        tabix -f -p bed %(outfile)s_%(suffix)s.bed.gz
        '''
        P.run()
        
    for suffix in ('peaks.xls', 'negative_peaks.xls'):
        statement = '''grep -v "^$" 
                       < %(outfile)s_%(suffix)s 
                       | bgzip > %(outfile)s_%(suffix)s.gz;
                       tabix -f -p bed %(outfile)s_%(suffix)s.gz;
                       checkpoint;
                       rm -f %(outfile)s_%(suffix)s
                    '''
        P.run()

############################################################
############################################################
############################################################
def getCounts( contig, start, end, samfiles, offsets = [] ):
    '''count reads per position.

    If offsets are given, shift tags by offset / 2 and extend
    by offset / 2.
    '''
    assert len(offsets) == 0 or len(samfiles) == len(offsets)

    length = end - start
    counts = numpy.zeros( length )

    nreads = 0

    if offsets:
        # if offsets are given, shift tags. 
        for samfile, offset in zip(samfiles,offsets):

            shift = offset / 2
            # for peak counting I follow the MACS protocoll,
            # see the function def __tags_call_peak in PeakDetect.py
            # In words
            # Only take the start of reads (taking into account the strand)
            # add d/2=offset to each side of peak and start accumulate counts.
            # for counting, extend reads by offset
            # on + strand shift tags upstream
            # i.e. look at the downstream window
            xstart, xend = max(0, start - shift), max(0, end + shift)
            
            for read in samfile.fetch( contig, xstart, xend ):
                pos = read.pos
                if read.is_reverse:
                    # offset = 2 * shift
                    rstart = read.pos + read.alen - offset
                else: 
                    rstart = read.pos + shift

                rend = rstart + shift
                rstart = max( 0, rstart - start )
                rend = min( length, rend - start )
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
def runZinba( infile, outfile, controlfile, action = "full" ):
    '''run Zinba for peak detection.'''

    to_cluster = True

    job_options= "-l mem_free=16G -pe dedicated %i -R y" % PARAMS["zinba_threads"]

    mappability_dir = os.path.join( PARAMS["zinba_mappability_dir"], 
                                    PARAMS["genome"],
                                    "%i" % PARAMS["zinba_read_length"],
                                    "%i" % PARAMS["zinba_alignability_threshold"],
                                    "%i" % PARAMS["zinba_fragment_size"])

    if not os.path.exists( mappability_dir ):
        raise OSError("mappability not found, expected to be at %s" % mappability_dir )

    bit_file = os.path.join( PARAMS["zinba_index_dir"], 
                             PARAMS["genome"] ) + ".2bit"
    if not os.path.exists( bit_file):
        raise OSError("2bit file not found, expected to be at %s" % bit_file )

    options = []
    if controlfile:
        options.append( "--control-filename=%(controlfile)s" % locals() )

    options = " ".join(options)

    statement = '''
    python %(scriptsdir)s/WrapperZinba.py
           --input-format=bam
           --fdr-threshold=%(zinba_fdr_threshold)f
           --fragment-size=%(zinba_fragment_size)s
           --threads=%(zinba_threads)i
           --bit-file=%(bit_file)s
           --mappability-dir=%(mappability_dir)s
           --improvement=%(zinba_improvement)f
           --action=%(action)s
           %(zinba_options)s
           %(options)s
    %(infile)s %(outfile)s
    >& %(outfile)s
    '''

    P.run()


############################################################
############################################################
############################################################
def loadZinba( infile, outfile, bamfile, 
               tablename = None,
               controlfile = None ):
    '''load Zinba results in *tablename*

    This method loads only positive peaks. It filters peaks by p-value,
    q-value and fold change and loads the diagnostic data and
    re-calculates peakcenter, peakval, ... using the supplied bamfile.

    If *tablename* is not given, it will be :file:`<track>_intervals`
    where track is derived from ``infile`` and assumed to end
    in :file:`.zinba`.

    If no peaks were predicted, an empty table is created.

    This method creates :file:`<outfile>.tsv.gz` with the results
    of the filtering.

    This method uses the refined peak locations.

    Zinba peaks can be overlapping. This method does not merge
    overlapping intervals.

    Zinba calls peaks in regions where there are many reads inside
    the control. Thus this method applies a filtering step 
    removing all intervals in which there is a peak of
    more than readlength / 2 height in the control.
    '''

    track = P.snip( os.path.basename(infile), ".zinba" )
    folder = os.path.dirname(infile)

    infilename = infile + ".peaks"

    outtemp = P.getTempFile()
    tmpfilename = outtemp.name

    outtemp.write( "\t".join( ( \
                "interval_id", 
                "contig", "start", "end",
                "npeaks", "peakcenter", 
                "length", 
                "avgval", 
                "peakval",
                "nprobes",
                "pvalue", "fold", "qvalue",
                "macs_summit", "macs_nprobes",
                )) + "\n" )

    counter = E.Counter()
    
    if not os.path.exists(infilename):
        E.warn("could not find %s" % infilename )
    elif P.isEmpty( infile ):
        E.warn("no data in %s" % filename )
    else:
        # filter peaks
        shift = getPeakShiftFromZinba( infile )
        assert shift != None, "could not determine peak shift from Zinba file %s" % infile

        E.info( "%s: found peak shift of %i" % (track, shift ))

        samfiles = [ pysam.Samfile( bamfile, "rb" ) ]
        offsets = [ shift / 2 ]

        if controlfile:
            controlfiles =  [ pysam.Samfile( controlfile, "rb" ) ]
            readlength = PipelineMapping.getReadLengthFromBamfile( controlfile )
            control_max_peakval = readlength // 2
            E.info( "removing intervals in which control has peak higher than %i reads" % control_max_peakval )
        else:
            controlfiles = None

        id = 0

        ## get thresholds
        max_qvalue = float(PARAMS["zinba_fdr_threshold"])

        with IOTools.openFile( infilename, "r" ) as ins:
            for peak in WrapperZinba.iteratePeaks( ins ):

                # filter by qvalue
                if peak.fdr > max_qvalue:
                    counter.removed_qvalue += 1
                    continue

                assert peak.refined_start < peak.refined_end

                # filter by control
                if controlfiles:
                    npeaks, peakcenter, length, avgval, peakval, nreads = countPeaks( peak.contig, 
                                                                                      peak.refined_start, 
                                                                                      peak.refined_end, 
                                                                                      controlfiles, 
                                                                                      offsets )
                    
                    if peakval > control_max_peakval: 
                        counter.removed_control += 1
                        continue

                # output peak
                npeaks, peakcenter, length, avgval, peakval, nreads = countPeaks( peak.contig, 
                                                                                  peak.refined_start, 
                                                                                  peak.refined_end, 
                                                                                  samfiles, 
                                                                                  offsets )

                outtemp.write ( "\t".join( map(str, ( \
                                id, peak.contig, peak.refined_start, peak.refined_end, 
                                npeaks, peakcenter, length, avgval, peakval, nreads,
                                1.0 - peak.posterior, 1.0, peak.fdr,
                                peak.refined_start + peak.summit - 1, 
                                peak.height) ) ) + "\n" )
                id += 1                        
                counter.output += 1

    outtemp.close()

    # output filtering summary
    outf = IOTools.openFile( "%s.tsv.gz" % outfile, "w" )
    outf.write( "category\tcounts\n" )
    outf.write( "%s\n" % counter.asTable() )
    outf.close()

    E.info( "%s filtering: %s" % (track, str(counter)))
    if counter.output == 0:
        E.warn( "%s: no peaks found" % track )

    # load data into table
    if tablename == None:
        tablename = "%s_intervals" % track

    statement = '''
    python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              --allow-empty
              --index=interval_id 
              --index=contig,start
              --table=%(tablename)s 
    < %(tmpfilename)s 
    > %(outfile)s
    '''

    P.run()

    os.unlink( tmpfilename )

############################################################
############################################################
############################################################
##
############################################################
def makeIntervalCorrelation( infiles, outfile, field, reference ):
    '''compute correlation of interval properties between sets
    '''

    dbhandle = sqlite3.connect( PARAMS["database"] )

    tracks, idx = [], []
    for infile in infiles:
        track = P.snip( infile, ".bed.gz" )
        tablename = "%s_intervals" % P.quote( track )
        cc = dbhandle.cursor()
        statement = "SELECT contig, start, end, %(field)s FROM %(tablename)s" % locals()
        cc.execute( statement )
        ix = IndexedGenome.IndexedGenome()
        for contig, start, end, peakval in cc:
            ix.add( contig, start, end, peakval )        
        idx.append( ix )
        tracks.append( track )
    outs = IOTools.openFile( outfile, "w" )
    outs.write( "contig\tstart\tend\tid\t" + "\t".join( tracks ) + "\n" )

    for bed in Bed.iterator( infile = IOTools.openFile( reference, "r") ):
        
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


############################################################
############################################################
############################################################
def buildIntervalCounts( infile, outfile, track, fg_replicates, bg_replicates ):
    '''count read density in bed files comparing stimulated versus unstimulated binding.
    '''
    samfiles_fg, samfiles_bg = [], []

    # collect foreground and background bam files
    for replicate in fg_replicates:
        samfiles_fg.append( "%s.call.bam" % replicate.asFile() )

    for replicate in bg_replicates:
        samfiles_bg.append( "%s.call.bam" % replicate.asFile())
        
    samfiles_fg = [ x for x in samfiles_fg if os.path.exists( x ) ]
    samfiles_bg = [ x for x in samfiles_bg if os.path.exists( x ) ]

    samfiles_fg = ",".join(samfiles_fg)
    samfiles_bg = ",".join(samfiles_bg)

    tmpfile1 = P.getTempFilename( os.getcwd() ) + ".fg"
    tmpfile2 = P.getTempFilename( os.getcwd() ) + ".bg"

    # start counting
    to_cluster = True

    statement = """
    zcat < %(infile)s 
    | python %(scriptsdir)s/bed2gff.py --as-gtf 
    | python %(scriptsdir)s/gtf2table.py 
                --counter=read-coverage 
                --log=%(outfile)s.log 
                --bam-file=%(samfiles_fg)s 
    > %(tmpfile1)s"""
    P.run()

    if samfiles_bg:
        statement = """
        zcat < %(infile)s 
        | python %(scriptsdir)s/bed2gff.py --as-gtf 
        | python %(scriptsdir)s/gtf2table.py 
                    --counter=read-coverage 
                    --log=%(outfile)s.log 
                    --bam-file=%(samfiles_bg)s 
        > %(tmpfile2)s"""
        P.run()

        statement = '''
        python %(toolsdir)s/combine_tables.py 
               --add-file-prefix 
               --regex-filename="[.](\S+)$" 
        %(tmpfile1)s %(tmpfile2)s > %(outfile)s
        '''

        P.run()

        os.unlink( tmpfile2 )
        
    else:
        statement = '''
        python %(toolsdir)s/combine_tables.py 
               --add-file-prefix 
               --regex-filename="[.](\S+)$" 
        %(tmpfile1)s > %(outfile)s
        '''

        P.run()

    os.unlink( tmpfile1 )


def loadIntervalsFromBed( bedfile, track, outfile, 
                          bamfiles, offsets ):
    '''load intervals from :term:`bed` formatted files into database.
    
    Re-evaluate the intervals by counting reads within
    the interval. In contrast to the initial pipeline, the
    genome is not binned. In particular, the meaning of the
    columns in the table changes to:

    nProbes: number of reads in interval
    PeakCenter: position with maximum number of reads in interval
    AvgVal: average coverage within interval

    '''

    tmpfile = P.getTempFile()

    headers = ("AvgVal","DisttoStart","GeneList","Length","PeakCenter","PeakVal","Position","interval_id","nCpGs","nGenes","nPeaks","nProbes","nPromoters", "contig","start","end" )

    tmpfile.write( "\t".join(headers) + "\n" )

    avgval,contig,disttostart,end,genelist,length,peakcenter,peakval,position,start,interval_id,ncpgs,ngenes,npeaks,nprobes,npromoters = \
        0,"",0,0,"",0,0,0,0,0,0,0,0,0,0,0,

    mlength = int(PARAMS["calling_merge_min_interval_length"])

    c = E.Counter()

    # count tags
    for bed in Bed.iterator( IOTools.openFile(infile, "r") ): 

        c.input += 1

        if "name" not in bed:
            bed.name = c.input
        
        # remove very short intervals
        if bed.end - bed.start < mlength: 
            c.skipped_length += 1
            continue

        if replicates:
            npeaks, peakcenter, length, avgval, peakval, nprobes = \
                PipelineChipseq.countPeaks( bed.contig, bed.start, bed.end, samfiles, offsets )

            # nreads can be 0 if the intervals overlap only slightly
            # and due to the binning, no reads are actually in the overlap region.
            # However, most of these intervals should be small and have already be deleted via 
            # the merge_min_interval_length cutoff.
            # do not output intervals without reads.
            if nprobes == 0:
                c.skipped_reads += 1

        else:
            npeaks, peakcenter, length, avgval, peakval, nprobes = ( 1, 
                                                                     bed.start + (bed.end - bed.start) // 2, 
                                                                     bed.end - bed.start, 
                                                                     1, 
                                                                     1,
                                                                     1 )
            
        c.output += 1
        tmpfile.write( "\t".join( map( str, (avgval,disttostart,genelist,length,
                                             peakcenter,peakval,position, bed.name,
                                             ncpgs,ngenes,npeaks,nprobes,npromoters, 
                                             bed.contig,bed.start,bed.end) )) + "\n" )

    if c.output == 0:
        E.warn( "%s - no intervals" )
 
    tmpfile.close()

    tmpfilename = tmpfile.name
    tablename = "%s_intervals" % track.asTable()
    
    statement = '''
    python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --index=interval_id 
              --table=%(tablename)s
    < %(tmpfilename)s 
    > %(outfile)s
    '''

    P.run()
    os.unlink( tmpfile.name )

    L.info( "%s\n" % str(c) )
