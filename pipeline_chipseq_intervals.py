'''ChIP-Seq tasks associated with intervals.
'''

import sys, tempfile, optparse, shutil, itertools, csv, math, random, re, glob, os, shutil, collections
import sqlite3
import cStringIO

import Experiment as E
import Pipeline as P

import csv
import IndexedFasta, IndexedGenome, FastaIterator, Genomics
import IOTools
import MAST, GTF, GFF, Bed
# import Stats

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
def getPeakShift( infile ):
    '''get peak shift for filename infile (.macs output file).

    returns None if no shift found'''

    shift = None
    with open(infile, "r") as ins:
        rx = re.compile("#2 predicted fragment length is (\d+) bps")
        r2 = re.compile("#2 Use (\d)+ as shiftsize, \d+ as fragment length" )
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
def getMappedReads( infile ):
    '''return number of reads mapped. '''
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
def getExonLocations(filename):
    '''return a list of exon locations as (contig,start,end) tuples 
    from a file contain a one ensembl gene ID per line
    '''
    fh = open(filename,"r")
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
        contig, start, end = result
        region_list.append( (contig,int(start),int(end)) )
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

    region_list = getBedLocations(bed_file)
    n_regions = 0
    for region in region_list:
        contig, start, end = region
        #it is neccessary to extend the region to make an accurate mask
        mask.add(contig,(start-1),(end+1),1)
        n_regions += 1

    E.info("Built Quicksect mask for %i regions" % n_regions)

    return(mask)

############################################################
############################################################
############################################################
def buildBAMforPeakCalling( infiles, outfile, dedup, mask):
    ''' Make a BAM file suitable for peak calling.
        Infiles are merged and unmapped reads removed. If specificied
        duplicate reads are removed. If a mask is specified, reads falling within
        the mask are filtered out. The mask is a quicksect object containing
        the regions from which reads are to be excluded.
    '''
    #open the infiles, if more than one merge and sort first using samtools.

    samfiles = []
    num_reads = 0
    nfiles = 0

    if len(infiles) > 1 and isinstance(infiles,str)==0:

        merge = True
        tmpfilename = P.getTempFilename()
        statement = '''
        samtools merge %s %s 
        ''' % (tmpfilename, infiles.join(" "))
        P.run()
        pysam_in = pysam.Samfile(tmpfilename,"rb")
        pysam_in.sort()

    else:
        merge = False
        pysam_in = pysam.Samfile(infiles,"rb")

    pysam_out = pysam.Samfile( outfile, "wb", template = pysam_in )

    #make a simple counter object to pass through the generators
    class read_counter:
        def __init__(self):
            self.n_in, self.n_out, self.unmapped, self.duplicate, self.filtered = 0, 0, 0, 0, 0
    counts = read_counter() 

    ## build the generator
    it = pysam_in

    # discard unmapped reads
    def drop_unmapped(i,c):
        for read in i:
            c.n_in += 1
            if read.is_unmapped:
                c.unmapped += 1 
                continue
            else:
                yield read

    it = drop_unmapped(it,counts)

    # discard duplicate reads if requested
    if dedup:
        def check_for_dup(i,c):

            last_contig, last_start = None, None
            for read in i:
                if read.rname == last_contig and read.pos == last_start:
                    c.duplicate += 1
                    continue
                else:
                    last_contig, last_start = read.rname, read.pos
                    yield read

        it = check_for_dup( it, counts )

    # filter reads if requested
    if mask:
        def filter_reads(i,c,mask):

            for read in i:
                contig = pysam_in.getrname(read.tid)
                if mask.contains(contig,read.pos,read.pos):
                    counts.filtered += 1
                    continue
                else:
                    yield read

        it = filter_reads( it,counts,mask )

    # call the built generator and write out emmitted reads
    for read in it:
        counts.n_out += 1
        pysam_out.write(read)

    pysam_out.close()
    pysam.index( outfile )

    #E.info("Screened %i reads, output %i reads - %5.2f%% from %s" %\
    #      (counts.n_in, counts.n_out, 100.0*counts.n_out/counts.n_in, "+".join(infiles)) )
    #E.info("No. unmapped reads %i, no. duplicates: %i,  no.filtered: %i" %\
    #      (counts.unmapped, counts.duplicate, counts.filtered) ) 

    if merge == True:
        os.unlink( tmpfilename )

############################################################
############################################################
############################################################
def buildSimpleNormalizedBAM( infiles, outfile, nreads ):
    '''normalize a bam file to given number of counts
       by random sampling
    '''
    infile,countfile = infiles

    pysam_in = pysam.Samfile (infile,"rb")
    
    fh = open(countfile,"r")
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
def exportMacsAsBed( infile, outfile ):
    '''export sequences for all intervals.'''

    dbhandle = sqlite3.connect( PARAMS["database"] )

    track = P.toTable( infile )
    assert track.endswith("_macs")
    track = track[:-len("_macs")]

    cc = dbhandle.cursor()
    statement = "SELECT contig, start, end, interval_id, peakval FROM %s_macs_intervals ORDER by contig, start" % track
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
    '''export sequences for all peaks.'''

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
def mergeBedFiles( infiles, outfile ):
    '''generic method for merging bed files. '''

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
    '''generic method for merging bed files.

    Bed files are normalized (overlapping intervals within 
    a file are merged) before intersection. 

    Intervals are renumbered.
    '''

    if len(infiles) == 1:
        shutil.copyfile( infiles[0], outfile )

    elif len(infiles) == 2:
        
        statement = '''
        intersectBed -a %s -b %s 
        | cut -f 1,2,3,4,5 
        | awk 'BEGIN { OFS="\\t"; } {$4=++a; print;}'
        > %%(outfile)s 
        ''' % (infiles[0], infiles[1])

        P.run()
        
    else:

        tmpfile = P.getTempFilename(".")

        # need to merge incrementally
        fn = infiles[0]
        statement = '''mergeBed -i %(fn)s > %(tmpfile)s'''
        P.run()
        
        for fn in infiles[1:]:
            statement = '''mergeBed -i %(fn)s | intersectBed -a %(tmpfile)s -b stdin > %(tmpfile)s.tmp; mv %(tmpfile)s.tmp %(tmpfile)s'''
            P.run()

        statement = '''cat %(tmpfile)s
        | cut -f 1,2,3,4,5 
        | awk 'BEGIN { OFS="\\t"; } {$4=++a; print;}'
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
    statement = '''
        intersectBed -v -a %(infile)s -b %(subtractfile)s |\
        cut -f 1,2,3,4,5 |\
        awk 'BEGIN { OFS="\\t"; } {$4=++a; print;}'
        > %%(outfile)s 
        ''' % locals()

    P.run()

############################################################
############################################################
############################################################
def summarizeMACS( infiles, outfile ):
    '''run MACS for peak detection.'''
    def __get( line, stmt ):
        x = line.search(stmt )
        if x: return x.groups() 

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

###################################################################
###################################################################
###################################################################
###################################################################
def loadMACSSummary( infile, outfile ):
    '''load regions of interest.'''
    
    table = P.snip( os.path.basename(outfile), ".load" )
    statement = '''python %(scriptsdir)s/csv2db.py %(csv2db_options)s
                      --index=track 
                      --table=%(table)s
                   < %(infile)s > %(outfile)s'''
    P.run()

############################################################
############################################################
############################################################
def loadMACS( infile, outfile, bamfile ):
    '''load MACS results.

    Loads only positive peaks. It filters peaks by p-value,
    q-value and fold change and loads the diagnostic data.

    Does re-counting of peakcenter, peakval, ... using
    bamfile.
    '''

    #track = infile[:-len(".macs")]    
    track = P.snip( os.path.basename(infile), ".macs" )
    folder = os.path.dirname(infile)
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

    samfiles = [ pysam.Samfile( bamfile, "rb" ) ]
    offsets = [ shift / 2 ]
    outtemp = open(infile + '.tmp', 'w') #P.getTempFile()

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
            # skip empty lines
            if line.startswith("\n"): continue
            counter.input += 1
            data = line[:-1].split("\t")
            if len(data) == 9:
                contig,start,end,length,summit,ntags,pvalue,fold,qvalue = data
            elif len(data) == 8:
                contig,start,end,length,summit,ntags,pvalue,fold = data
                qvalue = 0.0
            else:
                raise ValueError( "could not parse line %s" % line )
            
            # qvalue is in percent, divide by 100.
            pvalue, qvalue, summit, fold = float(pvalue), float(qvalue) / 100, int(summit), float(fold)

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

    tablename = "%s_macs_intervals" % track
    tmpfilename = outtemp.name

    statement = '''python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                       --index=interval_id 
                       --index=contig,start
                       --table=%(tablename)s 
                   < %(tmpfilename)s > %(outfile)s'''
    P.run()

    # load diagnostic data
    if os.path.exists( filename_diag ):

        tablename = "%s_macsdiag" % track

        statement = '''
        sed "s/FC range.*/fc\\tnpeaks\\tp90\\tp80\\tp70\\tp60\\tp50\\tp40\\tp30\\tp20/" < %(filename_diag)s |\
       python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
                  --map=fc:str \
                  --table=%(tablename)s \
        > %(outfile)s
        '''

        P.run()

    # create plot
    if os.path.exists( filename_r ):

        target_path = os.path.join( os.getcwd(), "export", "MACS" )
        try:
            os.makedirs( target_path )
        except OSError: 
            # ignore "file exists" exception
            pass

        statement = '''R --vanilla < %(folder)s/%(track)s.macs_model.r > %(outfile)s '''
        
        P.run()

        #shutil.copyfile("%s.macs_model.pdf" % track, os.path.join( target_path, "%s_model.pdf" % track) )
        
    os.unlink( tmpfilename )

    E.info("%s: %s" % (track, str(counter))) 

############################################################
############################################################
############################################################
def runMACS( infile, outfile, controlfile = None ):
    '''run MACS for peak detection from BAM files.

    The output bed files contain the P-value as their score field.
    '''
    to_cluster = True

    if controlfile: control = "-c %s" % controlfile
    else: control = ""
        
    statement = '''
    macs -t %(infile)s %(control)s \
    --diag \
    --name=%(outfile)s \
    --format=%(format)s \
    %(macs_options)s >& %(outfile)s''' 
    
    P.run() 

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
def makeIntervalCorrelation( infiles, outfile, field, reference ):
    '''compute correlation of interval properties between sets
    '''

    dbhandle = sqlite3.connect( PARAMS["database"] )

    tracks, idx = [], []
    for infile in infiles:
        track = P.snip( infile, ".bed" )
        tablename = "%s_intervals" % P.quote( track )
        cc = dbhandle.cursor()
        statement = "SELECT contig, start, end, %(field)s FROM %(tablename)s" % locals()
        cc.execute( statement )
        ix = IndexedGenome.IndexedGenome()
        for contig, start, end, peakval in cc:
            ix.add( contig, start, end, peakval )        
        idx.append( ix )
        tracks.append( track )
    outs = open( outfile, "w" )
    outs.write( "contig\tstart\tend\tid\t" + "\t".join( tracks ) + "\n" )

    for bed in Bed.iterator( infile = open( reference, "r") ):
        
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

