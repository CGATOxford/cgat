'''PipelinePeakCalling.py - ChIP-Seq tasks associated with calling peaks
========================================================================

def runXXX():
    run a certain peak caller

def loadXXX():
    load results from peak caller. The minimum information required are:
        contig, start, end
'''

import shutil
import random
import re
import glob
import os
import collections
import sqlite3
import numpy
import pysam

# pybedtools recompilation can fail causing
# an import error when importing this script
try:
    import pybedtools
except ImportError:
    pass

##########################
import CGAT.Experiment as E
import CGAT.Pipeline as P
import CGAT.IndexedGenome as IndexedGenome
import CGAT.IOTools as IOTools
import CGAT.BamTools as BamTools
import CGAT.Bed as Bed
import CGAT.WrapperMACS as WrapperMACS

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################
P.getParameters(
    ["%s.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini"])

PARAMS = P.PARAMS

############################################################
############################################################
############################################################
##
############################################################


def getPeakShiftFromMacs(infile):
    '''get peak shift for filename infile (.macs output file).

    returns None if no shift found'''

    shift = None
    with IOTools.openFile(infile, "r") as ins:
        rx = re.compile("#2 predicted fragment length is (\d+) bps")
        r2 = re.compile("#2 Use (\d+) as shiftsize, \d+ as fragment length")
        r3 = re.compile("#1 fragment size = (\d+)")
        # when fragment length is set explicitely
        r4 = re.compile("#2 Use (\d+) as fragment length")

        for line in ins:
            x = rx.search(line)
            if x:
                shift = int(x.groups()[0])
                break
            x = r3.search(line)
            if x:
                shift = int(x.groups()[0])
                break
            x = r4.search(line)
            if x:
                shift = int(x.groups()[0])
                break
            x = r2.search(line)
            if x:
                shift = int(x.groups()[0])
                E.warn("shift size was set automatically - see MACS logfiles")
                break

    return shift

############################################################
############################################################
############################################################
##
############################################################


def getPeakShiftFromZinba(infile):
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
                shift = int(lines[i + 1].split()[1])
                break

    return shift

############################################################
############################################################
############################################################
##
############################################################


def getPeakShiftFromSPP(infile):
    '''get peak shift for filename infile (.spp output file).

    returns None if no shift found
    '''

    shift = None

    # search for
    # shift\t125
    with IOTools.openFile(infile, "r") as ins:
        lines = ins.readlines()
        for i, line in enumerate(lines):
            if line.startswith("shift\t"):
                shift = int(re.match("shift\t(\d+)\n", line).groups()[0])
                break

    return shift

############################################################
############################################################
############################################################


def getPeakShift(track):
    '''get peak shift for a track or filename.'''

    if track.endswith(".macs"):
        return getPeakShiftFromMacs(track)
    elif track.endswith(".zinba"):
        return getPeakShiftFromZinba(track)
    elif track.endswith(".spp"):
        return getPeakShiftFromSPP(track)
    elif os.path.exists("%s.macs" % track):
        return getPeakShiftFromMacs("%s.macs" % track)
    elif os.path.exists("%s.zinba" % track):
        return getPeakShiftFromZinba("%s.zinba" % track)
    elif os.path.exists("%s.spp" % track):
        return getPeakShiftFromSPP("%s.spp" % track)


############################################################
############################################################
############################################################
def getMappedReads(infile):
    '''return number of reads mapped. '''
    for lines in IOTools.openFile(infile, "r"):
        data = lines[:-1].split("\t")
        if data[1].startswith("without duplicates"):
            return int(data[0])
    return

############################################################
############################################################
############################################################


def getMinimumMappedReads(infiles):
    '''find the minimum number of mapped reads in infiles.'''
    v = []
    for infile in infiles:
        x = getMappedReads(infile)
        if x:
            v.append(x)
    if len(v) == 0:
        raise P.PipelineError(
            "could not find mapped reads in files %s" % (str(infiles)))
    return min(v)

############################################################
############################################################
############################################################


def getExonLocations(filename):
    '''return a list of exon locations as Bed entries
    from a file contain a one ensembl gene ID per line
    '''
    fh = IOTools.openFile(filename, "r")
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

    statement = "select contig,start,end from geneset_cds_gtf where " + \
        " OR ".join(gene_ids)

    cc.execute(statement)

    region_list = []
    n_regions = 0
    for result in cc:
        b = Bed.Bed()
        b.contig, b.start, b.end = result
        region_list.append(b)
        n_regions += 1

    cc.close()

    E.info("Retrieved exon locations for %i genes. Got %i regions" %
           (n_ids, n_regions))

    return(region_list)

############################################################
############################################################
############################################################


def getBedLocations(filename):
    '''return a list of regions as (contig,start,end) tuples
    from a bed file'''
    fh = open(filename, "r")
    region_list = []
    n_regions = 0

    for line in fh:
        if line.strip() != "" and line[0] != "#":
            fields = line.split("\t")
            contig, start, end = fields[0], int(fields[1]), int(fields[2])
            region_list.append((contig, start, end))
            n_regions += 1

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
    for bed in Bed.iterator(IOTools.openFile(bed_file)):
        # it is neccessary to extend the region to make an accurate mask
        mask.add(bed.contig, (bed.start - 1), (bed.end + 1), 1)
        n_regions += 1

    E.info("Built Quicksect mask for %i regions" % n_regions)

    return(mask)

############################################################
############################################################
############################################################


def buildBAMforPeakCalling(infiles, outfile, dedup, mask):
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

    # open the infiles, if more than one merge and sort first using samtools.

    samfiles = []
    num_reads = 0
    nfiles = 0

    to_cluster = True

    statement = []

    tmpfile = P.getTempFilename()

    if len(infiles) > 1 and isinstance(infiles, str) == 0:
        # assume: samtools merge output is sorted
        # assume: sam files are sorted already
        statement.append( '''samtools merge @OUT@ %s''' % (infiles.join(" ")) )
        statement.append( '''samtools sort @IN@ @OUT@''')

    if dedup:
        job_options = "-l mem_free=16G -l picard=1"
        statement.append( '''MarkDuplicates
                                       INPUT=@IN@
                                       ASSUME_SORTED=true 
                                       REMOVE_DUPLICATES=true
                                       QUIET=true
                                       OUTPUT=@OUT@
                                       METRICS_FILE=%(outfile)s.picardmetrics
                                       VALIDATION_STRINGENCY=SILENT 
                   >& %(outfile)s.picardlog ''' )

    if mask:
        statement.append(
            '''intersectBed -abam @IN@ -b %(mask)s -wa -v > @OUT@''' )

    statement.append('''mv @IN@ %(outfile)s''' )
    statement.append('''samtools index %(outfile)s''' )

    statement = P.joinStatements(statement, infiles)
    P.run()

############################################################
############################################################
############################################################


def buildSimpleNormalizedBAM(infiles, outfile, nreads):
    '''normalize a bam file to given number of counts
       by random sampling
    '''
    infile, countfile = infiles

    pysam_in = pysam.Samfile(infile, "rb")

    fh = IOTools.openFile(countfile, "r")
    readcount = int(fh.read())
    fh.close()

    threshold = float(nreads) / float(readcount)

    pysam_out = pysam.Samfile(outfile, "wb", template=pysam_in)

    # iterate over mapped reads thinning by the threshold
    ninput, noutput = 0, 0
    for read in pysam_in.fetch():
        ninput += 1
        if random.random() <= threshold:
            pysam_out.write(read)
            noutput += 1

    pysam_in.close()
    pysam_out.close()
    pysam.index(outfile)

    E.info("buildNormalizedBam: %i input, %i output (%5.2f%%), should be %i" %
           (ninput, noutput, 100.0 * noutput / ninput, nreads))

############################################################
############################################################
############################################################
####### Depreciate this function? ##########################
############################################################


def buildNormalizedBAM(infiles, outfile, normalize=True):
    '''build a normalized BAM file.

    Infiles are merged and duplicated reads are removed. 
    If *normalize* is set, reads are removed such that all 
    files will have approximately the same number of reads.

    Note that the duplication here is wrong as there
    is no sense of strandedness preserved.
    '''

    min_reads = getMinimumMappedReads(glob.glob("*.readstats"))

    samfiles = []
    num_reads = 0
    for infile, statsfile in infiles:
        samfiles.append(pysam.Samfile(infile, "rb"))
        num_reads += getMappedReads(statsfile)

    threshold = float(min_reads) / num_reads

    E.info("%s: min reads: %i, total reads=%i, threshold=%f" %
           (infiles, min_reads, num_reads, threshold))

    pysam_out = pysam.Samfile(outfile, "wb", template=samfiles[0])

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
                pysam_out.write(read)
                noutput += 1

            last_contig, last_pos = read.rname, read.pos

        pysam_in.close()

    pysam_out.close()

    logs = IOTools.openFile(outfile + ".log", "w")
    logs.write("# min_reads=%i, threshold= %5.2f\n" %
               (min_reads, threshold))
    logs.write("set\tcounts\tpercent\n")
    logs.write("ninput\t%i\t%5.2f%%\n" % (ninput, 100.0))
    nwithout_dups = ninput - nduplicates
    logs.write("duplicates\t%i\t%5.2f%%\n" %
               (nduplicates, 100.0 * nduplicates / ninput))
    logs.write("without duplicates\t%i\t%5.2f%%\n" %
               (nwithout_dups, 100.0 * nwithout_dups / ninput))
    logs.write("target\t%i\t%5.2f%%\n" %
               (min_reads, 100.0 * min_reads / nwithout_dups))
    logs.write("noutput\t%i\t%5.2f%%\n" %
               (noutput, 100.0 * noutput / nwithout_dups))

    logs.close()

    # if more than one samfile: sort
    if len(samfiles) > 1:
        tmpfilename = P.getTempFilename()
        pysam.sort(outfile, tmpfilename)
        shutil.move(tmpfilename + ".bam", outfile)
        os.unlink(tmpfilename)

    pysam.index(outfile)

    E.info("buildNormalizedBam: %i input, %i output (%5.2f%%), should be %i" %
           (ninput, noutput, 100.0 * noutput / ninput, min_reads))


############################################################
############################################################
############################################################
def buildBAMStats(infile, outfile):
    '''calculate bamfile statistics - currently only single-ended
    duplicates.
    '''

    # no bedToBigBed
    # to_cluster = True
    outs = IOTools.openFile(outfile, "w")
    outs.write("reads\tcategory\n")
    for line in pysam.flagstat(infile):
        data = line[:-1].split(" ")
        outs.write("%s\t%s\n" % (data[0], " ".join(data[1:])))

    pysam_in = pysam.Samfile(infile, "rb")

    outs_dupl = IOTools.openFile(outfile + ".duplicates", "w")
    outs_dupl.write("contig\tpos\tcounts\n")

    outs_hist = IOTools.openFile(outfile + ".histogram", "w")
    outs_hist.write("duplicates\tcounts\tcumul\tfreq\tcumul_freq\n")

    last_contig, last_pos = None, None
    ninput, nduplicates = 0, 0

    duplicates = collections.defaultdict(int)
    counts = collections.defaultdict(int)
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
            outs_dupl.write("%s\t%i\t%i\n" % (last_contig, last_pos, count))
            counts[count] += 1

        count = 1
        last_contig, last_pos = read.rname, read.pos

    outs.write("%i\tduplicates (%5.2f%%)\n" %
               (nduplicates, 100.0 * nduplicates / ninput))
    outs.write("%i\twithout duplicates (%5.2f%%)\n" % (ninput - nduplicates,
                                                       100.0 * (ninput - nduplicates) / ninput))
    pysam_in.close()
    outs.close()
    outs_dupl.close()

    keys = counts.keys()
    # count per position (not the same as nduplicates, which is # of reads)
    c = 0
    total = sum(counts.values())
    for k in sorted(keys):
        c += counts[k]
        outs_hist.write("%i\t%i\t%i\t%f\t%f\n" % (k, counts[k], c,
                                                  100.0 * counts[k] / total,
                                                  100.0 * c / total))
    outs_hist.close()


############################################################
############################################################
############################################################
def exportIntervalsAsBed(infile, outfile,
                         tablename,
                         bedfilter=None,
                         merge=False):
    '''export intervals from database as bed files.

    If *bedfilter* is set, remove intervals overlapping
    any of the intervals in the :term:`bed` formatted file.

    If *merge* is set, merge overlapping intervals.

    Returns a counter object.
    '''

    dbhandle = sqlite3.connect(PARAMS["database"])

    if outfile.endswith(".gz"):
        compress = True
        track = P.snip(outfile, ".bed.gz")
    else:
        compress = False
        track = P.snip(outfile, ".bed")

    cc = dbhandle.cursor()
    statement = """SELECT contig, max(0, start), end
    FROM %s""" % tablename
    cc.execute(statement)

    bd = pybedtools.BedTool(list(cc.fetchall()))
    tmpfile = P.getTempFilename()
    bd.saveas(tmpfile)
    cc.close()

    c = E.Counter()
    bd = pybedtools.BedTool(tmpfile)
    c.input = len(bd)
    bd = pybedtools.BedTool(tmpfile)
    latest = c.input

    if bedfilter:

        bd = bd.intersect(pybedtools.BedTool(bedfilter), v=True, wa=True)
        bd.saveas(tmpfile)
        bd = pybedtools.BedTool(tmpfile)

        c.after_bedfilter = len(bd)
        bd = pybedtools.BedTool(tmpfile)
        
        c.removed_bedfilter = latest - c.after_bedfilter
        latest = c.after_bedfilter

    if merge and latest > 0:
        # empty bedfiles cause an error
        # pybedtools not very intuitive.
        bd = bd.sort()
        bd.saveas(tmpfile)
        bd = pybedtools.BedTool(tmpfile)

        bd = bd.merge()
        bd.saveas(tmpfile)
        bd = pybedtools.BedTool(tmpfile)

        c.after_merging = len(bd)
        bd = pybedtools.BedTool(tmpfile)

        c.removed_merging = latest - c.after_merging
        latest = c.after_merging

    c.output = latest        

    # one final sort
    bd = bd.sort()
    bd.saveas(tmpfile)
    bd = pybedtools.BedTool(tmpfile)
    bd.saveas(track + '.bed')

    if compress:
        E.info("compressing and indexing %s" % outfile)
        statement = 'bgzip -f %(track)s.bed; tabix -f -p bed %(outfile)s'
        P.run()

    os.unlink(tmpfile)
    return c

############################################################
############################################################
############################################################


def exportMacsIntervalsAsBed(infile, outfile, foldchange):
    '''export sequences for all intervals.'''

    dbhandle = sqlite3.connect(PARAMS["database"])

    track = P.toTable(os.path.basename(infile))
    assert track.endswith("_macs")
    track = track[:-len("_macs")]

    cc = dbhandle.cursor()
    statement = """SELECT contig, start, end, interval_id, fold
    FROM %(track)s_macs_intervals where fold >= %(foldchange)s
    ORDER by contig, start""" % locals()
    cc.execute(statement)

    outs = open(outfile, "w")

    for result in cc:
        contig, start, end, interval_id, fold = result
        outs.write("%s\t%i\t%i\t%s\t%d\n" %
                   (contig, start, end, str(interval_id), fold))
    cc.close()
    outs.close()

############################################################
############################################################
############################################################


def exportPeaksAsBed(infile, outfile):
    '''export peaks as bed files.'''

    dbhandle = sqlite3.connect(PARAMS["database"])

    if infile.endswith("_macs.load"):
        track = infile[:-len("_macs.load")]
    else:
        track = infile[:-len("_intervals.load")]

    if track.startswith("control"):
        return

    peakwidth = PARAMS["peakwidth"]

    cc = dbhandle.cursor()
    statement = '''SELECT contig, peakcenter - %(peakwidth)i, peakcenter + %(peakwidth)i,
                          interval_id, peakval FROM %(track)s_intervals ORDER by contig, start''' % locals()
    cc.execute(statement)

    outs = IOTools.openFile(outfile, "w")

    for result in cc:
        contig, start, end, interval_id, peakval = result
        # peakval is truncated at a 1000 as this is the maximum permitted
        # score in a bed file.
        peakval = int(min(peakval, 1000))
        outs.write("%s\t%i\t%i\t%s\t%i\n" %
                   (contig, start, end, str(interval_id), peakval))

    cc.close()
    outs.close()

############################################################
############################################################
############################################################


def mergeBedFiles(infiles, outfile):
    '''generic method for merging bed files. '''

    if len(infiles) < 2:
        raise ValueError(
            "expected at least two files to merge into %s" % outfile)

    infile = " ".join(infiles)
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


def mergeIntervalsWithScores(infile, outfile, dist, method):
    '''merge adjacent intervals (within dist) and integrate scores (using method) from merged peaks.
       Assume bed file sorted by position and score in column 5.
       Methods: mean, max, length weighted mean'''

    intervals = open(infile, "r")
    merged = open(outfile, "w")

    topline = intervals.readline()
    last_contig, last_start, last_end, last_id, last_score = topline[
        :-1].split("\t")[:5]
    last_start = int(last_start)
    last_end = int(last_end)
    last_score = int(last_score)
    for line in intervals:
        data = line[:-1].split("\t")
        contig, start, end, interval_id, score = data[:5]
        start = int(start)
        end = int(end)
        score = int(score)
        if (contig == last_contig) and ((last_end + dist) >= start):
            if method == "mean":
                newscore = (score + last_score) / 2
            elif method == "length_weighted_mean":
                length1 = end - start
                length2 = last_end - last_start
                newscore = (
                    (score * length1) + (last_score * length2)) / (length1 + length2)
            elif method == "max":
                newscore = max(last_score, score)
            last_end = end
            last_score = newscore
        else:
            merged.write(
                "%(last_contig)s\t%(last_start)i\t%(last_end)i\t%(last_id)s\t%(last_score)s\n" % locals())
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


def intersectBedFiles(infiles, outfile):
    '''merge :term:`bed` formatted *infiles* by intersection
    and write to *outfile*.

    Only intervals that overlap in all files are retained.
    Interval coordinates are given by the first file in *infiles*.

    Bed files are normalized (overlapping intervals within 
    a file are merged) before intersection. 

    Intervals are renumbered starting from 1.
    '''

    if len(infiles) == 1:

        shutil.copyfile(infiles[0], outfile)

    elif len(infiles) == 2:

        if P.isEmpty(infiles[0]) or P.isEmpty(infiles[1]):
            P.touch(outfile)
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
        if P.isEmpty(infiles[0]):
            P.touch(outfile)
            return

        statement = '''mergeBed -i %(fn)s > %(tmpfile)s'''
        P.run()

        for fn in infiles[1:]:
            if P.isEmpty(infiles[0]):
                P.touch(outfile)
                os.unlink(tmpfile)
                return

            statement = '''mergeBed -i %(fn)s | intersectBed -u -a %(tmpfile)s -b stdin > %(tmpfile)s.tmp; mv %(tmpfile)s.tmp %(tmpfile)s'''
            P.run()

        statement = '''cat %(tmpfile)s
        | cut -f 1,2,3,4,5 
        | awk 'BEGIN { OFS="\\t"; } {$4=++a; print;}'
        | bgzip
        > %(outfile)s '''
        P.run()

        os.unlink(tmpfile)

############################################################
############################################################
############################################################


def subtractBedFiles(infile, subtractfile, outfile):
    '''subtract intervals in *subtractfile* from *infile*
    and store in *outfile*.
    '''

    if P.isEmpty(subtractfile):
        shutil.copyfile(infile, outfile)
        return
    elif P.isEmpty(infile):
        P.touch(outfile)
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


def summarizeMACS(infiles, outfile):
    '''run MACS for peak detection.

    This script parses the MACS logfile to extract 
    peak calling parameters and results.
    '''

    def __get(line, stmt):
        x = line.search(stmt)
        if x:
            return x.groups()

    # mapping patternts to values.
    # tuples of pattern, label, subgroups
    map_targets = [
        ("tags after filtering in treatment: (\d+)",
         "tag_treatment_filtered", ()),
        ("total tags in treatment: (\d+)", "tag_treatment_total", ()),
        ("tags after filtering in control: (\d+)", "tag_control_filtered", ()),
        ("total tags in control: (\d+)", "tag_control_total", ()),
        ("#2 number of paired peaks: (\d+)", "paired_peaks", ()),
        ("#2   min_tags: (\d+)", "min_tags", ()),
        ("#2   d: (\d+)", "shift", ()),
        ("#2   scan_window: (\d+)", "scan_window", ()),
        ("#3 Total number of candidates: (\d+)",
         "ncandidates", ("positive", "negative")),
        ("#3 Finally, (\d+) peaks are called!",  "called", ("positive", "negative"))]

    mapper, mapper_header = {}, {}
    for x, y, z in map_targets:
        mapper[y] = re.compile(x)
        mapper_header[y] = z

    keys = [x[1] for x in map_targets]

    outs = IOTools.openFile(outfile, "w")

    headers = []
    for k in keys:
        if mapper_header[k]:
            headers.extend(["%s_%s" % (k, x) for x in mapper_header[k]])
        else:
            headers.append(k)
    outs.write("track\t%s" % "\t".join(headers) + "\n")

    for infile in infiles:
        results = collections.defaultdict(list)
        with IOTools.openFile(infile) as f:
            for line in f:
                if "diag:" in line:
                    break
                for x, y in mapper.items():
                    s = y.search(line)
                    if s:
                        results[x].append(s.groups()[0])
                        break

        row = [P.snip(os.path.basename(infile), ".macs")]
        for key in keys:
            val = results[key]
            if len(val) == 0:
                v = "na"
            else:
                c = len(mapper_header[key])
                # append missing data (no negative peaks without control files)
                v = "\t".join(map(str, val + ["na"] * (c - len(val))))
            row.append(v)
            # assert len(row) -1 == len( headers )
        outs.write("\t".join(row) + "\n")

    outs.close()

############################################################
############################################################
############################################################


def summarizeMACSsolo(infiles, outfile):
    '''run MACS for peak detection.'''
    def __get(line, stmt):
        x = line.search(stmt)
        if x:
            return x.groups()

    map_targets = [
        ("total tags in treatment: (\d+)", "tag_treatment_total", ()),
        ("#2 number of paired peaks: (\d+)", "paired_peaks", ()),
        ("#2   min_tags: (\d+)", "min_tags", ()),
        ("#2   d: (\d+)", "shift", ()),
        ("#2   scan_window: (\d+)", "scan_window", ()),
        ("#3 Total number of candidates: (\d+)", "ncandidates", ("positive",)),
        ("#3 Finally, (\d+) peaks are called!",  "called", ("positive",))]

    mapper, mapper_header = {}, {}
    for x, y, z in map_targets:
        mapper[y] = re.compile(x)
        mapper_header[y] = z

    keys = [x[1] for x in map_targets]

    outs = open(outfile, "w")

    headers = []
    for k in keys:
        if mapper_header[k]:
            headers.extend(["%s_%s" % (k, x) for x in mapper_header[k]])
        else:
            headers.append(k)
    outs.write("track\t%s" % "\t".join(headers) + "\n")

    for infile in infiles:
        results = collections.defaultdict(list)
        with open(infile) as f:
            for line in f:
                if "diag:" in line:
                    break
                for x, y in mapper.items():
                    s = y.search(line)
                    if s:
                        results[x].append(s.groups()[0])
                        break

        row = [P.snip(os.path.basename(infile), ".macs")]
        for key in keys:
            val = results[key]
            if len(val) == 0:
                v = "na"
            else:
                c = len(mapper_header[key])
                if c >= 1:
                    assert len(val) == c, "key=%s, expected=%i, got=%i, val=%s, c=%s" %\
                        (key,
                         len(val),
                            c,
                            str(val), mapper_header[key])
                v = "\t".join(val)
            row.append(v)
        outs.write("\t".join(row) + "\n")

    outs.close()

############################################################
############################################################
############################################################


def summarizeMACSFDR(infiles, outfile):
    '''compile table with peaks that would remain after filtering
    by fdr.
    '''

    fdr_thresholds = numpy.arange(0, 1.05, 0.05)

    outf = IOTools.openFile(outfile, "w")
    outf.write("track\t%s\n" % "\t".join(map(str, fdr_thresholds)))

    for infile in infiles:
        called = []
        track = P.snip(os.path.basename(infile), ".macs")
        infilename = infile + "_peaks.xls.gz"
        inf = IOTools.openFile(infilename)
        peaks = list(WrapperMACS.iterateMacsPeaks(inf))

        for threshold in fdr_thresholds:
            called.append(len([x for x in peaks if x.fdr <= threshold]))

        outf.write("%s\t%s\n" % (track, "\t".join(map(str, called))))

    outf.close()

############################################################
############################################################
############################################################


def runMACS(infile, outfile,
            controlfile=None,
            tagsize=None):
    '''run MACS for peak detection from BAM files.

    The output bed files contain the P-value as their score field.
    Output bed files are compressed and indexed.
    '''
    job_options = "-l mem_free=8G"

    options = []
    if controlfile:
        options.append("--control=%s" % controlfile)
    if tagsize is not None:
        options.append("--tsize %i" % tagsize)
    
    options = " ".join(options)

    statement = '''
    macs14
    -t %(infile)s
    --diag
    --verbose=10
    --name=%(outfile)s
    --format=BAM
    %(options)s
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


def summarizeMACS2(infiles, outfile):
    '''run MACS2 for peak detection.

    This script parses the MACS2 logfile to extract
    peak calling parameters and results.

    TODO: doesn't report peak numbers...
    '''

    def __get(line, stmt):
        x = line.search(stmt)
        if x:
            return x.groups()

    # mapping patternts to values.
    # tuples of pattern, label, subgroups
    map_targets = [
        ("tags after filtering in treatment:\s+(\d+)",
         "fragment_treatment_filtered", ()),
        ("total tags in treatment:\s+(\d+)",
         "fragment_treatment_total", ()),
        ("tags after filtering in control:\s+(\d+)",
         "fragment_control_filtered", ()),
        ("total tags in control:\s+(\d+)",
         "fragment_control_total", ()),
        ("predicted fragment length is (\d+) bps",
         "fragment_length", ()),
        # Number of peaks doesn't appear to be reported!.
        ("#3 Total number of candidates: (\d+)",
         "ncandidates", ("positive", "negative")),
        ("#3 Finally, (\d+) peaks are called!",
         "called",
         ("positive", "negative"))
    ]

    mapper, mapper_header = {}, {}
    for x, y, z in map_targets:
        mapper[y] = re.compile(x)
        mapper_header[y] = z

    keys = [x[1] for x in map_targets]

    outs = IOTools.openFile(outfile, "w")

    headers = []
    for k in keys:
        if mapper_header[k]:
            headers.extend(["%s_%s" % (k, x) for x in mapper_header[k]])
        else:
            headers.append(k)
    outs.write("track\t%s" % "\t".join(headers) + "\n")

    for infile in infiles:
        results = collections.defaultdict(list)
        with IOTools.openFile(infile) as f:
            for line in f:
                if "diag:" in line:
                    break
                for x, y in mapper.items():
                    s = y.search(line)
                    if s:
                        results[x].append(s.groups()[0])
                        break

        row = [P.snip(os.path.basename(infile), ".macs2")]
        for key in keys:
            val = results[key]
            if len(val) == 0:
                v = "na"
            else:
                c = len(mapper_header[key])
                # append missing data (no negative peaks without control files)
                v = "\t".join(map(str, val + ["na"] * (c - len(val))))
            row.append(v)
            # assert len(row) -1 == len( headers )
        outs.write("\t".join(row) + "\n")

    outs.close()

############################################################
############################################################
############################################################


def summarizeMACS2FDR(infiles, outfile):
    '''compile table with peaks that would remain after filtering
    by fdr.
    '''
    fdr_threshold = PARAMS["macs2_max_qvalue"]  # numpy.arange( 0, 1.05, 0.05 )

    outf = IOTools.openFile(outfile, "w")
    outf.write("track\t%s\n" % str(fdr_threshold))

    for infile in infiles:
        called = []
        track = P.snip(os.path.basename(infile), ".macs2")
        infilename = infile + "_peaks.xls.gz"
        inf = IOTools.openFile(infilename)
        peaks = list(WrapperMACS.iterateMacs2Peaks(inf))

        # for threshold in fdr_thresholds:
        called.append(len([x for x in peaks if x.qvalue <= fdr_threshold]))

        outf.write("%s\t%s\n" % (track, "\t".join(map(str, called))))

    outf.close()


############################################################
############################################################
############################################################
def bedGraphToBigwig(infile, contigsfile, outfile, remove=True):
    '''convert a bedgraph file to a bigwig file.

    The bedgraph file is deleted on success.
    '''

    if not os.path.exists(infile):
        raise OSError("bedgraph file %s does not exist" % infile)

    if not os.path.exists(contigsfile):
        raise OSError("contig size file %s does not exist" % infile)

    statement = '''
         bedGraphToBigWig %(infile)s %(contigsfile)s %(outfile)s
    '''
    P.run()

    if os.path.exists(outfile) and not P.isEmpty(outfile):
        os.remove(infile)

############################################################
############################################################
############################################################


def runMACS2(infile, outfile,
             controlfile=None,
             force_single_end=False,
             tagsize=None):
    '''run MACS for peak detection from BAM files.

    The output bed files contain the P-value as their score field.
    Output bed files are compressed and indexed.

    If *force_single_end* is set, do not use paired-ended
    mode for paired ended BAM files.

    Build bedgraph files and convert to bigwig files.
    '''
    options = []

    if controlfile:
        options.append("--control %s" % controlfile)
    if tagsize is not None:
        options.append("--tsize %i" % tagsize)
    
    options = " ".join(options)

    job_options = "-l mem_free=8G"

    # example statement: macs2 callpeak -t R1-paupar-R1.call.bam -c
    # R1-lacZ-R1.call.bam -f BAMPE -g 2.39e9 --verbose 5 --bw 150 -q
    # 0.01 -m 10 100000 --name test

    # used to set the option --format=bampe
    # removed to let macs2 detect the format.

    # format bam needs to be set explicitely, autodetection does not
    # work. Use paired end mode to detect tag size
    if not force_single_end and BamTools.isPaired(infile):
        format_options = '--format=BAMPE'
    else:
        format_options = '--format=BAM'

    # --bdg --SPMR: ask macs to create a bed-graph file with
    # fragment pileup per million reads
    statement = '''
    macs2 callpeak
    %(format_options)s
    --treatment %(infile)s
    --verbose=10
    --name=%(outfile)s
    --qvalue=%(macs2_max_qvalue)s
    --bdg
    --SPMR
    %(options)s
    %(macs2_options)s
    >& %(outfile)s
    '''
    P.run()

    # compress macs bed files and index with tabix
    for suffix in ('peaks', 'summits'):
        bedfile = outfile + "_" + suffix + ".bed"
        if os.path.exists(bedfile):
            statement = '''
                 bgzip -f %(bedfile)s
                 tabix -f-p bed %(bedfile)s.gz
            '''
        P.run()

    # convert normalized bed graph to bigwig
    # saves 75% of space
    # compressing only saves 60%
    if os.path.exists(outfile + "_treat_pileup.bdg"):
        bedGraphToBigwig(outfile + "_treat_pileup.bdg",
                         os.path.join(PARAMS["annotations_dir"],
                                      "contigs.tsv"),
                         outfile + "_treat_pileup.bw")
    if os.path.exists(outfile + "_control_lambda.bdg"):
        bedGraphToBigwig(outfile + "_control_lambda.bdg",
                         os.path.join(PARAMS["annotations_dir"],
                                      "contigs.tsv"),
                         outfile + "_control_lambda.bw")

    # index and compress peak file
    suffix = 'peaks.xls'
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


def runZinba(infile,
             outfile,
             controlfile,
             action="full",
             fragment_size=None,
             tag_size=None):
    '''run Zinba for peak detection.'''
    E.info("zinba: running action %s" % (action))

    job_threads = PARAMS["zinba_threads"]
    job_options = "-l mem_free=32G"

    # TODO: use closest size or build mapability file
    if not 40 <= tag_size < 60:
        E.warn("tag size out of range of 40-60: %i" % tag_size)

    tag_size = 50
    fragment_size = 200

    mappability_dir = os.path.join(PARAMS["zinba_mappability_dir"],
                                   PARAMS["genome"],
                                   "%i" % tag_size,
                                   "%i" % PARAMS[
                                       "zinba_alignability_threshold"],
                                   "%i" % fragment_size)

    if not os.path.exists(mappability_dir):
        raise OSError(
            "mappability not found, expected to be at %s" % mappability_dir)

    bit_file = os.path.join(PARAMS["zinba_index_dir"],
                            PARAMS["genome"]) + ".2bit"
    if not os.path.exists(bit_file):
        raise OSError("2bit file not found, expected to be at %s" % bit_file)

    options = []
    if controlfile:
        options.append("--control-filename=%(controlfile)s" % locals())

    options = " ".join(options)

    statement = '''
    python %(scriptsdir)s/runZinba.py
           --input-format=bam
           --fdr-threshold=%(zinba_fdr_threshold)f
           --fragment-size=%(fragment_size)i
           --threads=%(zinba_threads)i
           --bit-file=%(bit_file)s
           --zinba-mappability-dir=%(mappability_dir)s
           --zinba-improvement=%(zinba_improvement)f
           --action=%(action)s
           --min-insert-size=%(calling_min_insert_size)i
           --max-insert-size=%(calling_max_insert_size)i
           %(zinba_options)s
           %(options)s
    %(infile)s %(outfile)s
    >& %(outfile)s
    '''

    P.run()

############################################################
############################################################
############################################################


def loadMACS(infile, outfile, bamfile, controlfile=None):
    '''load MACS results into database.

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

    track = P.snip(os.path.basename(infile), ".macs")
    filename_bed = infile + "_peaks.xls.gz"
    filename_diag = infile + "_diag.xls"
    filename_r = infile + "_model.r"
    filename_rlog = infile + ".r.log"
    filename_pdf = infile + "_model.pdf"
    filename_subpeaks = P.snip(infile, ".macs", ) + ".subpeaks.macs_peaks.bed"

    if not os.path.exists(filename_bed):
        E.warn("could not find %s" % filename_bed)
        P.touch(outfile)
        return

    exportdir = os.path.join(PARAMS['exportdir'], 'macs')
    if not os.path.exists(exportdir):
        os.mkdir(exportdir)
    ###############################################################
    # create plot by calling R
    if os.path.exists(filename_r):
        statement = '''R --vanilla < %(filename_r)s > %(filename_rlog)s; mv %(filename_pdf)s %(exportdir)s'''
        P.run()

    ###############################################################
    # filter peaks
    # get thresholds
    max_qvalue = float(PARAMS["macs_max_qvalue"])
    # min, as it is -10log10
    min_pvalue = float(PARAMS["macs_min_pvalue"])

    outtemp = P.getTempFile(".")
    tmpfilename = outtemp.name

    id = 0

    counter = E.Counter()
    with IOTools.openFile(filename_bed, "r") as ins:
        for peak in WrapperMACS.iterateMacsPeaks(ins):

            if peak.fdr > max_qvalue:
                counter.removed_qvalue += 1
                continue
            elif peak.pvalue < min_pvalue:
                counter.removed_pvalue += 1
                continue

            assert peak.start < peak.end

            outtemp.write("\t".join(map(str, (
                peak.contig, peak.start, peak.end,
                id,
                peak.pvalue, peak.fold, peak.fdr,
                peak.start + peak.summit - 1,
                peak.tags))) + "\n")
            id += 1
            counter.output += 1

    outtemp.close()

    ###################################################################
    # output filtering summary
    outf = IOTools.openFile("%s.tsv.gz" % outfile, "w")
    outf.write("category\tcounts\n")
    outf.write("%s\n" % counter.asTable())
    outf.close()

    E.info("%s filtering: %s" % (track, str(counter)))
    if counter.output == 0:
        E.warn("%s: no peaks found" % track)

    ###############################################################
    # load peaks
    shift = getPeakShiftFromMacs(infile)
    assert shift is not None, "could not determine peak shift from MACS file %s" % infile

    E.info("%s: found peak shift of %i" % (track, shift))

    offset = shift * 2

    to_cluster = True

    headers = ",".join((
        "contig", "start", "end",
        "interval_id",
        "pvalue", "fold", "qvalue",
        "macs_summit", "macs_nprobes"))

    tablename = P.toTable(outfile) + "_peaks"

    if controlfile:
        control = "--control-bam-file=%(controlfile)s --control-offset=%(shift)i" % locals()
    else:
        control = ""

    statement = '''python %(scriptsdir)s/bed2table.py 
                           --counter=peaks
                           --bam-file=%(bamfile)s
                           --offset=%(shift)i
                           %(control)s
                           --all-fields 
                           --bed-header=%(headers)s
                           --log=%(outfile)s
                < %(tmpfilename)s
                | python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                       --add-index=contig,start
                       --add-index=interval_id
                       --table=%(tablename)s
                       --allow-empty-file 
                > %(outfile)s'''

    P.run()

    os.unlink(tmpfilename)

    ############################################################
    if os.path.exists(filename_subpeaks):

        headers = ",".join((
            "contig", "start", "end",
            "interval_id",
            "Height",
            "SummitPosition"))

        tablename = P.toTable(outfile) + "_summits"

        # add a peak identifier and remove header
        statement = '''
                    awk '/Chromosome/ {next; } {printf("%%s\\t%%i\\t%%i\\t%%i\\t%%i\\t%%i\\n", $1,$2,$3,++a,$4,$5)}'
                    < %(filename_subpeaks)s
                    | python %(scriptsdir)s/bed2table.py 
                               --counter=peaks
                               --bam-file=%(bamfile)s
                               --offset=%(shift)i
                               %(control)s
                               --all-fields 
                               --bed-header=%(headers)s
                               --log=%(outfile)s
                    | python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                           --add-index=contig,start
                           --add-index=interval_id
                           --table=%(tablename)s
                           --allow-empty-file 
                    > %(outfile)s'''

        P.run()

    ############################################################
    # load diagnostic data
    if os.path.exists(filename_diag):

        tablename = P.toTable(outfile) + "_diagnostics"
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
def loadMACS2(infile, outfile, bamfile, controlfile=None):
    '''load MACS 2 results.

    This method loads only positive peaks. It filters peaks by p-value,
    q-value and fold change and loads the diagnostic data and
    re-calculates peakcenter, peakval, ... using the supplied bamfile.

    It maps the following MACS2 output files to tables:
    * _peaks: MACS2 peaks
    * _summits: MACS2 sub-peaks
    * _regions: MACS2 broad peaks

    This method creates two optional additional files:

       * if the file :file:`<track>_model.r` is present, call R to
         create a MACS peak-shift plot and save it as :file:`<track>_model.pdf`
         in the :file:`export/MACS` directory.

    This method creates :file:`<outfile>.tsv.gz` with the results
    of the filtering.

    '''
    track = P.snip(os.path.basename(infile), ".macs2")
    filename_bed = infile + "_peaks.xls.gz"
    filename_r = infile + "_model.r"
    filename_rlog = infile + ".r.log"
    filename_pdf = infile + "_model.pdf"
    filename_broadpeaks = infile + "_broad_peaks.bed"
    filename_subpeaks = infile + "_summits.bed.gz"

    if not os.path.exists(filename_bed):
        E.warn("could not find %s" % filename_bed)
        P.touch(outfile)
        return

    # Jethro: os.mkdir can't create nested directories
    # AH: use os.makedirs
    exportdir = os.path.join(PARAMS['exportdir'], 'macs2')
    if not os.path.exists(exportdir):
        os.makedirs(exportdir)

    ###############################################################
    # create plot by calling R
    if os.path.exists(filename_r):
        statement = '''
        R --vanilla < %(filename_r)s > %(filename_rlog)s;
        mv %(filename_pdf)s %(exportdir)s'''
        P.run()

    ###############################################################
    # filter peaks - this isn't needed...
    # get thresholds
    max_qvalue = float(PARAMS["macs_max_qvalue"])
    # min, as it is -10log10
    min_pvalue = float(PARAMS["macs_min_pvalue"])

    outtemp = P.getTempFile(".")
    tmpfilename = outtemp.name

    id = 0

    counter = E.Counter()
    with IOTools.openFile(filename_bed, "r") as ins:
        for peak in WrapperMACS.iterateMacs2Peaks(ins):

            if peak.qvalue > max_qvalue:
                counter.removed_qvalue += 1
                continue
            elif peak.pvalue < min_pvalue:
                counter.removed_pvalue += 1
                continue

            assert peak.start < peak.end

            # deliberately not writing out the macs2 assigned peak name...
            outtemp.write("\t".join(map(str, (
                peak.contig, peak.start, peak.end,
                id,
                peak.pvalue, peak.fold, peak.qvalue,
                peak.pileup))) + "\n")
            id += 1
            counter.output += 1

    outtemp.close()

    ###################################################################
    # output filtering summary
    outf = IOTools.openFile("%s.tsv.gz" % outfile, "w")
    outf.write("category\tcounts\n")
    outf.write("%s\n" % counter.asTable())
    outf.close()

    E.info("%s filtering: %s" % (track, str(counter)))
    if counter.output == 0:
        E.warn("%s: no peaks found" % track)

    ###############################################################
    # load peaks
    shift = getPeakShiftFromMacs(infile)
    assert shift is not None,\
        "could not determine peak shift from MACS file %s" % infile

    E.info("%s: found peak shift of %i" % (track, shift))

    offset = shift * 2

    to_cluster = True

    headers = ",".join((
        "contig", "start", "end",
        "interval_id",
        "pvalue", "fold", "qvalue",
        "macs_nprobes"))

    tablename = P.toTable(outfile) + "_peaks"

    if controlfile:
        control = "--control-bam-file=%(controlfile)s --control-offset=%(shift)i" % locals()
    else:
        control = ""

    statement = '''python %(scriptsdir)s/bed2table.py
                           --counter=peaks
                           --bam-file=%(bamfile)s
                           --offset=%(shift)i
                           %(control)s
                           --all-fields
                           --bed-header=%(headers)s
                           --log=%(outfile)s
                < %(tmpfilename)s
                | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
                       --add-index=contig,start
                       --add-index=interval_id
                       --table=%(tablename)s
                       --allow-empty-file
                > %(outfile)s'''

    P.run()

    os.unlink(tmpfilename)

    ############################################################
    if os.path.exists(filename_subpeaks):

        headers = ",".join((
            "contig", "start", "end",
            "interval_id",
            "Height",
            "SummitPosition"))

        tablename = P.toTable(outfile) + "_summits"

        # add a peak identifier and remove header
        statement = '''
                    zcat %(filename_subpeaks)s
                    | awk '/Chromosome/ {next; } {printf("%%s\\t%%i\\t%%i\\t%%i\\t%%i\\t%%i\\n", $1,$2,$3,++a,$4,$5)}'
                    | python %(scriptsdir)s/bed2table.py 
                               --counter=peaks
                               --bam-file=%(bamfile)s
                               --offset=%(shift)i
                               %(control)s
                               --all-fields 
                               --bed-header=%(headers)s
                               --log=%(outfile)s
                    | python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                           --add-index=contig,start
                           --add-index=interval_id
                           --table=%(tablename)s
                           --allow-empty-file 
                    > %(outfile)s'''

        P.run()

    ############################################################
    if os.path.exists(filename_broadpeaks):

        headers = ",".join((
            "contig", "start", "end",
            "interval_id",
            "Height"))

        tablename = P.toTable(outfile) + "_regions"

        # add a peak identifier and remove header
        statement = '''
                    cat %(filename_broadpeaks)s
                    | awk '/Chromosome/ {next; } {printf("%%s\\t%%i\\t%%i\\t%%i\\t%%i\\n", $1,$2,$3,++a,$4)}'
                    | python %(scriptsdir)s/bed2table.py 
                               --counter=peaks
                               --bam-file=%(bamfile)s
                               --offset=%(shift)i
                               %(control)s
                               --all-fields 
                               --bed-header=%(headers)s
                               --log=%(outfile)s
                    | python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                           --add-index=contig,start
                           --add-index=interval_id
                           --table=%(tablename)s
                           --allow-empty-file 
                    > %(outfile)s'''

        P.run()


def loadZinba(infile, outfile, bamfile,
              tablename=None,
              controlfile=None):
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

    track = P.snip(os.path.basename(infile), ".zinba")
    folder = os.path.dirname(infile)

    infilename = infile + ".peaks"

    #######################################################################
    if not os.path.exists(infilename):
        E.warn("could not find %s" % infilename)
    elif P.isEmpty(infile):
        E.warn("no data in %s" % filename)
    else:

        # filter peaks
        offset = getPeakShiftFromZinba(infile)
        assert offset is not None, \
            ("could not determine peak shift from Zinba file %s" %
             infile)

        E.info("%s: found peak shift of %i" % (track, offset))

        if controlfile:
            control = ("--control-bam-file=%(controlfile)s "
                       "--control-offset=%(offset)i" % locals())

        # Steve - Guessing these are actually "peak calls"
        tablename = P.toTable(outfile) + "_peaks"

        headers = "contig,start,end,interval_id,sig,maxloc,maxval,median,qvalue"

        statement = '''cat %(infilename)s
                    | python %(scriptsdir)s/csv_cut.py
                         Chrom Start Stop Sig Maxloc Max Median qValue
        | awk -v FS='\\t' -v OFS='\\t' \
        '/Chrom/ {next; } \
        {$4=sprintf("%%i\\t%%s", ++a, $4); print}'
                    | python %(scriptsdir)s/bed2table.py 
                               --counter=peaks
                               --bam-file=%(bamfile)s
                               --offset=%(offset)i
                               %(control)s
                               --all-fields 
                               --bed-header=%(headers)s
                               --log=%(outfile)s
                    | python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                           --add-index=contig,start
                           --add-index=interval_id
                           --table=%(tablename)s
                           --allow-empty-file 
                    > %(outfile)s'''

        P.run()

        tablename = P.toTable(outfile) + "_summits"

        statement = '''cat %(infilename)s
                    | python %(scriptsdir)s/csv_cut.py Chrom pStart pStop Sig Maxloc Max Median qValue
        | awk -v FS='\\t' -v OFS='\\t' \
        '/Chrom/ {next; } \
        {$4=sprintf("%%i\\t%%s", ++a, $4); print}'
                    | python %(scriptsdir)s/bed2table.py 
                               --counter=peaks
                               --bam-file=%(bamfile)s
                               --offset=%(offset)i
                               %(control)s
                               --all-fields 
                               --bed-header=%(headers)s
                               --log=%(outfile)s
                    | python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                           --add-index=contig,start
                           --add-index=interval_id
                           --table=%(tablename)s
                           --allow-empty-file 
                    > %(outfile)s'''

        P.run()


############################################################
############################################################
############################################################
def runSICER(infile,
             outfile,
             controlfile=None,
             mode="narrow",
             fragment_size=100):
    '''run sicer on infile.'''

    job_options = "-l mem_free=8G"

    workdir = outfile + ".dir"

    try:
        os.mkdir(workdir)
    except OSError:
        pass

    if BamTools.isPaired(infile):
        # output strand as well
        statement = ['''cat %(infile)s
        | python %(scriptsdir)s/bam2bed.py
              --merge-pairs
              --min-insert-size=%(calling_min_insert_size)i
              --max-insert-size=%(calling_max_insert_size)i
              --log=%(outfile)s.log
              --bed-format=6
        > %(workdir)s/foreground.bed''']
    else:
        statement = ['bamToBed -i %(infile)s > %(workdir)s/foreground.bed']

    outfile = os.path.basename(outfile)

    if mode == "narrow":
        window_size = PARAMS["sicer_narrow_window_size"]
        gap_size = PARAMS["sicer_narrow_gap_size"]
    elif mode == "broad":
        window_size = PARAMS["sicer_broad_window_size"]
        gap_size = PARAMS["sicer_broad_gap_size"]
    else:
        raise ValueError("SICER mode unrecognised")

    if controlfile:
        statement.append(
            'bamToBed -i %(controlfile)s > %(workdir)s/control.bed')
        statement.append('cd %(workdir)s')
        statement.append( '''SICER.sh . foreground.bed control.bed . %(genome)s 
                    %(sicer_redundancy_threshold)i
                    %(window_size)i
                    %(fragment_size)i
                    %(sicer_effective_genome_fraction)f
                    %(gap_size)i
                    %(sicer_fdr_threshold)f
                    >& ../%(outfile)s''' )
    else:
        statement.append('cd %(workdir)s')
        statement.append( '''SICER-rb.sh . foreground.bed . %(genome)s 
                    %(sicer_redundancy_threshold)i
                    %(window_size)i
                    %(fragment_size)i
                    %(sicer_effective_genome_fraction)f
                    %(gap_size)i
                    %(sicer_evalue_threshold)f
                    >& ../%(outfile)s''' )

    statement.append('rm -f foreground.bed background.bed')
    statement = '; '.join(statement)

    P.run()

############################################################
############################################################
############################################################


def loadSICER(infile, outfile, bamfile, controlfile=None, mode="narrow"):
    '''load Sicer results.'''

    to_cluster = True

    # build filename of input bedfile
    track = P.snip(os.path.basename(infile), ".sicer")
    sicerdir = infile + ".dir"
    window = PARAMS["sicer_" + mode + "_window_size"]
    gap = PARAMS["sicer_" + mode + "_gap_size"]
    fdr = "%8.6f" % PARAMS["sicer_fdr_threshold"]
    offset = PARAMS["sicer_fragment_size"]

    # taking the file islands-summary-FDR, which contains
    # 'summary file of significant islands with requirement of FDR=0.01'

    bedfile = os.path.join(
        sicerdir,
        "foreground" + "-W" + str(window) +
        "-G" + str(gap) + "-islands-summary-FDR" + fdr)

    assert os.path.exists(bedfile)

    if controlfile:
        control = "--control-bam-file=%(controlfile)s --control-offset=%(offset)i" % locals()

    tablename = P.toTable(outfile) + "_regions"

    headers = "contig,start,end,interval_id,chip_reads,control_reads,pvalue,fold,fdr"

    # add new interval id at fourth column
    statement = '''cat < %(bedfile)s
               | awk '{printf("%%s\\t%%s\\t%%s\\t%%i", $1,$2,$3,++a);
                          for (x = 4; x <= NF; ++x) {printf("\\t%%s", $x)}; printf("\\n" ); }' 
               | python %(scriptsdir)s/bed2table.py
                           --counter=peaks
                           --bam-file=%(bamfile)s
                           --offset=%(offset)i
                           %(control)s
                           --all-fields
                           --bed-header=%(headers)s
                           --log=%(outfile)s
                | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
                       --add-index=contig,start
                       --add-index=interval_id
                       --table=%(tablename)s
                       --allow-empty-file
                > %(outfile)s'''

    P.run()

############################################################


def summarizeSICER(infiles, outfile):
    '''summarize sicer results.'''

    def __get(line, stmt):
        x = line.search(stmt)
        if x:
            return x.groups()

    map_targets = [
        ("Window average: (\d+)", "window_mean", ()),
        ("Fragment size: (\d+) ", "fragment_size", ()),
        ("Minimum num of tags in a qualified window:  (\d+)",
         "window_min", ()),
        ("The score threshold is:  (\d+)", "score_threshold", ()),
        ("Total number of islands:  (\d+)", "total_islands", ()),
        ("chip library size   (\d+)", "chip_library_size", ()),
        ("control library size   (\d+)", "control_library_size", ()),
        ("Total number of chip reads on islands is:  (\d+)",
         "chip_island_reads", ()),
        ("Total number of control reads on islands is:  (\d+)",
         "control_island_reads", ()),
        ("Given significance 0.01 ,  there are (\d+) significant islands",  "significant_islands", ())]

    # map regex to column
    mapper, mapper_header, mapper2pos = {}, {}, {}
    for x, y, z in map_targets:
        mapper[y] = re.compile(x)
        mapper_header[y] = z
        # positions are +1 as first column in row is track
        mapper2pos[y] = len(mapper_header)

    keys = [x[1] for x in map_targets]

    outs = IOTools.openFile(outfile, "w")

    # build headers
    headers = []
    for k in keys:
        if mapper_header[k]:
            headers.extend(["%s_%s" % (k, x) for x in mapper_header[k]])
        else:
            headers.append(k)
    headers.append("shift")

    outs.write("track\t%s" % "\t".join(headers) + "\n")

    for infile in infiles:
        results = collections.defaultdict(list)
        with IOTools.openFile(infile) as f:
            for line in f:
                if "diag:" in line:
                    break
                for x, y in mapper.items():
                    s = y.search(line)
                    if s:
                        results[x].append(s.groups()[0])
                        break

        row = [P.snip(os.path.basename(infile), ".sicer")]
        for key in keys:
            val = results[key]
            if len(val) == 0:
                v = "na"
            else:
                c = len(mapper_header[key])
                if c >= 1:
                    assert len(val) == c, "key=%s, expected=%i, got=%i, val=%s, c=%s" %\
                        (key,
                         len(val),
                            c,
                            str(val), mapper_header[key])
                v = "\t".join(val)
            row.append(v)
        fragment_size = int(row[mapper2pos["fragment_size"]])
        shift = fragment_size / 2

        outs.write("%s\t%i\n" % ("\t".join(row), shift))

    outs.close()

############################################################
############################################################
############################################################


def runPeakRanger(infile, outfile, controlfile):
    '''run peak ranger
    '''

    job_options = "-l mem_free=8G"

    assert controlfile is not None, "peakranger requires a control"

    statement = '''peakranger ranger
              --data <( python %(scriptsdir)s/bam2bam.py -v 0 --method=set-sequence < %(infile)s)
              --control <( python %(scriptsdir)s/bam2bam.py -v 0 --method=set-sequence < %(controlfile)s)
              --output %(outfile)s
              --format bam
              --pval %(peakranger_pvalue_threshold)f
              --FDR %(peakranger_fdr_threshold)f
              --ext_length %(peakranger_extension_length)i
              --delta %(peakranger_delta)f
              --bandwidth %(peakranger_bandwidth)i
              --thread %(peakranger_threads)i
              %(peakranger_options)s
              >& %(outfile)s
    '''

    P.run()

    # usually there is no output
    P.touch(outfile)

############################################################
############################################################
############################################################


def loadPeakRanger(infile, outfile, bamfile, controlfile=None, table_suffix="peaks"):
    '''load peakranger results.'''

    to_cluster = True

    offset = PARAMS["peakranger_extension_length"] // 2

    if controlfile:
        control = "--control-bam-file=%(controlfile)s --control-offset=%(offset)i" % locals()

    # Steve - This was set to _details, but _details = regions (peaks) + summits. Hence changed.
    # Note that Peak ranger reports peaks even when the given fdr cut-off has failed and labels them
    # "fdrFailed" - here, such peaks are explicitely not loaded.
    # AFAIK, Peakranger ranger is optimised to detect peaks arising from point source binding
    # where as Peakranger ccat is optimised to detect regions arising from
    # more diffuse binding events.
    bedfile = infile + "_region.bed"
    headers = "contig,start,end,interval_id,qvalue,strand"
    tablename = P.toTable(outfile) + "_" + table_suffix
    statement = '''python %(scriptsdir)s/bed2table.py
                           --counter=peaks
                           --bam-file=%(bamfile)s
                           --offset=%(offset)i
                           %(control)s
                           --all-fields
                           --bed-header=%(headers)s
                           --log=%(outfile)s
                < <( grep -v "fdrFailed" %(bedfile)s )
                | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
                       --add-index=contig,start
                       --add-index=interval_id
                       --table=%(tablename)s
                       --allow-empty-file
                > %(outfile)s'''
    P.run()

    bedfile = infile + "_summit.bed"
    headers = "contig,start,end,interval_id,qvalue,strand"
    tablename = P.toTable(outfile) + "_summits"
    statement = '''python %(scriptsdir)s/bed2table.py
                           --counter=peaks
                           --bam-file=%(bamfile)s
                           --offset=%(offset)i
                           %(control)s
                           --all-fields
                           --bed-header=%(headers)s
                           --log=%(outfile)s
                < <( grep -v "fdrFailed" %(bedfile)s )
                | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
                       --add-index=contig,start
                       --add-index=interval_id
                       --table=%(tablename)s
                       --allow-empty-file
                > %(outfile)s'''

    P.run()


############################################################
def summarizePeakRanger(infiles, outfile):
    '''summarize peakranger results.'''

    def __get(line, stmt):
        x = line.search(stmt)
        if x:
            return x.groups()

    map_targets = [
        ("# FDR cut off:\s*(\S+)", "fdr_cutoff", ()),
        ("# P value cut off:\s*(\S+)", "pvalue_cutoff", ()),
        ("# Read extension length:\s*(\d+)", "fragment_size", ()),
        ("# Smoothing bandwidth:\s*(\d+)", "smoothing_bandwidth", ()),
    ]

    # map regex to column
    mapper, mapper_header, mapper2pos = {}, {}, {}
    for x, y, z in map_targets:
        mapper[y] = re.compile(x)
        mapper_header[y] = z
        # positions are +1 as first column in row is track
        mapper2pos[y] = len(mapper_header)

    keys = [x[1] for x in map_targets]

    outs = IOTools.openFile(outfile, "w")

    # build headers
    headers = []
    for k in keys:
        if mapper_header[k]:
            headers.extend(["%s_%s" % (k, x) for x in mapper_header[k]])
        else:
            headers.append(k)
    headers.append("shift")

    outs.write("track\t%s" % "\t".join(headers) + "\n")

    for infile in infiles:
        results = collections.defaultdict(list)
        with IOTools.openFile(infile + "_details") as f:
            for line in f:
                if "#region_chr" in line:
                    break
                for x, y in mapper.items():
                    s = y.search(line)
                    if s:
                        results[x].append(s.groups()[0])
                        break

        row = [P.snip(os.path.basename(infile), ".peakranger")]
        for key in keys:
            val = results[key]
            if len(val) == 0:
                v = "na"
            else:
                c = len(mapper_header[key])
                if c >= 1:
                    assert len(val) == c, "key=%s, expected=%i, got=%i, val=%s, c=%s" %\
                        (key,
                         len(val),
                            c,
                            str(val), mapper_header[key])
                v = "\t".join(val)
            row.append(v)
        fragment_size = int(row[mapper2pos["fragment_size"]])
        shift = fragment_size / 2

        outs.write("%s\t%i\n" % ("\t".join(row), shift))

    outs.close()

############################################################
############################################################
############################################################


def runPeakRangerCCAT(infile, outfile, controlfile):
    '''run peak ranger
    '''
    job_options = "-l mem_free=8G"

    assert controlfile is not None, "peakranger requires a control"

    statement = '''peakranger ccat
              --data %(infile)s
              --control %(controlfile)s
              --output %(outfile)s
              --format bam
              --FDR %(peakranger_fdr_threshold)f
              --ext_length %(peakranger_extension_length)i
              --win_size %(ccat_winsize)i
              --win_step %(ccat_winstep)i
              --min_count %(ccat_mincount)i
              --min_score %(ccat_minscore)i
              --thread %(peakranger_threads)i
              %(peakranger_options)s
              >& %(outfile)s
    '''

    P.run()

    # usually there is no output
    P.touch(outfile)

############################################################
############################################################
############################################################


def runSPP(infile, outfile, controlfile):
    '''run spp for peak detection.'''

    job_threads = PARAMS["spp_threads"]
    assert controlfile is not None, "spp requires a control"

    statement = '''
    python %(scriptsdir)s/runSPP.py
           --input-format=bam
           --control-filename=%(controlfile)s
           --fdr-threshold=%(spp_fdr_threshold)f
           --threads=%(spp_threads)i
           --spp-srange-min=%(spp_srange_min)i
           --spp-srange-max=%(spp_srange_max)i
           --bin=%(spp_bin)i
           --spp-z-threshold=%(spp_z_threshold)f
           --window-size=%(spp_window_size)s
           %(spp_options)s
    %(infile)s %(outfile)s
    >& %(outfile)s
    '''

    P.run()

############################################################
############################################################
############################################################


def loadSPP(infile, outfile, bamfile, controlfile=None):
    '''load spp results.'''

    to_cluster = True

    offset = getPeakShift(infile) * 2

    if controlfile:
        control = "--control-bam-file=%(controlfile)s --control-offset=%(offset)i" % locals()

    #
    # Now commented out - the broadpeaks file records arbitrary broad regions of enrichment
    # not controlled by p or q value, it is a preprocessing step in the spp pipeline.
    #
    #bedfile = infile + ".broadpeak.txt"
    # headers="contig,start,end,interval_id"
    #tablename = P.toTable( outfile ) + "_regions"
    # statement = '''
    #            awk '{printf("%%s\\t%%i\\t%%i\\t%%s\\n", $1,$2,$3,++a);}'
    #            < %(bedfile)s
    #            | python %(scriptsdir)s/bed2table.py
    #                       --counter=peaks
    #                       --bam-file=%(bamfile)s
    #                       --offset=%(offset)i
    #                       %(control)s
    #                       --all-fields
    #                       --bed-header=%(headers)s
    #                       --log=%(outfile)s
    #            | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
    #                   --add-index=contig,start
    #                   --add-index=interval_id
    #                   --table=%(tablename)s
    #                   --allow-empty-file
    #            > %(outfile)s'''
    #
    # P.run()

    bedfile = infile + ".narrowpeak.txt"
    headers = "contig,start,end,interval_id,peakval1,qvalue,peakpos"
    tablename = P.toTable(outfile) + "_peaks"
    statement = '''awk '{printf("%%s\\t%%i\\t%%i\\t%%s\\t%%f\\t%%f\\t%%i\\n", $1,$2,$3,++a,$7,$9,$1+$10);}' 
                < %(bedfile)s
                | python %(scriptsdir)s/bed2table.py 
                           --counter=peaks
                           --bam-file=%(bamfile)s
                           --offset=%(offset)i
                           %(control)s
                           --all-fields 
                           --bed-header=%(headers)s
                           --log=%(outfile)s
                | python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                       --add-index=contig,start
                       --add-index=interval_id
                       --table=%(tablename)s
                       --allow-empty-file 
                > %(outfile)s'''

    #
    #  TODO - spp does calculate summit positions, these should be loaded
    #
    #bedfile = infile + ".summits.txt"
    # headers="contig,start,end,interval_id,peakval1,qvalue,peakpos"
    #tablename = P.toTable( outfile ) + "_peaks"
    # statement = '''awk '{printf("%%s\\t%%i\\t%%i\\t%%s\\t%%f\\t%%f\\t%%i\\n", $1,$2,$3,++a,$7,$9,$1+$10);}'
    #            < %(bedfile)s
    #            | python %(scriptsdir)s/bed2table.py
    #                       --counter=peaks
    #                       --bam-file=%(bamfile)s
    #                       --offset=%(offset)i
    #                       %(control)s
    #                       --all-fields
    #                       --bed-header=%(headers)s
    #                       --log=%(outfile)s
    #            | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
    #                   --add-index=contig,start
    #                   --add-index=interval_id
    #                   --table=%(tablename)s
    #                   --allow-empty-file
    #            > %(outfile)s'''

    P.run()

############################################################
############################################################
############################################################
# summarize SPP
############################################################


def summarizeSPP(infiles, outfile):
    '''summarize SPP results by parsing spp output file.'''

    outf = IOTools.openFile(outfile, "w")

    outf.write(
        "track\ttreatment\ttreatment_nreads\tcontrol\tcontrol_nreads\tshift\tfdr\tthreshold\tnpeaks\n")

    for infile in infiles:

        track = P.snip(os.path.basename(infile), ".spp")

        with IOTools.openFile(infile) as inf:
            files, reads = [], []
            for line in inf:
                if line.startswith("opened"):
                    files.append(re.match("opened (\S+)", line).groups()[0])

                elif line.startswith("done. read"):
                    reads.append(
                        re.match("done. read (\d+) fragments", line).groups()[0])

                elif line.startswith("shift\t"):
                    shift = int(re.match("shift\t(\d+)\n", line).groups()[0])

                elif line.startswith("FDR"):
                    fdr, threshold = re.match(
                        "FDR\s+(\S+)\s+threshold=\s+(\S+)", line).groups()

                elif line.startswith("detected_peaks\t"):
                    npeaks = int(
                        re.match("detected_peaks\t(\d+)\n", line).groups()[0])

        outf.write("\t".join(map(str, (track, files[0], reads[0], files[1], reads[1],
                                       shift, fdr, threshold, npeaks))) + "\n")

    outf.close()

############################################################
############################################################
############################################################
# Run broadpeak
############################################################


def createGenomeWindows(genome, outfile, windows):
    to_cluster = True
    statement = '''
        python %(scriptsdir)s/windows2gff.py
            --genome=%(genome)s
            --output-format=bed
            --fixed-width-windows=%(windows)s
            --log=%(outfile)s.log |
        bedtools sort -i stdin |
        gzip > %(outfile)s
        '''
    P.run()


def normalize(infile, larger_nreads, outfile, smaller_nreads):
    threshold = float(smaller_nreads) / float(larger_nreads)

    pysam_in = pysam.Samfile(infile, "rb")
    pysam_out = pysam.Samfile(outfile, "wb", template=pysam_in)

    for read in pysam_in.fetch():
        if random.random() <= threshold:
            pysam_out.write(read)

    pysam_in.close()
    pysam_out.close()
    pysam.index(outfile)

    return outfile


def normalizeFileSize(sample_file, input_file, sample_outfile, input_outfile):
    sample_sam = pysam.Samfile(sample_file, "rb")
    sample_nreads = 0
    for read in sample_sam:
        sample_nreads += 1

    input_sam = pysam.Samfile(input_file, "rb")
    input_nreads = 0
    for read in input_sam:
        input_nreads += 1

    if input_nreads > sample_nreads:
        P.info("INPUT bam has %s reads, SAMPLE bam has %s reads"
               % (input_nreads, sample_nreads))
        P.info("INPUT being downsampled to match SAMPLE")
        input_outfile = normalize(input_file,
                                  input_nreads,
                                  input_outfile,
                                  sample_nreads)
        shutil.copyfile(sample_file, sample_outfile)
        pysam.index(sample_outfile)

        return sample_outfile, input_outfile

    elif sample_nreads > input_nreads:
        P.info("SAMPLE bam has %s reads, INPUT bam has %s reads"
               % (sample_nreads, input_nreads))
        P.info("SAMPLE being downsampled to match INPUT")
        sample_outfile = normalize(sample_file,
                                   sample_nreads,
                                   sample_outfile,
                                   input_nreads)
        shutil.copyfile(input_file, input_outfile)
        pysam.index(input_outfile)

        return sample_outfile, input_outfile

    else:
        P.warn("WARNING: input and sample bamfiles are the same size!!")
        shutil.copyfile(sample_file, sample_outfile)
        pysam.index(sample_outfile)
        shutil.copyfile(input_file, input_outfile)
        pysam.index(input_outfile)

        return sample_outfile, input_outfile


def buildBedFile(infile, outfile):
    statement = ("python %(scriptsdir)s/bam2bed.py %(infile)s"
                 " | sortBed -i stdin"
                 " > %(outfile)s")
    P.run()
    return(outfile)


def createBedgraphFile(infile, outfile, windows_file, overlap):
    statement = ("intersectBed"
                 " -c"
                 " -f %(overlap)s"
                 " -a %(windows_file)s"
                 " -b %(infile)s"
                 " | sortBed -i stdin"
                 " > %(outfile)s")
    P.run()
    return(outfile)


def removeBackground(sample_bedgraph, input_bedgraph, outfile):
    inf = IOTools.openFile(sample_bedgraph, "r").readlines()
    contf = IOTools.openFile(input_bedgraph, "r").readlines()
    outf = IOTools.openFile(outfile, "w")

    for i in range(0, len(inf)):
        line_inf = inf[i].split()
        line_contf = contf[i].split()
        if line_inf[0] == line_contf[0] and line_inf[1] == line_contf[1]:
            if int(line_inf[3]) >= int(line_contf[3]):
                outf.write(str(line_inf[0])
                           + "\t" + str(line_inf[1])
                           + "\t" + str(line_inf[2])
                           + "\t" + str(int(line_inf[3])
                                        - int(line_contf[3])) + "\n")
            else:
                outf.write(str(line_inf[0])
                           + "\t" + str(line_inf[1])
                           + "\t" + str(line_inf[2])
                           + "\t0" + "\n")
        else:
            raise Exception("Windows in '%s' aren't aligned with windows in %s"
                            % (sample_bedgraph, input_bedgraph))
    return outfile


def removeEmptyBins(infile, outfile):
    statement = '''cat %(infile)s
    | awk '$4!=0 {print $1"\t"$2"\t"$3"\t"$4}'
    > %(outfile)s '''
    P.run()


def createBroadPeakBedgraphFile(infiles, outfile, params):
    tmpdir = P.getTempDir("/scratch")
    P.info("Creating tempdir: %s" % tmpdir)
    sample_bed = P.getTempFilename(tmpdir)
    sample_bedgraph = P.getTempFilename(tmpdir)
    input_bed = P.getTempFilename(tmpdir)
    input_bedgraph = P.getTempFilename(tmpdir)
    bgremoved_bedgraph = P.getTempFilename(tmpdir)

    overlap, remove_background = params

    if remove_background == "true":
        sample_file, windows_file, input_file = infiles
        P.info("Background removal:"
               " Input reads will be deducted from sample reads\n"
               "Normalizing size of bamfiles by down-sampling larger file")

        sample_norm = P.snip(sample_file, ".call.bam") + "_normalized.bam"
        input_norm = P.snip(input_file, ".call.bam") + "_normalized.bam"
        sample_norm, input_norm = normalizeFileSize(sample_file,
                                                    input_file,
                                                    sample_norm,
                                                    input_norm)
        P.info("Created normalized SAMPLE bam:\t %s "
               "\nCreated normalized INPUT bam:\t %s"
               % (sample_norm, input_norm))

        P.info("Creating bed file from SAMPLE bam file")
        sample_bed = buildBedFile(sample_norm, sample_bed)

        P.info("Created SAMPLE bed file:\t %s \nCreating SAMPLE bedgraph file.\n"
               "Reads must have coverage >= %s*bin width in order to be binned"
               % (sample_bed, overlap))
        sample_bedgraph = createBedgraphFile(sample_bed,
                                             sample_bedgraph,
                                             windows_file,
                                             overlap)

        P.info("Created SAMPLE bedgraph file:\t %s\nCreating INPUT bed file"
               % sample_bedgraph)
        input_bed = buildBedFile(input_norm, input_bed)

        P.info("Created INPUT bed file:\t %s \nCreating INPUT bedgraph file\n"
               "Reads must have coverage >= %s*bin width in order to be binned"
               % (input_bed, overlap))
        input_bedgraph = createBedgraphFile(input_bed,
                                            input_bedgraph,
                                            windows_file,
                                            overlap)

        P.info("Created INPUT bedgraph file:\t %s "
               "\nSubtracting INPUT bedgraph from SAMPLE bedgraph"
               % input_bedgraph)
        bgremoved_bedgraph = removeBackground(sample_bedgraph,
                                              input_bedgraph,
                                              bgremoved_bedgraph)
        removeEmptyBins(bgremoved_bedgraph, outfile)

        P.info("Created SAMPLE bedgraph with INPUT removed:\t %s" % outfile)
        shutil.rmtree(os.path.abspath(tmpdir))

    else:
        sample_file, windowsfile = infiles
        P.info("Background ignored: bedgraph file does not account for INPUT")
        P.info("Creating bed file from SAMPLE bam file")
        sample_bed = buildBedFile(sample_file, sample_bed)

        P.info("Created SAMPLE bed file:\n %s \nCreating SAMPLE bedgraph file"
               % sample_bed)
        sample_bedgraph = createBedgraphFile(sample_bed,
                                             sample_bedgraph,
                                             windows_file,
                                             overlap)
        removeEmptyBins(sample_bedgraph, outfile)
        P.info("Created SAMPLE bedgraph (without accounting for INPUT):\n"
               "%s" % outfile)

        shutil.rmtree(os.path.abspath(tmpdir))


def runBroadPeak(infile, stub, logfile, genome_size, training_set=False):
    to_cluster = True
    job_options = "-l mem_free=10G"

    if training_set:
        statement = ("BroadPeak -i %(infile)s"
                     " -m %(stub)s"
                     " -b 200"
                     " -g %(genome_size)s"
                     " -t supervised"
                     " -r (training_set)s"
                     " &> %(logfile)s")
    else:
        statement = ("BroadPeak -i %(infile)s"
                     " -m %(stub)s"
                     " -b 200"
                     " -g %(genome_size)s"
                     " -t unsupervised"
                     " &> %(logfile)s")

    P.run()


def summarizeBroadPeak(infiles, outfile, intervals=False):
    if intervals:
        P.warn("Pipeline for summarizing supervised broadpeak runs"
               " has not yet been written... summary file will be"
               " empty")
        IOTools.openFile(oufitle, "w").close()

    else:
        outf = IOTools.openFile(outfile, "w")
        outf.write("track\tnpeaks\tp\tq\ts1\ts2\n")
        for infile in infiles:
            for root, dirs, filenames in os.walk(infile):
                for f in filenames:
                    if f.endswith("_broadpeak_broad_peak_unsupervised.bed"):
                        track = P.snip(
                            f, "_broadpeak_broad_peak_unsupervised.bed")
                        os.symlink(os.path.join(root, f),
                                   track + "_broadpeak_unsupervised.bed")
                        npeaks = len(
                            open(os.path.join(root, f), "r").readlines())
                    elif f.endswith("_parameter_score.txt"):
                        params = open(os.path.join(root, f), "r").readlines()
                        p = params[0].split()[1]
                        q = params[1].split()[1]
                        s1 = params[2].split()[1]
                        s2 = params[3].split()[1]
                    else:
                        continue
            outf.write(track + "\t" + str(npeaks) + "\t" + p +
                       "\t" + q + "\t" + s1 + "\t" + s2 + "\n")
        outf.close()


############################################################
############################################################
############################################################
##
############################################################
def makeIntervalCorrelation(infiles, outfile, field, reference):
    '''compute correlation of interval properties between sets
    '''

    dbhandle = sqlite3.connect(PARAMS["database"])

    tracks, idx = [], []
    for infile in infiles:
        track = P.snip(infile, ".bed.gz")
        tablename = "%s_intervals" % P.quote(track)
        cc = dbhandle.cursor()
        statement = "SELECT contig, start, end, %(field)s FROM %(tablename)s" % locals(
        )
        cc.execute(statement)
        ix = IndexedGenome.IndexedGenome()
        for contig, start, end, peakval in cc:
            ix.add(contig, start, end, peakval)
        idx.append(ix)
        tracks.append(track)
    outs = IOTools.openFile(outfile, "w")
    outs.write("contig\tstart\tend\tid\t" + "\t".join(tracks) + "\n")

    for bed in Bed.iterator(infile=IOTools.openFile(reference, "r")):

        row = []
        for ix in idx:
            try:
                intervals = list(ix.get(bed.contig, bed.start, bed.end))
            except KeyError:
                row.append("")
                continue

            if len(intervals) == 0:
                peakval = ""
            else:
                peakval = str((max([x[2] for x in intervals])))
            row.append(peakval)

        outs.write(str(bed) + "\t" + "\t".join(row) + "\n")

    outs.close()

############################################################
############################################################
############################################################


def buildIntervalCounts(infile, outfile, track, fg_replicates, bg_replicates):
    '''count read density in bed files comparing stimulated versus unstimulated binding.
    '''
    samfiles_fg, samfiles_bg = [], []

    # collect foreground and background bam files
    for replicate in fg_replicates:
        samfiles_fg.append("%s.call.bam" % replicate.asFile())

    for replicate in bg_replicates:
        samfiles_bg.append("%s.call.bam" % replicate.asFile())

    samfiles_fg = [x for x in samfiles_fg if os.path.exists(x)]
    samfiles_bg = [x for x in samfiles_bg if os.path.exists(x)]

    samfiles_fg = ",".join(samfiles_fg)
    samfiles_bg = ",".join(samfiles_bg)

    tmpfile1 = P.getTempFilename(os.getcwd()) + ".fg"
    tmpfile2 = P.getTempFilename(os.getcwd()) + ".bg"

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

        os.unlink(tmpfile2)

    else:
        statement = '''
        python %(toolsdir)s/combine_tables.py 
               --add-file-prefix 
               --regex-filename="[.](\S+)$" 
        %(tmpfile1)s > %(outfile)s
        '''

        P.run()

    os.unlink(tmpfile1)


############################################################
############################################################
############################################################
def loadIntervalsFromBed(bedfile, track, outfile,
                         bamfiles, offsets):
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

    headers = ("AvgVal", "DisttoStart", "GeneList", "Length", "PeakCenter", "PeakVal", "Position",
               "interval_id", "nCpGs", "nGenes", "nPeaks", "nProbes", "nPromoters", "contig", "start", "end")

    tmpfile.write("\t".join(headers) + "\n")

    avgval, contig, disttostart, end, genelist, length, peakcenter, peakval, position, start, interval_id, ncpgs, ngenes, npeaks, nprobes, npromoters = \
        0, "", 0, 0, "", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

    mlength = int(PARAMS["calling_merge_min_interval_length"])

    c = E.Counter()

    # count tags
    for bed in Bed.iterator(IOTools.openFile(infile, "r")):

        c.input += 1

        if "name" not in bed:
            bed.name = c.input

        # remove very short intervals
        if bed.end - bed.start < mlength:
            c.skipped_length += 1
            continue

        if replicates:
            npeaks, peakcenter, length, avgval, peakval, nprobes = \
                PipelineChipseq.countPeaks(
                    bed.contig, bed.start, bed.end, samfiles, offsets)

            # nreads can be 0 if the intervals overlap only slightly
            # and due to the binning, no reads are actually in the overlap region.
            # However, most of these intervals should be small and have already be deleted via
            # the merge_min_interval_length cutoff.
            # do not output intervals without reads.
            if nprobes == 0:
                c.skipped_reads += 1

        else:
            npeaks, peakcenter, length, avgval, peakval, nprobes = (1,
                                                                    bed.start +
                                                                    (bed.end -
                                                                     bed.start) // 2,
                                                                    bed.end -
                                                                    bed.start,
                                                                    1,
                                                                    1,
                                                                    1)

        c.output += 1
        tmpfile.write("\t".join(map(str, (avgval, disttostart, genelist, length,
                                          peakcenter, peakval, position, bed.name,
                                          ncpgs, ngenes, npeaks, nprobes, npromoters,
                                          bed.contig, bed.start, bed.end))) + "\n")

    if c.output == 0:
        E.warn("%s - no intervals")

    tmpfile.close()

    tmpfilename = tmpfile.name
    tablename = "%s_intervals" % track.asTable()

    statement = '''
    python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty-file
              --add-index=interval_id 
              --table=%(tablename)s
    < %(tmpfilename)s 
    > %(outfile)s
    '''

    P.run()
    os.unlink(tmpfile.name)

    L.info("%s\n" % str(c))


def makeReproducibility(infiles, outfile):
    '''compute overlap between intervals.

    Compute pairwise overlap between all sets in a group
    of :term:`bed` formatted files.
    '''

    to_cluster = True

    if os.path.exists(outfile):
        # note: update does not work due to quoting
        os.rename(outfile, outfile + ".orig")
        options = "--update=%s.orig" % outfile
    else:
        options = ""

    infiles = " ".join(infiles)

    # note: need to quote track names
    statement = '''
        python %(scriptsdir)s/diff_bed.py --pattern-identifier='([^/]+).bed.gz' %(options)s %(infiles)s 
        | awk -v OFS="\\t" '!/^#/ { gsub( /-/,"_", $1); gsub(/-/,"_",$2); } {print}'
        > %(outfile)s
        '''

    P.run()


############################################################
############################################################
############################################################
def runScripture(infile, outfile,
                 contig_sizes,
                 mode="narrow"):
    '''run scripture on infile.'''

    job_options = "-l mem_free=8G"

    samfile = pysam.Samfile(infile, "rb")
    contigs = samfile.references

    s = '''scripture
                   -task chip 
                   -trim
                   -minMappingQuality %%(scripture_min_mapping_quality)f
                   -windows 200
                   -fullScores
                   -sizeFile %%(contig_sizes)s
                   -alpha %%(scripture_fdr)f
                   -alignment %%(infile)s
                   -chr %(contig)s
                   -out %%(outfile)s.data.%(contig)s
                   >& %%(outfile)s.log.%(contig)s
    '''

    statements = [s % {'contig': x} for x in contigs]
    P.run()

    statements = None

    # collect all results into a single bed file
    statement = '''cat %(outfile)s.data.*.scores | gzip > %(outfile)s.bed.gz'''
    P.run()

    statement = '''cat %(outfile)s.log.* > %(outfile)s'''
    P.run()

    statement = '''rm -f %(outfile)s.data.*'''
    P.run()

    statement = '''rm -f %(outfile)s.log.*'''
    P.run()


def loadScripture(infile, outfile, bamfile, controlfile=None):
    '''load scripture peaks.'''

    # Note: not sure if the following will work for
    #       paired end data.

    # no offset
    #    offset = getPeakShift( infile ) * 2
    offset = 0

    if controlfile:
        control = "--control-bam-file=%(controlfile)s --control-offset=%(offset)i" % locals()

    bedfile = infile + ".bed.gz"

    headers = "contig,start,end,interval_id,score,pvalue,score2,score3,score4"
    tablename = P.toTable(outfile) + "_peaks"
    statement = '''zcat %(bedfile)s 
                | awk '{printf("%%s\\t%%i\\t%%i\\t%%s\\t%%f\\t%%f\\t%%f\\t%%f\\t%%f\\n", $1,$2,$3,++a,$5,$7,$8,$9,$10);}' 
                | python %(scriptsdir)s/bed2table.py 
                           --counter=peaks
                           --bam-file=%(bamfile)s
                           --offset=%(offset)i
                           %(control)s
                           --all-fields 
                           --bed-header=%(headers)s
                           --log=%(outfile)s
                | python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                       --add-index=contig,start
                       --add-index=interval_id
                       --table=%(tablename)s
                       --allow-empty-file 
                > %(outfile)s'''

    P.run()
