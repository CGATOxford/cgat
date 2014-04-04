'''
cgat_script_template.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------
Convert BAM files into BED files

Usage
-----

Example::

   python bam2bed.py <input.bam> --options

Type::

   python cgat_script_template.py --help

for command line help.

Documentation
-------------
This script takes a :term:`bam` file as input and counts the number of reads 
in windows across each chromsome.  Output is a :term:`bed` formatted file 
for each interval with the relevant counts of reads mapping to each interval.  
Reads can be selected on the basis of :term:`tags`, such as NH or NM and a 
threshold for the proportion of reads mapping to an interval before 
it is output.  Read counting is done using mapping position within an interval.

The :term:`bam` file must have an :term:`index` file (.bai) in the same 
directory. A contig file can be provided with the --contig-file option, or a 
new one can be generated using a specific genome and directory with the 
--genome and --genome-dir.
The default action is to create a dictionary of contig sizes from the bam file
header info, which requires no additional user input.

The original bam2bed.py functionality was to convert BAM files to BED files 
over a user-specified interval.  This functionality has been preserved in 
this version.

For example::

  samtools view example.bam

  READ1    163    1     13040   15      76M     =       13183   219     ...
  READ1    83     1     13183   7       76M     =       13040   -219    ...
  READ2    147    1     13207   0       76M     =       13120   -163    ...

  python bam2bed.py example.bam 

   1       13039   13115   READ1     15      +
   1       13119   13195   READ2     0       +
   1       13182   13258   READ1     7       -
   1       13206   13282   READ2     0       -
Options
^^^^^^^

+-------------------+--------------------------------------------------------+
|--genome-dir       |location of genome .fasta files to construct contig     |
|                   |dictionary                                              | 
+-------------------+--------------------------------------------------------+
|--genome           |genome to extract contig data for                       |
+-------------------+--------------------------------------------------------+
|--contig-file      |contig file containing contig IDs and sizes             |
+-------------------+--------------------------------------------------------+
|--condition        |used to specify whether a range or categorical filter   |
|                   |is used to select reads for counting.  The default      |
|                   |behaviour with `cat` is to output by contig.            |
+-------------------+--------------------------------------------------------+
|--interval         | size of intervals over which to count reads            |
+-------------------+--------------------------------------------------------+
|--tag              |SAM format tag used to count reads, e.g. to only ouput  |
|                   |reads that are multi-mapping use the NH tag and supply  |
|                   |a float value to the '--threshold' option               |
+-------------------+--------------------------------------------------------+
|--threshold        |The cut-off value above which regions with this number  |
|                   |of reads mapping to it will be ouput.  e.g. for regions |
|                   |where > 1% of reads are multimapping use --tag NH and   |
|                   |--threshold 0.1                                         |
+-------------------+--------------------------------------------------------+
|--region, -r       |output read intervals that overlap a specified region   |
+-------------------+--------------------------------------------------------+
|--merge-pairs, -m  |merge paired-end reads and output interval for entire   |
|                   |fragment                                                |
+-------------------+--------------------------------------------------------+
|--max-insert-size  |only merge if insert size is less than specified no. of |
|                   |bases                                                   |
+-------------------+--------------------------------------------------------+
|--min-insert-size  |only merge if insert size is greater than specified     |
|                   |no. of bases                                            |
+-------------------+--------------------------------------------------------+


For example to output only 1000bp intervals that have >10% of reads that 
are multi-mapping::

  python bam2bed.py example.bam --condition range --interval 1000 --tag NH --threshold 0.1
  
To just count reads within 100bp intervals you can supply the command thus::

  python bam2bed.py example.bam --condition range --interval 100

The output files are a :term:bedGraph format file with each interval that has 
>0 reads aligned to it and a tab-delimtted per-contig :term:summary file that 
contains::

  Contig: Intervals: Total: Min: Max: Mean:   Median:  25th_Q:  75th_Q:
  chr1    40         100189 4    6971 2504.72 2942.00  1040.00  3521.00
  chr10   48         7358   1    1136 153.29  21.00    9.00     47.00
  chr11   16         3502   1    3312 218.88  4.00     1.00     10.00

To use the original function of bam2bed.py use the `--region` option.

The default behaviour is to output the sum of reads across regions in 10kbp windows.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import pysam
import itertools
import math
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Bed as Bed
try:
    import pyximport
    pyximport.install(build_in_temp=False)
    import  _bam2bed
except ImportError:
    import CGAT._bam2bed as _bam2bed

def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--genome-dir", dest="gDir", type="string",
                      help = "genome directory containing fasta format files")
    
    parser.add_option("--genome", dest="genome", type="string",
                      help = "the genome to which reads were aligned")

    parser.add_option("--contig-file", dest="contig", type="string",
                      help = "list of contigs and their sizes in tsv format")

    parser.add_option("--condition", dest="condition", type="string", 
                      help="either `range` or `cat`, `range` is used with --interval.")
    
    parser.add_option("--interval", dest="interval", type="int",
                      help = "the interval size to count reads in. [default=%default]")
    
    parser.add_option("--tag", dest="tag", type = "string",
                      help = "outputs regions with >:term:`threshold` reads containing `tag`")

    parser.add_option("--threshold", dest="threshold", type="float",
                      help = "a threshold for reads matching :term:`tag` in an interval. [default=%default]")
    
    parser.add_option("-r", "--region", dest="region", type="string", 
                      help="output read intervals that overlap samtools region striing [default=%default]. ")

    parser.add_option("-m", "--merge-pairs", dest="merge_pairs", action="store_true",
                      help="merge paired-end reads and output interval for entire fragment [default=%default]. ")

    parser.add_option("--max-insert-size", dest="max_insert_size", type="int",
                      help="only merge paired-end reads if they are less than # bases apart."
                      " 0 turns off this filter [default=%default]. ")

    parser.add_option("--min-insert-size", dest="min_insert_size", type="int",
                     help="only merge paired-end reads if theya re at least # bases apart. "
                     " 0 turns off this filter [default=%default]. ")

    parser.add_option("--bed-format", dest="bed_format", type="choice",
                      choices=('3', '4', '5', '6'),
                      help="bed format to output. Only used with --merge-pairs option. "
                      " [default=%default]")

    parser.set_defaults(gDir = None,
                        genome = None, 
                        contig = None,
                        condition = None,
                        threshold = 0.01,
                        interval = 10000,
                        tag = None,
                        region = None,
                        call_peaks = None,
                        merge_pairs = None,
                        min_insert_size = 0,
                        max_insert_size = 0,
                        bed_format = '6')
    
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if len(args) == 0:
        args.append("-")
        
    options.bed_format = int(options.bed_format)
    
    # test presence of input files and bam index

    if sys.argv[1] and not os.path.exists(sys.argv[1]):
        raise IOError("file %s does not exist" % sys.argv[1])
    else:
        pass

    bam_index = "%s.bai" % sys.argv[1]
    file_prefix = str(sys.argv[1]).rstrip(".bam")

    if not os.path.exists(bam_index):
        raise IOError("there is no index file for %s" % sys.argv[1])

    if options.interval and options.interval < 100:
        E.warn("Minimum interval size is 100bp! Overriding input to 100bp")
        options.interval = 100      
    elif options.condition == "cat":
        options.interval = None

    # if interval is given, assume range is also required for condition    
    
    if options.interval and not options.condition:
        options.condition = "range"       
    
    # currently assumes that the input file is bam format
    # will need to change to either accept sam format
    # or automatically detect format

    # read in whole bam file 
        
    E.info("reading input file %s" % sys.argv[1])    
    samfile = pysam.Samfile(sys.argv[1], "rb")

       
    def contigFileMaker(genome_dir, genome):
        '''Make a contig file unless a file is passed through options.contig.
        This code is ripped from pipeline_annotations function buildContigSizes'''

        genome_file = open("%s/%s.fasta" % (genome_dir, genome), "r")
        prefix = "%s/%s" % (genome_dir, genome)
        fasta  = IndexedFasta.IndexedFasta(prefix)
        
        E.info("Creating contig file")
        with open("contig.tsv", "w") as outs:
            
            for contig, size in fasta.getContigSizes(with_synonyms = False).iteritems():
                outs.write("%s\t%i\n" % (contig, size))
        
        genome_file.close
        

    def summaryInfo(bedfile):
        '''Collects and summarises information from the bedfile output'''
        
        res_dict = {}
        contig_dict = {}
        cur_chr = 'chr1'
        for interval in Bed.iterator(bedfile):
            if interval.contig == cur_chr:
                res_dict["%s:%i-%i" % (interval.contig, interval.start, interval.end)] = float(interval.fields[0])
                contig_dict[cur_chr] = res_dict
            elif interval.contig != cur_chr:
                res_dict = {}
                res_dict[("%s:%i-%i" % (interval.contig, interval.start, interval.end))] = float(interval.fields[0])
                contig_dict[interval.contig] = res_dict
                cur_chr = interval.contig
            
        # summary info:
        # Number of intervals per contig, highest number of reads across intervals
        # lowest number of reads across intervals, average number of reads across intervals
        # Total number of reads across intervals, median number of reads across intervals
        # 25th and 75th quartiles across intervals


        for contig in contig_dict:
            counts = contig_dict[contig].values()
            num_regions = len(counts)
            max_in_contig = max(counts)
            min_in_contig = min(counts)
            mean_in_contig = sum(counts)/float(len(counts))
            total_in_contig = sum(counts)
            counts.sort()
            median_in_contig = counts[int(len(counts)*0.5)]
            u75 = counts[int(len(counts)*0.75)]
            l75 = counts[int(len(counts)*0.25)]
            yield (contig, num_regions, total_in_contig, min_in_contig, 
                   max_in_contig, mean_in_contig, median_in_contig, l75, u75)
            
        
    # check for contig file or genome directory
    # make a contig file if genome and directory are specified
    # otherwise use the file header information

    if options.contig == None and options.gDir != None:
        if options.gDir and not os.path.exists(options.gDir):
            raise IOError("Genome directory %s does not exist" % options.gDir)
        else:
            contigFileMaker(options.gDir, options.genome)
            E.info("Generating contig dictionary")
            contig_dict = {}
            with open("contig.tsv","r") as contig_file:
                for line in contig_file.readlines():
                    contig_dict[(line.rstrip("\n").split("\t"))[0]] = int((line.rstrip("\n").split("\t"))[1]) 

    else:
        if options.contig and not os.path.exists(options.contig):
            raise IOError("contig file %s does not exist" % options.contig)
        else: 
            E.info("Making contig dictionary from input file %s header" % sys.argv[1])
            contig_it = itertools.izip(samfile.references, samfile.lengths)
            contig_dict = {}
            for i  in contig_it:
                contig_dict[i[0]] = i[1]

    # This version of the counter uses a Profiler class to iterate over a bam file and check
    # each read against a condition to see if it should be counted.
    # The initial simple implementation is counting per contig.
    # Will generalize out to include user defined options

    class Result(object):
        ''' The result class wraps up useful read information for the Aggregator class functions '''

        def __init__(self, contig, start, end, idx, strand, tag = None):
            self.contig, self.start, self.end = contig, start, end
            self.idx, self.strand, self.tag =  idx, strand, tag

    class ReadCounter(object):
        ''' A class that determines if a read should be counted '''

        def __init__(self, read):
            self.qname, self.pos, self.aend, self.inferred_length, self.tid =\
                read.qname, read.pos, read.aend, read.inferred_length, read.tid

        def __getattr__(self, item):
            return item

        def checkCount(self, item, condition= None):
            ''' Checks whether the read meets a categorical condition, returns True or False '''

            if condition == None:
                return True
            elif self.__getattr__(item) == condition:
                return True
            else:
                return False

        def checkRange(self, item, condition = None):
            ''' Checks whether a read attribute falls within a given range
             Range is defined by passing in a tuple of (`lower limit`, `upper limit`)
             Returns True or False
            '''

            if condition == None:
                return True
            elif self.__getattr__(item) >= condition[0] and self.__getattr__(item) <= condition[1]:
                return True
            else:
                return False
                   

    class Profiler(object):
        ''' A class that takes a list of ReadCounter instances and a sam file and
        yields profiles of counts'''

        def __init__(self, samfile, condition = None, k = 10000, tag = None, threshold = None):
            # The samfile is a standard pysam samfile object
            
            self.samfile = samfile
            self.condition = condition
            self.k = k
            self.tag = tag
            self.threshold = threshold

        def __iter__(self):
            # Make me an interator
            return self.TheCounter(self.samfile, self.condition, self.tag, self.k, self.threshold)
                
        def TheCounter(self, samfile, condition = None, tag = None, k = 10000, threshold = None):

            contig_it = itertools.izip(samfile.references, samfile.lengths)
            contig_dict = {}

            for i in contig_it:
                contig_dict[i[0]] = i[1]


            # Count reads in a bamfile
            # initialise the variables

            running_total = 0
            last_tid = 0
            res = []
            for read in samfile:
                current_tid = read.tid
                current_pos = read.pos

            # assign strandedness to each read
            # not used currently - implemented for future use

                if read.is_reverse:
                    strand = '+'
                else:
                    strand = '-'

            # count reads based on input - only using range at the moment
            # need specific implementations for `cat`

                if self.condition == "cat":
                    tocount = ReadCounter(read).checkCount(read.tid)
                elif self.condition == "range":
                    tocount = ReadCounter(read).checkRange(read.pos)
                else:
                    tocount = ReadCounter(read).checkCount(read.tid, current_tid)
            
                # check reads are on the same contig, if not change the contig variable and
                # reinitialise the variables
                # pass results to Aggregator class functions

                if current_tid != last_tid:
                    if self.condition == "range" and self.tag == None:
                        yield Aggregator(res).sumInterval(res, self.k, 
                                                          contig_dict[samfile.getrname(last_tid)])
                    elif self.condition == "cat":
                        yield (samfile.getrname(current_tid), 1, contig_dict[samfile.getrname(current_tid)], 
                               running_total)
                    elif self.condition == "range" and self.tag != None:
                        E.info("%s: Checking for intervals with > %0.2f reads with tag %s" % 
                               (samfile.getrname(current_tid), threshold, tag))
                        yield Aggregator(res).tagsInterval(res, self.k, 
                                                           contig_dict[samfile.getrname(last_tid)], 
                                                           self.threshold)

                    last_tid = current_tid
                    last_pos = current_pos
                    running_total = 0 
                    res = []
                else:
                    pass

        # begin counting reads from the list of Counters 
        # the running total becomes the read index - currently not used

                if tocount:
                    running_total += 1
                    if tag == None:
                        res.append(Result(read.tid, read.pos, read.aend, running_total, strand))
                    else:
                        res.append(Result(read.tid, read.pos, read.aend, running_total, strand, read.opt(tag)))
                else:
                    pass

    class Aggregator(object):
        ''' A class to perform functions to be called by Profiler.
        Options are: sum over interval and tagged reads over intervals
        '''

        def __init__(self, results, func = None):
            self.results = results
            self.func = func

        def sumInterval(self, results, k, contig_size):
            ''' Sums read counts over each interval of size k'''

            interval = 1
            last_interval = 1
            intsum = 0

            # check if read alignment positions falls within in interval
            # increment interval by k until it does, then count
            for record in results:
                if record.start <= interval:
                    intsum += 1
                    last_interval = interval
                else:
                    while record.start > interval:
                        if interval <= contig_size:
                            last_interval = interval
                            interval += k
                        # check the incremented interval isn't larger than the contig size
                        # reset to contig size if it is
                        elif interval > contig_size:
                            interval = contig_size
                    if interval > contig_size:
                        interval = contig_size
                    else:
                        pass

                    # return total count for each new interval that contains aligned reads
                    if interval > last_interval and intsum > 0:
                        yield samfile.getrname(record.contig), last_interval, interval, intsum
                        intsum = 0

                      
        def tagsInterval(self, results, k,  contig_size, threshold):
            '''Counts the number of reads that match a tag within a region and highlights
            regions > threshold tagged reads
            E.g. output regions with a large proportion of multi-mapping reads. 
            Default is intervals with >1% reads :term:`tag`>1
            '''

            reads = 0
            counter = 0
            interval = 1
            last_interval = 1
            # use read mapping position to allocate to intervals
            # only count reads that have a tag value > 1
            # this will capture NH > 1 and NM > 1
            # TODO: implement this as a user-defined threshold
            for record in results:
                if record.start <= interval:
                    reads += 1
                    last_interval = interval
                    if record.tag > 1:
                        counter += 1
                    else:
                        pass
                else:
                    while record.start > interval:
                        if interval <= contig_size:
                            last_interval = interval
                            interval += k
                        elif interval > contig_size:
                            interval = contig_size

                if interval > contig_size:
                    interval = contig_size
                else:
                    pass
                # only yield results if the proportion of reads in an interval exceed a threshold
                if counter > 1 and interval > last_interval and float(counter)/float(reads) >= threshold:
                    yield samfile.getrname(record.contig), last_interval, interval, counter
                    reads = 0
                    counter = 0

    # The following code is the original implementation of bam2bed.py
    # It should be wrapped up in a function.  It is exectuted if the --region or
    # --merge-pairs options are provided

    if (options.region != None) or (options.merge_pairs != None):

        if options.merge_pairs is not None:
            counter = _bam2bed.merge_pairs(samfile,
                                           options.stdout,
                                           min_insert_size=options.min_insert_size,
                                           max_insert_size=options.max_insert_size,
                                           bed_format=options.bed_format)

            options.stdlog.write("category\tcounts\n%s\n" % counter.asTable())

        else:
            if options.region is not None:
                if args[0] == "-":
                    raise ValueError("can't use region with a file from stdin")
                it = samfile.fetch(region=options.region)
            else:
                # use until_eof. Files from stdin have no index
                it = samfile.fetch(until_eof=True)

            # more comfortable cigar parsing will
            # come with the next pysam release
            BAM_CMATCH = 0
            BAM_CDEL = 2
            BAM_CREF_SKIP = 3
            take = (BAM_CMATCH, BAM_CDEL, BAM_CREF_SKIP)
            outfile = options.stdout
                        
            for read in it:
                if read.is_unmapped:
                    continue

                t = 0
                for op, l in read.cigar:
                    if op in take:
                        t += l

                start = read.pos
                if read.is_reverse:
                    strand = "-"
                else:
                    strand = "+"
                # IMS: converted rname to reference name
                outfile.write("%s\t%d\t%d\t%s\t%d\t%c\n" %
                              (samfile.getrname(read.rname),
                               read.pos,
                               read.pos + t,
                               read.qname,
                               read.mapq,
                               strand))

    else:

        # Write out a summary file
        # TODO: have as a user defined option to switch on/off

        
        if options.interval and options.tag:
            bedfileOut = IOTools.openFile("%s-%s-interval.bed" % (file_prefix, options.tag), "w")
            bedfileOut.write('track type=bedGraph name="%s intervals" description="Intervals with %s > %0.2f%%" visibility=full priority=20\n' % (options.tag, options.tag, round(options.threshold*100, 2)))
            summaryfile = IOTools.openFile("%s-%s.summary.txt" % (file_prefix, options.tag), "w")        
            summaryfile.write("Contig:\tIntervals:\tTotal:\tMin:\tMax:\tMean:\tMedian:\t25th_quartile:\t75th_quartile:\n")

            E.info("Writing bedGraph file %s-%s-interval.bed" % (file_prefix, options.tag))
            prof1 = Profiler(samfile, options.condition, options.interval, options.tag, options.threshold)
            prof2 = prof1
            for res in prof1:
                for i in res:
                    bedfileOut.write("%s\t%i\t%i\t%0.2f\n" %(i[0], i[1], i[2], i[3]))                    

            bedfileOut.close()
            bedfileIn = IOTools.openFile("%s-%s-interval.bed" % (file_prefix, options.tag), "r")
            for result in summaryInfo(bedfileIn):
                summaryfile.write("%s\t%i\t%i\t%i\t%i\t" % (result[0], result[1], result[2], result[3], result[4]))
                summaryfile.write("%0.2f\t%0.2f\t%0.2f\t%0.2f\n" % (result[5], result[6], result[7], result[8]))
            bedfileIn.close()

        elif options.interval and not options.tag:
            bedfileOut = IOTools.openFile("%s-interval.bed" % (file_prefix), "w")
            bedfileOut.write('track type=bedGraph name="%s intervals" description="Intervals with %s > %0.2f%%" visibility=full priority=20\n' % (options.tag, options.tag, round(options.threshold*100, 2)))
            summaryfile = IOTools.openFile("%s.summary.txt" % (file_prefix), "w")        
            summaryfile.write("Contig:\tIntervals:\tTotal:\tMin:\tMax:\tMean:\tMedian:\t25th_quartile:\t75th_quartile:\n")
            E.info("Writing bedGraph file %s-interval.bed" % (file_prefix))
            prof1 = Profiler(samfile, options.condition, options.interval, threshold = options.threshold)
            prof2 = prof1
            for res in prof1:
                for i in res:
                    bedfileOut.write("%s\t%i\t%i\t%0.2f\n" %(i[0], i[1], i[2], i[3]))                    

            bedfileOut.close()
            bedfileIn = IOTools.openFile("%s-interval.bed" % (file_prefix), "r")
            for result in summaryInfo(bedfileIn):
                summaryfile.write("%s\t%i\t%i\t%i\t%i\t" % (result[0], result[1], result[2], result[3], result[4]))
                summaryfile.write("%0.2f\t%0.2f\t%0.2f\t%0.2f\n" % (result[5], result[6], result[7], result[8]))
            bedfileIn.close()
        
        else:
            bedfileOut = IOTools.openFile("%s-cat.bed" % (file_prefix), "w")
            bedfileOut.write('track type=bedGraph name="category" description="Categorical" visibility=full priority=20\n')
            for res in Profiler(samfile, options.condition):
                bedfileOut.write("%s\t%i\t%i\t%i\n" %(res[0], res[1], res[2], res[3]))

    # write footer and output benchmarkinformation.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
