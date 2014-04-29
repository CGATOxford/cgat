'''
bam2feature.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script takes a :term:`bam` file as input and counts the number of reads
in features.  Output is a :term:`bed` formatted file
for each interval with the relevant counts of reads mapping to each interval.
Reads can be selected on the basis of :term:`tags`, such as NH or NM and a
threshold for the proportion of reads mapping to an interval before
it is output.  Read counting is done using mapping position within an interval.

The :term:`bam` file must have an :term:`index` file (.bai) in the same
directory.

Options
------

+-------------------+--------------------------------------------------------+
|--condition        |used to specify whether a range or categorical filter   |
|                   |is used to select reads for counting.  The default      |
|                   |behaviour with `cat` is to output by contig.            |
+-------------------+--------------------------------------------------------+
|--interval         |size of intervals over which to count reads             |
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
|--soft-clip        |Counts reads that are soft-clipped based on the CIGAR   |
|                   |string. Combine with `--threshold` to select regions    |
|                   |with proportion reads soft-clipped > threshold          |
+-------------------+--------------------------------------------------------+
|--hard-clip        |Counts reads that hard-clipped analagously to the       |
|                   |`--soft-clipped` option.                                |
+-------------------+--------------------------------------------------------+


Usage
-----

For example to output only 1000bp intervals that have >10% of reads that
are multi-mapping::

  python bam2feature.py example.bam
         --condition range
         --interval 1000
         --tag NH
         --threshold 0.1
         > example.bed
  
To just count reads within 100bp intervals you can supply the command thus::

  python bam2feature.py example.bam
  --condition range
  --interval 100
  > example.bed

The output files are a :term:bedGraph format file with each interval that has
>0 reads aligned to it::

  track type=bedGraph name="NH intervals" ...
  chr1   24690001  24700001   7.00
  chr11  58110001  58120001   12.00
  chr12  63790001  6938000    12.00
  chr17  24860001  2487000    13.00
  chr18  67400001  6741000    13.00

The default behaviour is to output the sum of reads across regions in
10kbp windows

Type::

   python cgat_script_template.py --help

for command line help.


Command line options
--------------------

'''

import os
import sys
import pysam
import itertools
import CGAT.Experiment as E

    # This version of the counter uses a Profiler class to iterate
    # over a bam file and check each read against a condition
    # to see if it should be counted.
    # The initial simple implementation is counting per contig.
    # Will generalize out to include user defined options


class Result(object):

    '''
    The result class wraps up useful read information
    for the Aggregator class functions '''

    def __init__(self, contig,
                 start, end,
                 idx, strand,
                 tag=None,
                 cigar=None):

        self.contig, self.start, self.end = contig, start, end
        self.idx = idx
        self.strand = strand
        self.tag = tag
        self.cigar = cigar


class ReadCounter(object):

    ''' A class that determines if a read should be counted '''

    def __init__(self, read):
        self.qname,  self.pos,
        self.aend, self.inferred_length,
        self.tid, self.cigarstring = read.qname, read.pos,
        read.aend, read.inferred_length, read.tid, read.cigarstring

    def __getattr__(self, item):
        return item

    def checkCount(self, item, condition=None):
        '''
        Checks whether the read meets a categorical condition,
        returns True or False
        '''

        if condition is None:
            return True
        elif self.__getattr__(item) == condition:
            return True
        else:
            return False

    def checkRange(self,
                   item,
                   condition=None):
        ''' Checks whether a read attribute falls within a given range
        Range is defined by passing in
        a tuple of (`lower limit`, `upper limit`)
        Returns True or False
        '''

        if condition is None:
            return True
        elif (self.__getattr__(item) >= condition[0]
              and self.__getattr__(item) <= condition[1]):
            return True
        else:
            return False


class Profiler(object):

    ''' 
    A class that takes a list of ReadCounter instances
    and a sam file and yields profiles of counts
    '''

    def __init__(self,
                 samfile,
                 condition=None,
                 k=10000,
                 tag=None,
                 threshold=0.01):

        # The samfile is a standard pysam samfile object

        self.samfile = samfile
        self.condition = condition
        self.k = k
        self.tag = tag
        self.threshold = threshold

    def __iter__(self):
        # Make me an interator
        return self.TheCounter(samfile=self.samfile,
                               condition=self.condition,
                               k=self.k,
                               tag=self.tag,
                               threshold=self.threshold)

    def TheCounter(self,
                   samfile,
                   condition=None,
                   k=10000,
                   tag=None,
                   threshold=0.01):

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
            if read.tid != -1:
                current_tid = read.tid
                contig_id = samfile.getrname(current_tid)
                last_contig = samfile.getrname(last_tid)

                # assign strandedness to each read
                # not used currently - implemented for future use

                if read.is_reverse:
                    strand = '+'
                else:
                    strand = '-'

                # count reads based on input - only using range at the moment
                # need specific implementations for `cat`

                    if self.condition == "cat":
                        tocount = ReadCounter(
                            read=read).checkCount(item=read.tid)
                    elif self.condition == "range":
                        tocount = ReadCounter(read).checkRange(item=read.pos)
                    else:
                        tocount = ReadCounter(read).checkCount(
                            item=read.tid,
                            condition=current_tid)

            # check reads are on the same contig
            # if not change the contig variable and
            # reinitialise the variables
            # pass results to Aggregator class functions

                if current_tid != last_tid:
                    if (self.condition == "range" and tag is not None):
                        E.info("%s:Intervals with %0.2f tag: %s" % (contig_id,
                                                                    threshold,
                                                                    tag))
                        result = Aggregator(samfile=samfile,
                                            results=res)
                        contig_size = contig_dict[last_contig]
                        yield result.tagsInterval(results=res,
                                                  k=self.k,
                                                  contig_size=contig_size,
                                                  threshold=threshold)

                    elif (self.condition == "range"):
                        result = Aggregator(samfile=samfile,
                                            results=res)
                        contig_size = contig_dict[last_contig]
                        yield result.sumInterval(results=res,
                                                 k=self.k,
                                                 contig_size=contig_size)
                    elif self.condition == "cat":
                        yield (contig_id,
                               1,
                               contig_dict[contig_id],
                               running_total)

                    last_tid = current_tid
                    running_total = 0
                    res = []
                else:
                    pass

        # begin counting reads from the list of Counters
        # the running total becomes the read index - currently not used

                    if tocount:
                        running_total += 1

                        if tag is not None:
                            res.append(Result(contig=read.tid,
                                              start=read.pos,
                                              end=read.aend,
                                              idx=running_total,
                                              strand=strand,
                                              tag=read.opt(tag)))
                        else:
                            res.append(Result(contig=read.tid,
                                              start=read.pos,
                                              end=read.aend,
                                              idx=running_total,
                                              strand=strand))
                    else:
                        pass
            else:
                pass


class Aggregator(object):

    ''' A class to perform functions to be called by Profiler.
    Options are: sum over interval and tagged reads over intervals
    '''

    def __init__(self,
                 samfile,
                 results,
                 func=None):
        self.samfile = samfile
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
                        # check the incremented interval
                        # isn't larger than the contig size
                        # reset to contig size if it is
                    elif interval > contig_size:
                        interval = contig_size
                        if interval > contig_size:
                            interval = contig_size
                        else:
                            pass

                            # return total count for each new interval
                            # that contains aligned reads

            if interval > last_interval and intsum > 0:
                yield (self.samfile.getrname(record.contig),
                       last_interval,
                       interval,
                       intsum)
                intsum = 0

    def tagsInterval(self,
                     results,
                     k,
                     contig_size,
                     threshold):
        '''
        Counts the number of reads that match a tag within a region
        and highlights regions > threshold tagged reads
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
            # only yield results if the proportion of reads in
            # an interval exceed a threshold

            if counter > 1 and \
               interval > last_interval and \
               float(counter) / float(reads) >= threshold:
                yield (self.samfile.getrname(record.contig),
                       last_interval,
                       interval,
                       counter)
                reads = 0
                counter = 0


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--condition", dest="condition",
                      type="string", help="either `range`"
                      "or `cat`, `range` is used with --interval.")

    parser.add_option("--interval", dest="interval", type="int",
                      help="the interval size to count reads"
                      "in [default=%default]")

    parser.add_option("--tag", dest="tag", type="string",
                      help="outputs regions with >:term:`threshold`"
                      "reads containing `tag`")

    parser.add_option("--threshold", dest="threshold", type="float",
                      help="a threshold for reads matching :term:`tag`"
                      "in an interval. [default=%default]")

    parser.add_option("--output-filename-pattern", dest="filename",
                      type="string", help="file prefix pattern")

    parser.set_defaults(condition="range",
                        threshold=0.01,
                        interval=10000,
                        tag=None)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if len(args) == 0:
        args.append("-")

    # test presence of bam index file

    bam_index = "%s.bai" % args[0]
    if not os.path.exists(bam_index):
        raise IOError("there is no index file for %s" % args[0])

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

    E.info("reading input file %s" % args[0])
    samfile = pysam.Samfile(args[0], "rb")

    # use the file header information to generate a dictionary of contig ids
    # and sizes
    E.info("Making contig dictionary from input file %s header" % (args[0]))

    contig_it = itertools.izip(samfile.references, samfile.lengths)
    contig_dict = {}

    for i in contig_it:
        contig_dict[i[0]] = i[1]

        # Write out a summary file
        # TODO: have as a user defined option to switch on/off

    if (options.interval and options.tag):
        options.stdout.write('track type=bedGraph name="%s intervals" '
                             'description="Intervals with %s > %0.2f%%" '
                             'visibility=full priority=20\n'
                             % (options.tag,
                                options.tag,
                                round(options.threshold * 100, 2)))

        prof1 = Profiler(samfile=samfile,
                         condition=options.condition,
                         k=options.interval,
                         tag=options.tag,
                         threshold=options.threshold)

        for res in prof1:
            for i in res:
                options.stdout.write("%s\t%i\t%i\t%0.2f\n" % (i[0],
                                                              i[1],
                                                              i[2],
                                                              i[3]))
    elif (options.interval and not options.tag):
        options.stdout.write('track type=bedGraph name="%s'
                             'description=stdout" '
                             'visibility=full priority=20\n'
                             % (options.condition))

        prof1 = Profiler(samfile=samfile,
                         condition=options.condition,
                         k=options.interval,
                         threshold=options.threshold)

        for res in prof1:
            for i in res:
                options.stdout.write("%s\t%i\t%i\t%0.2f\n" % (i[0],
                                                              i[1],
                                                              i[2],
                                                              i[3]))

    elif (options.condition == "cat" and not options.tag):
        options.stdout.write('track type=bedGraph name="%s'
                             'description=stdout" '
                             'visibility=full priority=20\n'
                             % (options.condition))

        prof1 = Profiler(samfile=samfile,
                         condition=options.condition)

        for res in prof1:
            options.stdout.write("%s\t%i\t%i\t%i\n" % (res[0],
                                                       res[1],
                                                       res[2],
                                                       res[3]))

    # write footer and output benchmarkinformation.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
