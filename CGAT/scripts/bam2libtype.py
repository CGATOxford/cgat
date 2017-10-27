"""bam2libtype.py - determine the library type of a bam file
============================================================

Author: Adam Cribbs

Purpose
-------

This tool determines the library type of a BAM file. The naming
convention used is from the salmon documentation:
http://salmon.readthedocs.io/en/latest/library_type.html.

BAM files need to have a corresponding index file i.e. example.bam
and example.bam.bai


Usage
-----

    cat example.bam | cgat bam2libtype > out.tsv

options
-------

There are no options for this script, just pass the script a bam file
as the stdin and an outfile as the stdout.


Type::

   python bam2bed.py --help

for command line help.

Command line options
--------------------
"""

import sys
import pysam
import CGATCore.Experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(
        version="%prog version: $Id$", usage=globals()["__doc__"])

    (options, args) = E.Start(parser, argv=argv)

    samfile = pysam.AlignmentFile(options.stdin, "rb")
    outfile = options.stdout

    # initialise counts for each library type
    MSR = 0
    MSF = 0
    ISF = 0
    ISR = 0
    OSF = 0
    OSR = 0
    SR = 0
    SF = 0

    for read in samfile:

        # to handle paired end reads:
        if read.is_paired and read.is_proper_pair:

            # get attributes of read
            read_start = read.reference_start
            read_end = read.reference_end
            read_neg = read.is_reverse

            # specify which read is R1 and which is R2:
            # specify which read is R1 and which is R2:
            if read.is_read1 is True:
                R1_is_reverse = read.is_reverse
                R1_reference_start = read.reference_start

                R2_is_reverse = read.mate_is_reverse
                R2_reference_start = read.next_reference_start
            else:
                R1_is_reverse = read.mate_is_reverse
                R1_reference_start = read.next_reference_start

                R2_is_reverse = read.is_reverse
                R2_reference_start = read.reference_start

                # Decision tree to specify strandness:
                # potential to convert this to a machine learning
                # decision tree algorithm in the future:
            if R1_is_reverse is True:

                if R2_is_reverse is True:

                    MSF += 1
                else:
                    if R2_reference_start - R1_reference_start >= 0:
                        OSR += 1
                    else:
                        ISR += 1

            else:

                if R2_is_reverse is True:

                    if R1_reference_start - R2_reference_start >= 0:

                        OSF += 1
                    else:
                        ISF += 1
                else:
                    MSR += 1
        else:
            if read.is_reverse:
                SR += 1
            else:
                SF += 1

    total = MSR + ISR + OSR + ISF + MSF + OSF + SF + SR

    def total_percent(strand, total):
        return float(strand)/float(total)*100

    MSR_total = total_percent(MSR, total)
    ISR_total = total_percent(ISR, total)
    OSR_total = total_percent(OSR, total)
    ISF_total = total_percent(ISF, total)
    MSF_total = total_percent(MSF, total)
    OSF_total = total_percent(OSF, total)
    SF_total = total_percent(SF, total)
    SR_total = total_percent(SR, total)

    outfile.write("MSR\tISR\tOSR\tISF\tMSF\tOSF\tSF\tSR\n")
    outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                  (int(MSR_total), int(ISR_total), int(OSR_total),
                   int(ISF_total), int(MSF_total),
                   int(OSF_total), int(SF_total), int(SR_total)))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
