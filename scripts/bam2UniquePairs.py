'''
bam2UniquePairs.py - filter/report uniquely or bes mapped
                     read pairs from a (bwa!) bam-file
===============================================================

:Author: Steve Sansom
:Release: $Id$
:Date: |today|
:Tags: Genomics NGS

Purpose
-------

Utitily script to report and/or filter out "uniquely or best" mapped
 properly paired reads

Reports:
1. The percentage of properly mapped read pairs with at least
   one uniquely mapped (XT=U) read

2. The percentage of properly mapped read pairs with at least one best
   mapped (X0-1) read

3. The percentage of properly mapped read pairs with at least one
   uniquely or best mapped (X0-1) read

If outfile is specified, reads are emitted when they are properly
paired and the pair has at least one read that is either best or
uniquely mapped.

Duplication is ignored.

Only BWA and GSNAP are supported.

TODO: cache and emit reads rather than iterating over the samfile twice...

'''

import os
import sys
import re
import optparse
import collections
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import pysam


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id: "
                            "cgat_script_template.py "
                            "2871 2010-03-03 10:20:44Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-f", "--filename", dest="filename", type="string",
                      help="bamfile")

    parser.add_option("-a", "--aligner", dest="aligner", type="string",
                      help="Algorithm used to align the reads - "
                           "supported aligners are: "
                           "bwa and gsnap", default="must_specify")

    parser.add_option("-r", "--output-report", type="string", dest="report",
                      help="file name of report", default="")

    parser.add_option("-n", "--number-of-reads",
                      type="int",
                      dest="number_of_reads",
                      help="Number of reads to process",
                      default=-1)

    parser.add_option("-o", "--outfile", dest="outfile", type="string",
                      help="Optional: specify an outfile to write"
                           "the uniquely mapping pairs to",
                      default="")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv, add_output_options=True)

    # Check the aligner is specified supported
    supported_aligners = ["bwa", "gsnap"]
    if options.aligner == "must_specify":
        raise ValueError("please specify a supported aligner. "
                         "supported aligners: "
                         "%s" ", ".join(supported_aligners))

    if options.aligner not in supported_aligners:
        raise ValueError(
            "%s is not a supported as aligner" % options.aligner)

    E.info("Aligner specified is %s" % options.aligner)

    # Check that either a report or outfile name has been specified
    if options.report == "" and options.outfile == "":
        raise ValueError("A filename must be specified"
                         "for a report or alignments")

    # Analyse the bamfile
    samfile = pysam.Samfile(options.filename, "rb")
    uniq_map, best_map, uORb_map = {}, {}, {}
    properly_paired = 0

    if options.number_of_reads != -1:
        sample = True
        E.info("Sampling first %i reads" % options.number_of_reads)
    else:
        sample = False
        E.info("Processing all reads")

    nreads = 0
    for read in samfile.fetch():

        if sample:
            if nreads > options.number_of_reads:
                break
            else:
                nreads += 1

        if read.is_proper_pair:
            tagd = dict(read.tags)
            u, b, key = False, False, read.qname

            if options.aligner == "bwa":
                if tagd["XT"] == "U":
                    u = True
                    uniq_map[key] = 1

                if "X0" in tagd:
                    if tagd["X0"] == 1:
                        b = True
                        best_map[key] = 1

            elif options.aligner == "gsnap":
                if tagd["X2"] == 0:
                    u = True
                    uniq_map[key] = 1
                if tagd["XQ"] > tagd["X2"]:
                    b = True
                    best_map[key] = 1

            else:
                raise ValueError("Aligner not recognised")

            if u is True or b is True:
                uORb_map[key] = 1

            properly_paired += 1

    samfile.close()

    npp = properly_paired / 2

    E.info("No proper pairs considered: %s" % npp)
    E.info("No pairs uniquely or best placed: %s" % len(uORb_map.keys()))

    # Write a tabular report if report name given
    if options.report != "":

        E.info("Writing report on no. proper pairs with unique/best reads")

        def _row(x, npp=npp):
            name, d = x
            n = len(d.keys())
            pc = float(n) / npp * 100
            line = "%s\t%i\t%.2f" % (name, n, pc)
            return(line)

        header = "\t".join(["pair_criteria", "n_proper_pairs",
                            "percent_proper_pairs"])

        with open(options.report, "w") as report:
            report.write(header + "\n")
            for x in [("unique", uniq_map), ("best", best_map),
                      ("unique_or_best", uORb_map)]:
                report.write(_row(x) + "\n")

    # Create new bam containing uniquely mapping read pairs
    # if outfile specified
    if options.outfile != "":

        E.info("Writing proper pairs with unique or best read to %s" %
               options.outfile)

        samfile = pysam.Samfile(options.filename, "rb")
        outbam = pysam.Samfile(options.outfile, "wb", template=samfile)
        
        for read in samfile.fetch():

            if read.is_proper_pair:
                if read.qname in uORb_map:
                    outbam.write(read)
        samfile.close()
        outbam.close()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
