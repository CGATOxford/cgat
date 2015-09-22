'''bam2UniquePairs.py - filter/report uniquely mapped read pairs from a (bwa!) bam-file
======================================================================================

:Author: Steve Sansom
:Release: $Id$
:Date: |today|
:Tags: Genomics NGS

Purpose
-------

Utility script to report and/or filter out "uniquely mapped" properly
paired reads

Reports:

1. The percentage of properly mapped read pairs with at least one
   uniquely mapped (XT=U) read

2. The percentage of properly mapped read pairs with at least one best
   mapped (X0-1) read

3. The percentage of properly mapped read pairs with at least one
   uniquely or best mapped (X0-1) read

If outfile is specified, reads are emitted when they are properly
paired and the pair has at least one read that is either best or
uniquely mapped.

Duplication is ignored.

Only BWA is supported.

TODO: cache and emit reads rather than iterating over the samfile twice...

'''

import sys
import CGAT.Experiment as E
import pysam


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-f", "--filename", dest="filename", type="string",
                      help="bamfile")

    parser.add_option("-a", "--aligner", dest="aligner", type="string",
                      help="bamfile", default="bwa")

    parser.add_option("-r", "--output-report", type="string", dest="report",
                      help="bamfile", default="")

    parser.add_option("-o", "--outfile", dest="outfile", type="string",
                      help="bamfile", default="")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv, add_output_options=True)

    # Check the aligner is supported
    if options.aligner != "bwa":
        raise ValueError(
            "Currently only bwa is supported as aligner specific flags are used")

    # Check that either a report or outfile name has been specified
    if options.report == "" and options.outfile == "":
        raise ValueError("Nothing to do")

    # Analyse the bamfile
    samfile = pysam.Samfile(options.filename, "rb")
    uniq_map, best_map, uORb_map = {}, {}, {}
    properly_paired = 0

    for read in samfile.fetch():

        if read.is_proper_pair:
            tagd = dict(read.tags)
            u, b, key = False, False, read.qname

            if tagd["XT"] == "U":
                u = True
                uniq_map[key] = 1

            if "X0" in tagd:
                if tagd["X0"] == 1:
                    b = True
                    best_map[key] = 1

            if u is True or b is True:
                uORb_map[key] = 1

            properly_paired += 1

    samfile.close()

    npp = properly_paired / 2

    E.info("No proper pairs: %s" % npp)

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
