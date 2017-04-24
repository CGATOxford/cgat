"""
bam_vs_bam.py - compute coverage correlation between bam files
==============================================================

:Tags: Genomics NGS BAM Comparison

Purpose
-------

Compare per base coverage between two :term:`bam` formatted files.

Usage
-----

Example::

   python bam_vs_bam.py in1.bam in2.bam

This command generates a tab delimited output with columns chromosome,
base coordinate, number of overlapping reads in in1.bam, and number of
overlapping reads in in2.bam.

Type::

   python bam_vs_bam.py --help

for command line help.

Documentation
-------------

This tools allows users to compare the per base coverage between
two BAM files. The output includes all bases in the supplied reference
fasta except those with no coverage in the input BAMs.

At present the --interval or -i option has not been implemented.

Command line options
--------------------

``--regex-identifier``
    supply a regex to extract an identifier from the filenames.
    defualts to using the filename

"""

import sys
import re
import pysam
import CGAT.Experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option("-i", "--intervals-bed-file", dest="filename_intervals",
                      type="string",
                      help="filename with intervals to use "
                      "[%default].")

    parser.add_option("-e", "--regex-identifier", dest="regex_identifier",
                      type="string",
                      help="regular expression to extract identifier from "
                      "filename [%default].")

    parser.set_defaults(
        filename_intervals=None,
        regex_identifier="(.*)",
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if len(args) < 1:
        raise ValueError("please supply at least two BAM files.")

    samfiles = []
    for f in args:
        samfiles.append(pysam.AlignmentFile(f, "rb"))

    if options.filename_intervals:
        raise NotImplementedError(
            "It is not yet possible to specify intervals of interest.\
            Repeat command without intervals option.")

    titles = [re.search(options.regex_identifier, x).groups()[0] for x in args]

    options.stdout.write("contig\tpos\t%s\n" % "\t".join(titles))

    ninput, nskipped, noutput = 0, 0, 0
    contigs = samfiles[0].references

    for contig in contigs:

        missing_contig = False

        positions = {}

        # lazy way: use dictionary
        for x, f in enumerate(samfiles):
            try:
                i = f.pileup(contig)
            except ValueError:
                missing_contig = True
                break

            for v in i:
                vp = v.pos
                if vp in positions:
                    positions[vp].append(v.n)
                else:
                    positions[vp] = [0] * x + [v.n]

            # fill with 0 those not touched in this file
            for p in list(positions.keys()):
                if len(positions[p]) <= x:
                    positions[p].append(0)

        if missing_contig:
            nskipped += 1
            continue

        noutput += 1
        for pos in sorted(positions.keys()):
            vals = positions[pos]
            options.stdout.write("%s\t%i\t%s\n" % (contig, pos,
                                                   "\t".join(map(str, vals))))

    E.info("ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput, nskipped))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
