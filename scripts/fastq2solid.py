'''
fastq2solid.py - convert fastq to solid
=======================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Converts a fastq-formatted file to two solid
output files (.csfasta and .qual)

Usage
-----

Example::

   python cgat_script_template.py --help

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Fastq as Fastq


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(
        version="%prog version: $Id",
        usage=globals()["__doc__"])

    parser.add_option(
        "-f", "--target-format", dest="change_format", type="choice",
        choices=('sanger', 'solexa', 'phred64', 'integer'),
        help="set quality scores to format "
        "[default=%default].")

    parser.add_option(
        "--guess-format", dest="guess_format", type="choice",
        choices=('sanger', 'solexa', 'phred64', 'integer'),
        help="quality score format to assume if ambiguous [default=%default].")

    parser.add_option("--pattern-identifier", dest="pattern", type="string",
                      help="filename prefix [default=%default].")

    parser.set_defaults(
        change_format=None,
        guess_format=None,
        pattern="%s.gz"
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    c = E.Counter()

    outfile_seq = IOTools.openFile(options.pattern % "csfasta", "w")
    outfile_qual = IOTools.openFile(options.pattern % "qual", "w")

    if options.change_format:
        iter = Fastq.iterate_convert(options.stdin,
                                     format=options.change_format,
                                     guess=options.guess_format)
    else:
        iter = Fastq.iterate(options.stdin)

    for record in iter:
        c.input += 1
        outfile_seq.write(">%s\n%s\n" % (record.identifier, record.seq))
        outfile_qual.write(">%s\n%s\n" % (record.identifier, record.quals))
        c.output += 1

    outfile_seq.close()
    outfile_qual.close()

    # write footer and output benchmark information.
    E.info("%s" % str(c))
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
