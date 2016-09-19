'''fastq2summary.py - compute summary stats for a fastq file
======================================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Genomics NGS Sequences FASTQ Annotation

Purpose
-------

This script iterates over a fastq file and outputs
summary statistics for the complete file

The output is a tab-delimited text file with the some of following columns
depending on the option specified:


+----------------+-----------------------------------------------------------+
|*Column*        |*Content*                                                  |
+----------------+-----------------------------------------------------------+
|reads           |total reads in file                                        |
+----------------+-----------------------------------------------------------+
|bases           |total bases in file                                        |
+----------------+-----------------------------------------------------------+
|mean_length     |mean read length                                           |
+----------------+-----------------------------------------------------------+
|median_length   |median read length                                         |
+----------------+-----------------------------------------------------------+
|mean_quality    |mean read quality                                          |
+----------------+-----------------------------------------------------------+
|median_quality  |median read quality                                        |
+----------------+-----------------------------------------------------------+
|nfailed         |number of bases below quality threshold                    |
+----------------+-----------------------------------------------------------+


Usage
-----

Example::

   python fastq2summary.py --guess-format=sanger < in.fastq > out.tsv

In this example we know that our data have quality scores formatted as
sanger. Given that illumina-1.8 quality scores are highly overlapping
with sanger, this option defaults to sanger qualities. In default mode
the script may not be able to distinguish highly overlapping sets of
quality scores.

Type::

   python fastq2summary.py --help

for command line help.

Command line options
--------------------

'''

import sys
import numpy as np
import CGAT.Experiment as E
import CGAT.Fastq as Fastq


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "--guess-format", dest="guess_format", type="choice",
        choices=('sanger', 'solexa', 'phred64',
                 'illumina-1.8', 'integer'),
        help="The default behaviour of the script is to guess \
        the quality format of the input fastq file. The user \
        can specify the quality format of the input file using \
        the --format option. The script will use this format if \
        sequences qualities are ambiguous.[default=%default].")

    parser.add_option(
        "-f", "--target-format", dest="change_format",
        type="choice", choices=('sanger', 'solexa', 'phred64',
                                'illumina-1.8', 'integer'),
        help="The script guesses the quality format of the input \
        file and converts quality scores to the destination \
        format unless --format is specified [default=%default].")

    parser.set_defaults(
        change_format=None,
        guess_format=None,
        min_quality=10)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.change_format:
        iterator = Fastq.iterate_convert(options.stdin,
                                         format=options.change_format,
                                         guess=options.guess_format)
    else:
        iterator = Fastq.iterate_guess(options.stdin,
                                       guess=options.guess_format)

    min_quality = options.min_quality
    number_of_reads = 0
    number_of_bases = 0
    read_lengths = []
    read_qualities = []
    bases_below_min = 0

    for record in iterator:
        number_of_reads += 1
        quals = record.toPhred()
        length_read = len(quals)
        number_of_bases += length_read
        bases_below_min += len([x for x in quals if x < min_quality])
        read_lengths.append(length_read)
        read_qualities.append(np.mean(quals))

    mean_length = round(np.mean(read_lengths), 2)
    median_length = round(np.median(read_lengths), 2)
    mean_quality = round(np.mean(read_qualities), 2)
    median_quality = round(np.median(read_qualities), 2)

    options.stdout.write(
        "reads\tbases\tmean_length\tmedian_length\tmean_quality\tmedian_quality\tnfailed\n")

    options.stdout.write(
        "%i\t%i\t%s\t%s\t%s\t%s\t%i\n" % (number_of_reads, number_of_bases,
                                          str(mean_length),
                                          str(median_length),
                                          str(mean_quality),
                                          str(median_quality),
                                          bases_below_min))
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
