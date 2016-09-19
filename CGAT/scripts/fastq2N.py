'''fastq2N.py -
=============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::

   This script takes as input a fastq file and converts a position
   (base) to an N call. This is to keep all reads in line should there
   be an overrepresentation of N calls in any samples.

Usage
-----

Example::

   python fastq2N.py --help

Type::

   python fastq2N.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.Fastq as Fastq
import gzip

# define functions to be used in main


def replace(fastqfile, baseToReplace):
    '''replaces the specified base with N'''

    # use gzip as default to open the fastq file
    outf = gzip.open("replaced_" + fastqfile, "w")
    fastq = gzip.open(fastqfile)
    iterator = Fastq.iterate(fastq)
    for record in iterator:
        x = list(record.seq)
        x[int(baseToReplace)] = "N"
        record.seq = "".join(x)
        outf.write("@" + record.identifier + "\n" + record.seq +
                   "\n" + "+" + record.identifier + "\n" + record.quals + "\n")


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-i", "--infile", dest="infile", type="string",
                      help="Input filename")
    parser.add_option("-b", "--base-position", dest="base_position",
                      help="base position in which to replace with an N")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if not options.infile:
        raise ValueError("need to specify input fastq file")

    if not options.base_position:
        raise ValueError("need to specify which base to convert")

    # main function
    replace(options.infile, options.base_position)

    # write footer and output benchmark information.
    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
