'''
bed_vs_bed.py - compute enrichment between two bed files
========================================================

:Author:
:Tags: Python

Purpose
-------

This script is a wrapper around bits_test. It splits
one bedfile by name and computes the overlap of intervals
in the first :term:`bed` file against the other.

Usage
-----

.. Example use case

Example::

   python cgat_script_template.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import glob


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-g", "--filename-genome-sizes", 
        dest="filename_genome_sizes", type="string",
        help="Filename with chromosome sizes.")

    parser.add_option(
        "-i", "--num-iterations",
        dest="num_iterations", type="int",
        help="Number of iterations [%default].")

    parser.set_defaults(
        filename_genome_sizes=None,
        num_iterations=1000)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    tagfile, bedfile = args

    if options.filename_genome_sizes is None:
        raise ValueError("please specify a filename with genome sizes")

    E.info("splitting filename %s into sections" % bedfile)

    # ff = IOTools.FilePool(output_pattern="%s.bed.gz")
    # for bed in Bed.iterator(IOTools.openFile(bedfile)):
    #     ff.write(bed.name, str(bed) + "\n")
    # ff.close()
    ff = glob.glob("*.bed.gz")
    iterations = options.num_iterations
    genomefile = options.filename_genome_sizes
    for testfile in ff:
        E.info("working on %s" % testfile)
        statement = """
        bits_test -a %(testfile)s
                  -b %(tagfile)s
                  -n %(iterations)i
                  -g %(genomefile)s
        """ % locals()

        result = E.run(statement, return_stdout=True)
        print(testfile, result)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
