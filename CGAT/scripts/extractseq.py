'''
extractseq.py - extract sequences/sequence regions from a fasta file
====================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python extractseq.py --help

Type::

   python extractseq.py --help

for command line help.

Command line options
--------------------

'''
import sys
import CGAT.Experiment as E
import CGAT.IndexedFasta as IndexedFasta


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: extractseq.py 2861 2010-02-23 17:36:32Z andreas $")

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="pattern to look for sequence filename.")

    parser.add_option("-d", "--identifier", dest="identifier", type="string",
                      help="identifier(s).")

    parser.add_option("-o", "--output-coordinate-format", dest="output_coordinate_format", type="choice",
                      choices=("full", "long"),
                      help="""output format of coordinates. Output format is contig:strand:from:to in zero based
/forward/reverse strand coordinates in open/closed notation. 'long' includes the contig length as fifth field"""  )

    parser.add_option("--input-format", dest="input_format", type="choice",
                      choices=("list", "id-compressed"),
                      help="input format.")

    parser.add_option("-i", "--input-coordinate-format", dest="input_coordinate_format", type="choice",
                      choices=("zero-both", "zero-forward"),
                      help="coordinate format.")

    parser.add_option("-e", "--extend-region", dest="extend_region", type="int",
                      help="regions are extended by this margin at either end.")

    parser.add_option("-r", "--shorten-region", dest="shorten_region", type="int",
                      help="regions are shortened by this margin at either end.")

    parser.set_defaults(
        genome_file=None,
        identifier=None,
        input_coordinate_format="zero-both",
        output_coordinate_format="full",
        input_format="list",
        extend_region=0,
        shorten_region=0,
    )

    (options, args) = E.Start(parser)

    fasta = IndexedFasta.IndexedFasta(options.genome_file)

    lines = []
    if options.identifier:
        lines += [x.split(":") for x in options.identifier.split(",")]

    if args:
        lines += [x.split(":") for x in args]

    if len(lines) == 0:
        lines = [x[
                    :-1].split("\t") for x in ([x for x in options.stdin.readlines() if x[0] != "#"])]

    ninput, nskipped, noutput = 0, 0, 0
    for data in lines:

        if options.input_format == "list":
            if len(data) < 4:
                sbjct_token = data[0]
                sbjct_from, sbjct_to = "0", "0"
                sbjct_strand = "+"
            else:
                sbjct_token, sbjct_strand, sbjct_from, sbjct_to = data[:4]
            id = None
        elif options.input_format == "id-compressed":
            id = data[0]
            sbjct_token, sbjct_strand, sbjct_from, sbjct_to = data[
                1].split(":")

        ninput += 1

        try:
            sbjct_from, sbjct_to = int(sbjct_from), int(sbjct_to)
        except ValueError:
            E.warn("skipping line %s" % data)
            nskipped += 1
            continue

        sbjct_from -= (options.extend_region - options.shorten_region)
        sbjct_from = max(0, sbjct_from)
        lcontig = fasta.getLength(sbjct_token)
        if sbjct_to != 0:
            sbjct_to += (options.extend_region - options.shorten_region)
            sbjct_to = min(sbjct_to, lcontig)
        else:
            sbjct_to = lcontig

        if sbjct_to - sbjct_from <= 0:
            nskipped += 1
            continue

        sequence = fasta.getSequence(sbjct_token, sbjct_strand,
                                     sbjct_from, sbjct_to,
                                     converter=IndexedFasta.getConverter(options.input_coordinate_format))

        if options.output_coordinate_format == "full":
            coordinates = "%s:%s:%i:%i" % (sbjct_token,
                                           sbjct_strand,
                                           sbjct_from,
                                           sbjct_to)

        elif options.output_coordinate_format == "long":
            coordinates = "%s:%s:%i:%i:%i" % (sbjct_token,
                                              sbjct_strand,
                                              sbjct_from,
                                              sbjct_to,
                                              lcontig)

        if id:
            options.stdout.write(">%s %s\n%s\n" % (id, coordinates, sequence))
        else:
            options.stdout.write(">%s\n%s\n" % (coordinates, sequence))

        noutput += 1

    E.info("ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput, nskipped))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
