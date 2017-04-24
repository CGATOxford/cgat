'''
split_genome.py - 
======================================================

:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python split_genome.py --help

Type::

   python split_genome.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re

import CGAT.Experiment as E

USAGE = """python %s [OPTIONS] [genomic_sequence] [ < genomic sequence]

Version: $Id: split_genome.py 2781 2009-09-10 11:33:14Z andreas $

Split a genomic fasta file into smaller segments.
""" % sys.argv[0]


def Print(outfile, fragments, options):
    """print fragments to outfile and return overhang."""
    s = "".join(fragments)
    l = min(len(s), options.chunk_size)
    outfile.write("\n".join(s[x:x + options.width]
                            for x in range(0, l, options.width)))
    outfile.write("\n")
    outfile.flush()
    return s[l:]


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: split_genome.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-c", "--chunk-size", dest="chunk_size",
                      help="size of chunks in nucleotides.", type="int")
    parser.add_option("-o", "--filename-pattern-output", dest="filename_pattern_output",
                      help="filename for output (should contain one '%i').", type="string")

    parser.set_defaults(
        chunk_size=200000,
        filename_pattern_output="%i.fasta",
        width=100,
    )

    (options, args) = E.Start(parser, add_pipe_options=True)

    nchunk = 0
    chunksize = 0
    pos = 0
    fragments = []
    outfile = None

    for line in sys.stdin:

        is_header = line[0] == ">"

        if is_header or chunksize > options.chunk_size:

            if outfile:
                rest = Print(outfile, fragments, options)
                chunksize = len(rest)
                pos -= chunksize
                fragments = [rest]
                outfile.close()

            else:
                fragments = []
                chunksize = 0

            nchunk += 1

            outfile = IOTools.openFile(options.filename_pattern_output % nchunk, "w")

            if is_header:
                description = line[1:-1]
                id = re.split("\s", description)[0]
                pos = 0

            outfile.write(">%s|%i|%i %s\n" % (id, nchunk, pos, description))

            if is_header:
                continue

        s = re.sub("\s", "", line[:-1])
        l = len(s)
        pos += l
        chunksize += l
        fragments.append(s)

    if outfile:
        rest = Print(outfile, fragments, options)
        outfile.close()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
