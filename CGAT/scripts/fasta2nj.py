'''
fasta2nj.py - convert fasta file to nj input
============================================

:Tags: Python

Purpose
-------

convert a fasta file to NJ input.

This script translates identifiers like

species|transcripts|gene|class to transcript_species GENEID=gene

Usage
-----

Example::

   python fasta2nj.py --help

Type::

   python fasta2nj.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re
import CGAT.Experiment as E
import CGAT.IOTools as IOTools


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: fasta2nj.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-m", "--map", dest="filename_map", type="string",
                      help="filename with mapping of species ids to swissprot species ids.")

    parser.set_defaults(
        separator="|",
        filename_map=None,
    )

    (options, args) = E.Start(parser)

    if options.filename_map:
        map_species2sp = IOTools.ReadMap(open(options.filename_map, "r"))

    ninput, noutput, nerrors = 0, 0, 0
    for line in sys.stdin:
        if line[0] == ">":
            ninput += 1

            id = re.match(">([^/ \t]+)", line[:-1]).groups()[0]
            data = id.split(options.separator)

            species = data[0]

            if len(data) == 2:
                gene = data[1]
                transcript = None
            elif len(data) >= 3:
                gene = data[2]
                transcript = data[1]

            if map_species2sp:
                try:
                    species = map_species2sp[species]
                except IndexError:
                    nerrors += 1
                    if options.loglevel >= 1:
                        options.stdlog.write(
                            "# could not map species %s\n" % species)
            if transcript:
                options.stdout.write(
                    ">%s_%s GENEID=%s\n" % (transcript, species, gene))
            else:
                options.stdout.write(">%s_%s\n" % (species, gene))
            noutput += 1
        else:
            options.stdout.write(line)

    if options.loglevel >= 1:
        options.stdlog.write(
            "# ninput=%i, noutput=%i, nerrors=%i\n" % (ninput, noutput, nerrors))
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
