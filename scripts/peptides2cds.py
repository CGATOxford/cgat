'''peptides2cds.py - map peptide sequences onto nucleic acid sequences
===================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

build a map between peptides to cdnas. This script deals with
frameshifts but not with larger indels.

This script requires the :mod:`alignlib_lite.py_ library and :file:`muscle`
to be installed.

Usage
-----

Example::

   python peptides2cds.py --help

Type::

   python peptides2cds.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re
import CGAT.Experiment as E
import alignlib_lite
import CGAT.Genomics as Genomics
import CGAT.FastaIterator as FastaIterator
import CGAT.IOTools as IOTools
import CGAT.Peptides2Cds as Peptides2Cds


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id$")

    parser.add_option("-p", "--peptides-fasta-file", dest="filename_peptides", type="string",
                      help="filename with peptide sequences [%default].")

    parser.add_option("-c", "--cds-gtf-file", "--cdnas", dest="filename_cdna", type="string",
                      help="filename with cdna sequences [%default].")

    parser.add_option("-m", "--map", dest="filename_map", type="string",
                      help="filename with map of peptide identifiers to cdna identifiers [%default].")

    parser.add_option("--output-identifier", dest="output_identifier", type="choice",
                      choices=("cdna", "peptide"),
                      help="output identifier to use [%default].")

    parser.add_option("-f", "--output-format=", dest="output_format", type="choice",
                      choices=("alignment", "fasta"),
                      help="output format.")

    parser.set_defaults(
        peptides=None,
        filename_cdna=None,
        output_format="alignment",
        filename_map=None,
        stop_codons=("TAG", "TAA", "TGA"),
        output_identifier="peptide",
    )

    (options, args) = E.Start(parser, add_pipe_options=True)

    if not options.filename_cdna:
        raise ValueError("please supply filename with cds sequences.")

    if options.filename_peptides:
        infile = open(options.filename_peptides, "r")
        E.info("reading from %s" % options.filename_peptides)
    else:
        E.info("reading from stdin")
        infile = sys.stdin

    if options.filename_map:
        E.info("reading map")
        map_peptide2cds = IOTools.readMap(
            IOTools.openFile(options.filename_map, "r"))
        E.info("read map for %i identifiers" % len(map_peptide2cds))
    else:
        map_peptide2cds = {}

    E.info("reading cds sequences")

    cds_sequences = Genomics.ReadPeptideSequences(
        IOTools.openFile(options.filename_cdna, "r"))

    E.info("read %i cds sequences" % len(cds_sequences))

    ninput, noutput = 0, 0
    nskipped, nnosequence = 0, 0

    # iterate over peptide sequences
    iterator = FastaIterator.FastaIterator(infile)

    use_cds_id = options.output_identifier == "cds"

    for cur_record in iterator:

        ninput += 1

        peptide_identifier = re.split("\s+", cur_record.title)[0]
        cds_identifier = map_peptide2cds.get(
            peptide_identifier, peptide_identifier)

        if cds_identifier not in cds_sequences:
            nnosequence += 1
            continue

        p = cur_record.sequence
        c = cds_sequences[cds_identifier]

        E.debug("processing %s: laa=%i (without gaps=%i), lna=%i" %
                (peptide_identifier, len(p), len(re.sub("-", "", p)), len(c)))

        try:
            map_p2c = Peptides2Cds.getMapPeptide2Cds(p, c, options)
        except ValueError:
            nskipped += 1
            continue

        if use_cds_id:
            identifier = cds_identifier
        else:
            identifier = peptide_identifier

        if options.output_format == "alignment":
            options.stdout.write("\t".join(map(
                str,
                (identifier,
                 alignlib_lite.py_AlignmentFormatEmissions(map_p2c),
                 len(cur_record.sequence),
                 len(cds_sequences[identifier])))) + "\n")

        elif options.output_format == "fasta":

            map_p2c.switchRowCol()

            alignatum = alignlib_lite.py_makeAlignatum(c)

            alignatum.mapOnAlignment(map_p2c, len(p) * 3)

            s = alignatum.getString()
            if len(s) != len(p) * 3:
                raise ValueError(
                    "incomplete aligned string for %s: %s, cds=%s" %
                    (cur_record.title, s, c))

            options.stdout.write(">%s\n%s\n" % (identifier, s))

        noutput += 1
        sys.stdout.flush()

    E.info("ninput=%i, noutput=%i, nnosequence=%i, nskipped=%i" %
           (ninput, noutput, nnosequence, nskipped))

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
