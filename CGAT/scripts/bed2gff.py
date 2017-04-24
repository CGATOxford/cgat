"""bed2gff.py - convert bed to gff/gtf
===================================

:Tags: Genomics Intervals BED GFF Conversion

Purpose
-------

This script converts a :term:`bed`-formatted file to a :term:`gff` or
:term:`gtf`-formatted file.

It aims to populate the appropriate fields in the :term:`gff` file
with columns in the :term:`bed` file.

If ``--as-gtf`` is set and a name column in the :term:`bed` file is
present, its contents will be set as ``gene_id`` and
``transcript_id``. Otherwise, a numeric ``gene_id`` or
``transcript_id`` will be set according to ``--id-format``.

Usage
-----

Example::

   # Preview input bed file
   zcat tests/bed2gff.py/bed3/bed.gz | head
   # Convert BED to GFF format
   cgat bed2gff.py < tests/bed2gff.py/bed3/bed.gz > test1.gff
   # View converted file (excluding logging information)
   cat test1.gtf | grep -v "#" | head


+------+-----+------+-------+-------+---+---+---+---------------------------------------+
|chr1  |bed  |exon  |501    |1000   |.  |.  |.  |gene_id "None"; transcript_id "None";  |
+------+-----+------+-------+-------+---+---+---+---------------------------------------+
|chr1  |bed  |exon  |15001  |16000  |.  |.  |.  |gene_id "None"; transcript_id "None";  |
+------+-----+------+-------+-------+---+---+---+---------------------------------------+

Example::

   # Convert BED to GTF format
   cgat bed2gff.py --as-gtf < tests/bed2gff.py/bed3/bed.gz > test2.gtf
   # View converted file (excluding logging information)
   cat test2.gtf | grep -v "#" | head

+------+-----+------+-------+-------+---+---+---+-----------------------------------------------+
|chr1  |bed  |exon  |501    |1000   |.  |.  |.  |gene_id "00000001"; transcript_id "00000001";  |
+------+-----+------+-------+-------+---+---+---+-----------------------------------------------+
|chr1  |bed  |exon  |15001  |16000  |.  |.  |.  |gene_id "00000002"; transcript_id "00000002";  |
+------+-----+------+-------+-------+---+---+---+-----------------------------------------------+

Type::

   cgat bed2gff.py --help

for command line help.

Command line options
--------------------

"""
import sys
import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.Bed as Bed


def main(argv=sys.argv):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-a", "--as-gtf", dest="as_gtf", action="store_true",
                      help="output as gtf.")

    parser.add_option(
        "-f", "--id-format", dest="id_format", type="string",
        help="format for numeric identifier if --as-gtf is set and "
        "no name in bed file [%default].")

    parser.set_defaults(as_gtf=False,
                        id_format="%08i",
                        test=None)

    (options, args) = E.Start(parser, add_pipe_options=True)

    as_gtf = options.as_gtf
    id_format = options.id_format

    if as_gtf:
        gff = GTF.Entry()
    else:
        gff = GTF.Entry()

    gff.source = "bed"
    gff.feature = "exon"

    ninput, noutput, nskipped = 0, 0, 0

    id = 0
    for bed in Bed.iterator(options.stdin):

        ninput += 1

        gff.contig = bed.contig
        gff.start = bed.start
        gff.end = bed.end
        if bed.fields and len(bed.fields) >= 3:
            gff.strand = bed.fields[2]
        else:
            gff.strand = "."

        if bed.fields and len(bed.fields) >= 2:
            gff.score = bed.fields[1]

        if as_gtf:
            if bed.fields:
                gff.gene_id = bed.fields[0]
                gff.transcript_id = bed.fields[0]
            else:
                id += 1
                gff.gene_id = id_format % id
                gff.transcript_id = id_format % id
        else:
            if bed.fields:
                gff.source = bed.fields[0]

        options.stdout.write(str(gff) + "\n")

        noutput += 1

    E.info("ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput, nskipped))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
