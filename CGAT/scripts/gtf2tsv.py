'''gtf2tsv.py - convert gtf file to a tab-separated table
======================================================

:Tags: Genomics Genesets

Purpose
-------

convert a gtf formatted file to tab-separated table. The difference to
a plain :term:`gtf` formatted file is that column headers are added,
which can be useful when importing the gene models into a database.

Note that coordinates are converted to 0-based open/closed notation (all on
the forward strand).

By default, the gene_id and transcript_id are extracted from the
attributes field into separated columns.  If
``-f/--attributes-as-columns`` is set, all fields in the attributes
will be split into separate columns.

The script also implements the reverse operation, converting a tab-separated
table into a :term:`gtf` formatted file.

When using the ``-m, --map`` option, the script will output a table
mapping gene identifiers to transcripts or peptides.

USING GFF3 FILE:
The script also can convert gff3 formatted files to tsv files when
specifiying the option --is-gff3 and --attributes-as-columns. Currently only
the full GFF3 to task is implimented. Further improvements to this script can
be made to only output the attributes only, i.e. --output-only-attributes.




Usage
-----

Example::

   cgat gtf2tsv < in.gtf

+------+--------------------------------+-----------+------+------+-----+------+-----+---------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------+
|contig|source                          |feature    |start |end   |score|strand|frame|gene_id        |transcript_id  |attributes                                                                                                                           |
+------+--------------------------------+-----------+------+------+-----+------+-----+---------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------+
|chr19 |processed_transcript            |exon       |66345 |66509 |.    |-     |.    |ENSG00000225373|ENST00000592209|exon_number "1"; gene_name "AC008993.5"; gene_biotype "pseudogene"; transcript_name "AC008993.5-002"; exon_id "ENSE00001701708"      |
+------+--------------------------------+-----------+------+------+-----+------+-----+---------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------+
|chr19 |processed_transcript            |exon       |60520 |60747 |.    |-     |.    |ENSG00000225373|ENST00000592209|exon_number "2"; gene_name "AC008993.5"; gene_biotype "pseudogene"; transcript_name "AC008993.5-002"; exon_id "ENSE00002735807"      |
+------+--------------------------------+-----------+------+------+-----+------+-----+---------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------+
|chr19 |processed_transcript            |exon       |60104 |60162 |.    |-     |.    |ENSG00000225373|ENST00000592209|exon_number "3"; gene_name "AC008993.5"; gene_biotype "pseudogene"; transcript_name "AC008993.5-002"; exon_id "ENSE00002846866"      |
+------+--------------------------------+-----------+------+------+-----+------+-----+---------------+---------------+-------------------------------------------------------------------------------------------------------------------------------------+

To build a map between gene and transcrip identiers, type::

   cgat gtf2tsv --output-map=transcript2gene < in.gtf

+---------------+---------------+
|transcript_id  |gene_id        |
+---------------+---------------+
|ENST00000269812|ENSG00000141934|
+---------------+---------------+
|ENST00000318050|ENSG00000176695|
+---------------+---------------+
|ENST00000327790|ENSG00000141934|
+---------------+---------------+

To run the script to convert a gff3 formatted file to tsv, type::

   cat file.gff3.gz | cgat gtf3tsv --is-gff3 --attributes-as-columns
   > outfile.tsv

Type::

   cgat gtf2tsv --help

for command line help.

Command line options
---------------------

'''
import sys
import re
import CGAT.GTF as GTF
import CGAT.GFF3 as GFF3
import CGAT.Experiment as E


def main(argv=None):
    '''
    main function
    '''

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option(
        "-o", "--output-only-attributes", dest="only_attributes",
        action="store_true",
        help="output only attributes as separate columns "
        "[default=%default].")

    parser.add_option(
        "-f", "--attributes-as-columns", dest="output_full",
        action="store_true",
        help="output attributes as separate columns "
        "[default=%default].")

    parser.add_option("--is-gff3", dest="is_gtf", action="store_false",
                      help="input file is in gtf format [default=%default] ")

    parser.add_option(
        "-i", "--invert", dest="invert", action="store_true",
        help="convert tab-separated table back to gtf "
        "[default=%default].")

    parser.add_option(
        "-m", "--output-map", dest="output_map", type="choice",
        choices=(
            "transcript2gene",
            "peptide2gene",
            "peptide2transcript"),
        help="output a map mapping transcripts to genes "
        "[default=%default].")

    parser.set_defaults(
        only_attributes=False,
        output_full=False,
        invert=False,
        output_map=None,
        is_gtf=True
    )

    (options, args) = E.Start(parser, argv=argv)

    if options.output_full:
        # output full table with column for each attribute

        attributes = set()
        data = []
        if options.is_gtf:
            for gtf in GTF.iterator(options.stdin):
                data.append(gtf)
                attributes = attributes.union(set(gtf.keys()))

        else:
            for gff in GFF3.iterator_from_gff(options.stdin):
                data.append(gff)
                attributes = attributes.union(set(gff.attributes))

        # remove gene_id and transcript_id, as they are used
        # explicitely later
        attributes.difference_update(["gene_id", "transcript_id"])

        attributes = sorted(list(attributes))

        # Select whether gtf of gff for output columns
        if options.is_gtf:
            if options.only_attributes:
                header = ["gene_id", "transcript_id"] + attributes
            else:
                header = ["contig", "source", "feature",
                          "start", "end", "score", "strand",
                          "frame", "gene_id",
                          "transcript_id", ] + attributes
        else:
            if options.only_attributes:
                header = attributes
            else:
                header = ["contig", "source", "feature",
                          "start", "end", "score", "strand",
                          "frame"] + attributes

        attributes_new = header

        options.stdout.write("\t".join(header) + "\n")

        if options.is_gtf:
            for gtf in data:
                first = True
                for a in attributes_new:
                    try:
                        val = getattr(gtf, a)
                    except (AttributeError, KeyError):
                        val = ""
                    if first:
                        options.stdout.write("%s" % val)
                        first = False
                    else:
                        options.stdout.write("\t%s" % val)
                options.stdout.write("\n")
        else:
            for gff in data:
                options.stdout.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t") % (gff.contig,
                                                                             gff.source, gff.feature, gff.start, gff.end,
                                                                             gff.score, gff.strand, gff.frame))

                first = True
                for a in attributes:
                    try:
                        val = (gff.attributes[a])
                    except (AttributeError, KeyError):
                        val = ''
                    if first:
                        options.stdout.write("%s" % val)
                        first = False
                    else:
                        options.stdout.write("\t%s" % val)
                options.stdout.write("\n")

    elif options.invert:

        gtf = GTF.Entry()
        header = None
        for line in options.stdin:
            if line.startswith("#"):
                continue
            data = line[:-1].split("\t")
            if not header:
                header = data
                map_header2column = dict(
                    [(y, x) for x, y in enumerate(header)])
                continue

            # fill gtf entry with data
            try:
                gtf.contig = data[map_header2column["contig"]]
                gtf.source = data[map_header2column["source"]]
                gtf.feature = data[map_header2column["feature"]]
                # subtract -1 to start for 0-based coordinates
                gtf.start = int(data[map_header2column["start"]])
                gtf.end = int(data[map_header2column["end"]])
                gtf.score = data[map_header2column["score"]]
                gtf.strand = data[map_header2column["strand"]]
                gtf.frame = data[map_header2column["frame"]]
                gtf.gene_id = data[map_header2column["gene_id"]]
                gtf.transcript_id = data[map_header2column["transcript_id"]]
                gtf.parseInfo(data[map_header2column["attributes"]], line)
            except KeyError as msg:
                raise KeyError("incomplete entry %s: %s: %s" %
                               (str(data), str(map_header2column), msg))
            if gtf.frame is None:
                gtf.frame = "."
            # output gtf entry in gtf format
            options.stdout.write("%s\n" % str(gtf))

    elif options.output_map:

        if options.output_map == "transcript2gene":
            fr = lambda x: x.transcript_id
            to = lambda x: x.gene_id
            options.stdout.write("transcript_id\tgene_id\n")
        elif options.output_map == "peptide2gene":
            fr = lambda x: x.protein_id
            to = lambda x: x.gene_id
            options.stdout.write("peptide_id\tgene_id\n")
        elif options.output_map == "peptide2transcript":
            fr = lambda x: x.protein_id
            to = lambda x: x.transcript_id
            options.stdout.write("peptide_id\ttranscript_id\n")

        map_fr2to = {}
        for gtf in GTF.iterator(options.stdin):
            try:
                map_fr2to[fr(gtf)] = to(gtf)
            except (AttributeError, KeyError):
                pass

        for x, y in sorted(map_fr2to.items()):
            options.stdout.write("%s\t%s\n" % (x, y))
    else:
        header = ("contig", "source", "feature", "start", "end", "score",
                  "strand", "frame", "gene_id", "transcript_id", "attributes")
        options.stdout.write("\t".join(header) + "\n")

        for gtf in GTF.iterator(options.stdin):
            attributes = []
            for a in list(gtf.keys()):
                if a in ("gene_id", "transcript_id"):
                    continue
                attributes.append('%s %s' % (a, GTF.quote(gtf[a])))

            attributes = "; ".join(attributes)

            # Capture if None and set to . format
            if gtf.frame is None:
                gtf.frame = "."

            options.stdout.write(str(gtf) + "\n")

    E.Stop()

if __name__ == '__main__':
    sys.exit(main())
