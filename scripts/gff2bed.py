'''gff2bed.py - convert from gff/gtf to bed
=========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics Intervals GFF BED Conversion

Purpose
--------

This script converts GFF or GTF formatted files to BED formatted
files.

Documentation
--------------

Users can select the field from the GTF file to be used in the name
field of the BED file using ``--set-name``. Choices include "gene_id",
"transcript_id", "class", "family", "feature", "source" and "repName".
To specify the input is in GTF format use --is-gtf.

BED files can contain multiple tracks. If required, users can use the
"feature" or "source" fields in the input GFF file to specifiy
different tracks in the BED file (default none).

Usage
------

Example::

   # View input GTF file
   head tests/gff2bed.py/mm9_ens67_geneset_100.gtf

   # Convert GTF to bed format using gene_id as name and group by GTF feature
   cat tests/gff2bed.py/mm9_ens67_geneset_100.gtf | cgat gff2bed.py --is-gtf --set-name=gene_id --track=feature > mm9_ens67_geneset_100_feature.bed

+-------------------------------------------------------+
|track name=CDS                                         |
+------+---------+---------+--------------------+---+---+
|chr18 |3122494  |3123412  |ENSMUSG00000091539  |0  |-  |
+------+---------+---------+--------------------+---+---+
|chr18 |3327491  |3327535  |ENSMUSG00000063889  |0  |-  |
+------+---------+---------+--------------------+---+---+
|chr18 |3325358  |3325476  |ENSMUSG00000063889  |0  |-  |
+------+---------+---------+--------------------+---+---+

Command line options
--------------------

'''

import sys
import itertools
import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.Bed as Bed


def main(argv=sys.argv):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--is-gtf", dest="is_gtf", action="store_true",
                      help="input file is in gtf format [default=%default] ")

    parser.add_option(
        "--set-name", dest="name", type="choice",
        help="field from the GFF/GTF file to use as the "
        "name field in the BED file [%default]",
        choices=("gene_id", "transcript_id", "class", "family",
                 "feature", "source", "repName"))

    parser.add_option(
        "--track", dest="track", type="choice",
        choices=("feature", "source", None),
        help="use feature/source field to define BED tracks "
        "[default=%default]")

    parser.set_defaults(
        track=None,
        name="gene_id",
        is_gtf=False)

    (options, args) = E.Start(parser, add_pipe_options=True)

    ninput, noutput = 0, 0

    is_gtf, name = options.is_gtf, options.name
    iterator = GTF.iterator(options.stdin)

    if options.track:
        all_input = list(iterator)

        if options.track == "feature":
            grouper = lambda x: x.feature
        elif options.track == "source":
            grouper = lambda x: x.source

        all_input.sort(key=grouper)

        bed = Bed.Bed()
        for key, vals in itertools.groupby(all_input, grouper):
            options.stdout.write("track name=%s\n" % key)
            for gff in vals:
                ninput += 1
                bed.fromGTF(gff, is_gtf=is_gtf, name=name)
                options.stdout.write(str(bed) + "\n")
                noutput += 1

    else:
        bed = Bed.Bed()
        for gff in iterator:
            ninput += 1
            bed.fromGTF(gff, is_gtf=is_gtf, name=name)
            options.stdout.write(str(bed) + "\n")

            noutput += 1

    E.info("ninput=%i, noutput=%i" % (ninput, noutput))
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
