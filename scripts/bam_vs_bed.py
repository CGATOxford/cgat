'''
bam_vs_bed.py - count context that reads map to
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics NGS Intervals BAM BED Counting

Purpose
-------

This script takes as input a :term:`BAM` file from an RNASeq experiment
and a :term:`bed` formatted file. The :term:`bed` formatted file needs
at least four columns. The fourth (name) column is used to group counts.

It counts the number of alignments overlapping in the first input
file and that overlap each feature in the second file. Annotations in the
:term:`bed` file can be overlapping - they are counted independently.

This scripts requires bedtools to be installed.

Options
-------

-a, --bam-file / -b, --bed-file
    These are the input files. They can also be provided as provided as
    positional arguements, with the bam file being first and the (gziped
    or uncompressed) bed file coming second

-m, --min-overlap
    Using this option will only count reads if they overlap with a bed entry
    by a certain minimum fraction of the read.

Example
-------

Example::

   python bam_vs_bed.py in.bam in.bed.gz

Usage
-----

Type::

   cgat bam_vs_bed BAM BED [OPTIONS]
   cgat bam_vs_bed --bam-file=BAM --bed-file=BED [OPTIONS]

where BAM is either a bam or bed file and BED is a bed file.

Type::

   cgat bam_vs_bed --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import tempfile
import collections
import itertools

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import pysam
import CGAT.Bed as Bed


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-c", "--counting-mode", dest="counting_mode",
                      action="store_true",
                      help="execute script in counting mode ;) [%default]")

    parser.add_option("-m", "--min-overlap", dest="min_overlap",
                      type="float",
                      help="minimum overlap [%default]")

    parser.add_option("-s", "--sort-bed", dest="sort_bed",
                      action="store_true",
                      help="sort the bed file [%default]")

    parser.add_option("-a", "--bam-file", dest="filename_bam",
                      metavar="bam", type="string",
                      help="bam-file to use (required) [%default]")

    parser.add_option("-b", "--bed-file", dest="filename_bed",
                      metavar="bed", type="string",
                      help="bed-file to use (required) [%default]")

    parser.set_defaults(
        counting_mode=False,
        sort_bed=False,
        min_overlap=0.5,
        keep_temp=False,
        filename_bam=None,
        filename_bed=None,
    )

    # add common options (-h/--help, ...) and parse command line

    if "-c" in argv:
        (options, args) = E.Start(parser, argv=argv, quiet=True)

    else:
        (options, args) = E.Start(parser, argv=argv)

    filename_bam = options.filename_bam
    filename_bed = options.filename_bed

    if filename_bam is None and filename_bed is None:
        if len(args) != 2:
            raise ValueError(
                "please supply a bam and a bed file or two bed-files.")

        filename_bam, filename_bed = args

    if filename_bed is None:
        raise ValueError("please supply a bed file to compare to.")

    if filename_bam is None:
        raise ValueError("please supply a bam file to compare with.")

    E.info("intersecting the two files")

    min_overlap = options.min_overlap

    # get number of columns of reference bed file
    for bed in Bed.iterator(IOTools.openFile(filename_bed)):
        ncolumns_bed = bed.columns
        break

    E.info("assuming %s is bed%i format" % (filename_bed, ncolumns_bed))

    if ncolumns_bed < 4:
        raise ValueError("please supply a name attribute in the bed file")

    # get information about
    if filename_bam.endswith(".bam"):
        format = "-abam"
        samfile = pysam.Samfile(filename_bam, "rb")
        total = samfile.mapped
        # latest bedtools uses bed12 format when bam is input
        ncolumns_bam = 12
        # count per read
        sort_key = lambda x: x.name
    else:
        format = "-a"
        total = IOTools.getNumLines(filename_bam)
        # get bed format
        ncolumns_bam = 0
        for bed in Bed.iterator(IOTools.openFile(filename_bam)):
            ncolumns_bam = bed.columns
            break

        if ncolumns_bam > 0:
            E.info("assuming %s is bed%i fomat" % (filename_bam, ncolumns_bam))
            if ncolumns_bam == 3:
                # count per interval
                sort_key = lambda x: (x.contig, x.start, x.end)
            else:
                # count per interval category
                sort_key = lambda x: x.name

    # use fields for bam/bed file (regions to count with)
    data_fields = [
        "contig", "start", "end", "name",
        "score", "strand", "thickstart", "thickend", "rgb",
        "blockcount", "blockstarts", "blockends"][:ncolumns_bam]

    # add fields for second bed (regions to count in)
    data_fields.extend([
        "contig2", "start2", "end2", "name2",
        "score2", "strand2", "thickstart2", "thickend2", "rgb2",
        "blockcount2", "blockstarts2", "blockends2"][:ncolumns_bed])

    # add bases overlap
    data_fields.append("bases_overlap")

    data = collections.namedtuple("data", data_fields)

    if total == 0:
        if not options.counting_mode:
            E.warn("no data in %s" % filename_bam)
        return

    # IMS: newer versions of intersectBed have a very high memory
    #     requirement unless passed sorted bed files.

    if not options.counting_mode:

        # SNS sorting optional, off by default
        if options.sort_bed:
            sort_stat = "| sort -k1,1 -k2,2n"
        else:
            sort_stat = ""

        this_script = os.path.abspath(__file__)
        bam_path = os.path.abspath(filename_bam)
        bed_path = os.path.abspath(filename_bed)

        # note bedtools can handle gzipped files now.
        E.info("counting")
        statement = """intersectBed %(format)s %(filename_bam)s
        -b %(filename_bed)s %(sort_stat)s
        -sorted -bed -wo -f %(min_overlap)f
        | python %(this_script)s -a %(filename_bam)s -b %(filename_bed)s
                 -c True
        """ % locals()

        E.info("running %s" % statement)
        retcode = E.run(statement)

        if retcode != 0:
            raise ValueError("error while executing statement %s" % statement)

    else:

        options.stdout.write("category\talignments\n")
        options.stdout.write("total\t%i\n" % total)

        counts_per_alignment = collections.defaultdict(int)

        take_columns = len(data._fields)

        for line in sys.stdin:
            if not line.strip():
                continue
            entry = data._make(line[:-1].split()[:take_columns])
            counts_per_alignment[entry.name2] += 1

        for key, counts in counts_per_alignment.iteritems():
            options.stdout.write("%s\t%i\n" % (key, counts))

    # write footer and output benchmark information.
    if not options.counting_mode:
        E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
