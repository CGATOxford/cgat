'''bam_vs_bed.py - count context that reads map to
======================================================

:Tags: Genomics NGS Intervals BAM BED Counting

Purpose
-------

This script takes as input a :term:`BAM` file from an RNA-seq or
similar experiment and a :term:`bed` formatted file. The :term:`bed`
formatted file needs at least four columns. The fourth (name) column
is used to group counts.

The script counts the number of alignments overlapping in the first
input file that overlap each feature in the second file. Annotations
in the :term:`bed` file can be overlapping - they are counted
independently.

Note that duplicate intervals will be counted multiple times. This
situation can easily arise when building a set of genomic annotations
based on a geneset with alternative transcripts. For example::

   chr1     10000     20000     protein_coding            # gene1, transrcipt1
   chr1     10000     20000     protein_coding            # gene1, transcript2

Any reads overlapping the interval chr1:10000-20000 will be counted
twice into the protein_coding bin by bedtools. To avoid this, remove any
duplicates from the :term:`bed` file::

   zcat input_with_duplicates.bed.gz | cgat bed2bed --merge-by-name | bgzip > input_without_duplicates.bed.gz

This scripts requires bedtools_ to be installed.

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

import sys
import collections
import itertools
import subprocess
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

    parser.add_option("-m", "--min-overlap", dest="min_overlap",
                      type="float",
                      help="minimum overlap [%default]")

    parser.add_option("-a", "--bam-file", dest="filename_bam",
                      metavar="bam", type="string",
                      help="bam-file to use (required) [%default]")

    parser.add_option("-b", "--bed-file", dest="filename_bed",
                      metavar="bed", type="string",
                      help="bed-file to use (required) [%default]")

    parser.add_option(
        "-s", "--sort-bed", dest="sort_bed",
        action="store_true",
        help="sort the bed file by chromosomal location before "
        "processing. "
        "[%default]")

    parser.add_option(
        "--assume-sorted", dest="sort_bed",
        action="store_false",
        help="assume that the bed-file is sorted by chromosomal location. "
        "[%default]")

    parser.add_option(
        "--split-intervals", dest="split_intervals",
        action="store_true",
        help="treat split BAM intervals, for example spliced intervals, "
        "as separate intervals. Note that a single alignment might be "
        "counted several times as a result. "
        "[%default]")

    parser.set_defaults(
        min_overlap=0.5,
        filename_bam=None,
        filename_bed=None,
        sort_bed=True,
        split_intervals=False,
    )

    # add common options (-h/--help, ...) and parse command line
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

    options.stdout.write("category\talignments\n")

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
        samfile = pysam.AlignmentFile(filename_bam, "rb")
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

    options.stdout.write("total\t%i\n" % total)

    if total == 0:
        E.warn("no data in %s" % filename_bam)
        return

    # SNS: sorting optional, off by default
    if options.sort_bed:
        bedcmd = "<( gunzip < %s | sort -k1,1 -k2,2n)" % filename_bed
    else:
        bedcmd = filename_bed

    if options.split_intervals:
        split = "-split"
    else:
        split = ""

    # IMS: newer versions of intersectBed have a very high memory
    #      requirement unless passed sorted bed files.
    statement = """bedtools intersect %(format)s %(filename_bam)s
    -b %(bedcmd)s
    %(split)s
    -sorted -bed -wo -f %(min_overlap)f""" % locals()

    E.info("starting counting process: %s" % statement)
    proc = E.run(statement,
                 return_popen=True,
                 stdout=subprocess.PIPE)

    E.info("counting")
    counts_per_alignment = collections.defaultdict(int)
    take_columns = len(data._fields)

    def iterate(infile):
        for line in infile:
            if not line.strip():
                continue
            yield data._make(line[:-1].split()[:take_columns])

    for read, overlaps in itertools.groupby(
            iterate(IOTools.force_str(proc.stdout)), key=sort_key):
        annotations = [x.name2 for x in overlaps]
        for anno in annotations:
            counts_per_alignment[anno] += 1

    for key, counts in sorted(counts_per_alignment.items()):
        options.stdout.write("%s\t%i\n" % (key, counts))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
