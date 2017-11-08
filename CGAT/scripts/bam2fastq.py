'''bam2fastq.py - output fastq files from a bam-file
=================================================

:Tags: Genomics NGS Sequences BAM FASTQ Conversion

Purpose
-------

This script takes a :term:`bam` formatted file and converts to it to
one or two :term:`fastq` formatted files for single-end or paired-end
data, respectively.

For paired-end data, the first fastq file contains the first read of a
read pair and the other contains the second read of read pair.

Example
-------

For example::

   cat in.bam cgat bam2fastq out.1.fastq.gz out.2.fastq.gz

This command converts the :term:`bam` formatted file in.bam into
:term:`fastq` files containing forward reads (out.1.fastq.gz) and
reverse reads (out.2.fastq.gz).  The output files can alternatively
supplied via the option ``--output-pattern-filename``. The statement
below will create the same two output files::

   cat in.bam cgat bam2fastq --output-filename-pattern=out.%s.fastq.gz

Type::

   python bam2fastq.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import tempfile
import shutil
import CGAT.Experiment as E
import CGAT.IOTools as IOTools

import pysam


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.set_defaults(
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv, add_output_options=True)

    # do sth
    if len(args) == 1:
        fastqfile1 = args[0]
        fastqfile2 = options.output_filename_pattern % "2"
    elif len(args) == 2:
        fastqfile1, fastqfile2 = args
    else:
        fastqfile1 = options.output_filename_pattern % "1"
        fastqfile2 = options.output_filename_pattern % "2"

    # only output compressed data
    if not fastqfile1.endswith(".gz"):
        fastqfile1 += ".gz"
    if not fastqfile2.endswith(".gz"):
        fastqfile2 += ".gz"

    if options.stdin != sys.stdin:
        samfile = pysam.AlignmentFile(options.stdin.name, "rb")
    else:
        samfile = pysam.AlignmentFile("-", "rb")

    tmpdir = tempfile.mkdtemp()

    outtemp1 = os.path.join(tmpdir, "pair1.gz")
    outtemp2 = os.path.join(tmpdir, "pair2.gz")

    outstream1 = IOTools.openFile(outtemp1, "w")
    outstream2 = IOTools.openFile(outtemp2, "w")

    E.info('writing fastq files to temporary directory %s' % tmpdir)

    found1, found2 = set(), set()
    read1_qlen, read2_qlen = 0, 0

    c = E.Counter()
    for read in samfile.fetch(until_eof=True):
        c.input += 1
        if not read.is_paired:
            outstream1.write(
                "\t".join((read.qname, read.seq, read.qual)) + "\n")
            found1.add(read.qname)
            if not read1_qlen:
                read1_qlen = read.qlen
            c.unpaired += 1
        elif read.is_read1:
            outstream1.write(
                "\t".join((read.qname, read.seq, read.qual)) + "\n")
            found1.add(read.qname)
            if not read1_qlen:
                read1_qlen = read.qlen
            c.output1 += 1
        elif read.is_read2:
            if read.qname not in found2:
                outstream2.write(
                    "\t".join((read.qname, read.seq, read.qual)) + "\n")
                found2.add(read.qname)
                if not read2_qlen:
                    read2_qlen = read.qlen
                c.output2 += 1

    if c.unpaired == 0 and c.output1 == 0 and c.output2 == 0:
        E.warn("no reads were found")
        return

    sort_statement = '''gunzip < %s
    | sort -k1,1
    | awk '{printf("@%%s\\n%%s\\n+\\n%%s\\n", $1,$2,$3)}'
    | gzip > %s'''

    if c.output1 == 0 and c.output2 == 0:
        # single end data:
        outstream1.close()
        outstream2.close()
        E.info("sorting fastq files")
        E.run(sort_statement % (outtemp1, fastqfile1))

    else:
        # paired end data
        for qname in found2.difference(found1):
            outstream1.write(
                "\t".join((qname, "N" * read1_qlen, "B" * read1_qlen)) + "\n")
            c.extra1 += 1

        for qname in found1.difference(found2):
            outstream2.write(
                "\t".join((qname, "N" * read2_qlen, "B" * read2_qlen)) + "\n")
            c.extra2 += 1

        E.info("%s" % str(c))

        outstream1.close()
        outstream2.close()

        E.info("sorting fastq files")
        E.run(sort_statement % (outtemp1, fastqfile1))
        E.run(sort_statement % (outtemp2, fastqfile2))

    shutil.rmtree(tmpdir)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
