'''
fastq2fastq.py - manipulate fastq files
=============================================

:Tags: Genomics NGS Sequences FASTQ Manipulation

Purpose
-------

This script performs manipulations on :term:`fastq` formatted
files. For example it can be used to change the quality score format
or sample a subset of reads.

The script predominantly is used for manipulation of single fastq
files. However, for some of its functionality it will take paired data
using the ``--pair-fastq-file`` and ``--output-filename-pattern`` options.
This applies to the ``sample`` and ``sort`` methods.

Usage
-----

Example::
  In this example we randomly sample 50% of reads from paired data provided in
  two :term:`fastq` files.

   head in.fastq.1

   @SRR111956.1 HWUSI-EAS618:7:1:27:1582 length=36
   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   +SRR111956.1 HWUSI-EAS618:7:1:27:1582 length=36
   =@A@9@BAB@;@BABA?=;@@BB<A@9@;@2>@;??
   @SRR111956.2 HWUSI-EAS618:7:1:29:1664 length=36
   CCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCC
   +SRR111956.2 HWUSI-EAS618:7:1:29:1664 length=36
   =B@9@0>A<B=B=AAA?;*(@A>(@<=*9=9@BA>7
   @SRR111956.3 HWUSI-EAS618:7:1:38:878 length=36
   AGTGAGCAGGGAAACAATGTCTGTCTAAGAATTTGA

   head in.fastq.2

   +SRR111956.3 HWUSI-EAS618:7:1:38:878 length=36
   <?@BA?;A=@BA>;@@7###################
   @SRR111956.4 HWUSI-EAS618:7:1:38:1783 length=36
   ATTAGTATTATCCATTTATATAATCAATAAAAATGT
   +SRR111956.4 HWUSI-EAS618:7:1:38:1783 length=36
   ?ABBA2CCBBB2?=BB@C>=AAC@A=CBB#######
   @SRR111956.5 HWUSI-EAS618:7:1:39:1305 length=36
   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   +SRR111956.5 HWUSI-EAS618:7:1:39:1305 length=36
   AA>5;A>*91?=AAA@@BBA<B=?ABA>2>?A<BB@

   command-line::
     cat in.fastq.1 | python fastq2fastq.py
                      --method=sample --sample-size 0.5
                      --pair-fastq-file in.fastq.2
                      --output-filename-pattern out.fastq.2
                      > out.fastq.1

   head out.fastq.1
   @SRR111956.1 HWUSI-EAS618:7:1:27:1582 length=36
   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   +
   =@A@9@BAB@;@BABA?=;@@BB<A@9@;@2>@;??
   @SRR111956.2 HWUSI-EAS618:7:1:29:1664 length=36
   CCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCC
   +
   =B@9@0>A<B=B=AAA?;*(@A>(@<=*9=9@BA>7
   @SRR111956.3 HWUSI-EAS618:7:1:38:878 length=36
   AGTGAGCAGGGAAACAATGTCTGTCTAAGAATTTGA
   +
   <?@BA?;A=@BA>;@@7###################
   @SRR111956.4 HWUSI-EAS618:7:1:38:1783 length=36
   ATTAGTATTATCCATTTATATAATCAATAAAAATGT
   +
   ?ABBA2CCBBB2?=BB@C>=AAC@A=CBB#######

Options
-------

The following methods are implemented (``--method``).

``change-format``

    change the quality format to new format given as
    target-format. Options are ``sanger``,
  ``solexa``, ``phred64``, ``integer`` and ``illumina-1.8``

``sample``

    Sub-sample a fastq file. The size of the sample is set by
    --sample-size

``unique``

    Remove duplicate reads based on read name

``trim3``

    Trim a fixed number of nucleotides from the 3' end of reads.
    (see ``--num-bases``). Note that there are better tools for
   trimming.

``trim5``

    Trim a fixed number of nucleotides from the 5' end of reads.
    (see ``--num-bases``). Note that there are better tools for
   trimming.

``sort``

    Sort the fastq file by read name.

``renumber-reads``

    Rename the reads based on pattern given in ``--pattern-identifier``
    e.g. ``--pattern-identifier="read_%010i"``

Type::

   python fastq2fastq.py --help

for command line help.


Command line options
--------------------

'''

import os
import sys
import re
import random
import CGAT.IOTools as IOTools
import CGAT.Experiment as E
import CGAT.Fastq as Fastq
import CGAT.Genomics as Genomics


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=(
                          "apply",
                          "change-format",
                          "renumber-reads",
                          "sample",
                          "sort",
                          "trim3",
                          "trim5",
                          "unique",
                          "reverse-complement",
                          "grep"),
                      help="method to apply [%default]")

    parser.add_option(
        "--target-format", dest="target_format", type="choice",
        choices=('sanger', 'solexa', 'phred64', 'integer', 'illumina-1.8'),
        help="guess quality score format and set quality scores "
        "to format [default=%default].")

    parser.add_option(
        "--guess-format", dest="guess_format", type="choice",
        choices=('sanger', 'solexa', 'phred64', 'integer', 'illumina-1.8'),
        help="quality score format to assume if ambiguous [default=%default].")

    parser.add_option(
        "--sample-size", dest="sample_size", type="float",
        help="proportion of reads to sample. "
        "Provide a proportion of reads to sample, e.g. 0.1 for 10%, "
        "0.5 for 50%, etc [default=%default].")

    parser.add_option(
        "--pair-fastq-file", dest="pair", type="string",
        help="if data is paired, filename with second pair. "
        "Implemented for sampling [default=%default].")

    parser.add_option(
        "--map-tsv-file", dest="map_tsv_file", type="string",
        help="filename with tab-separated identifiers mapping for "
        "method apply [default=%default].")

    parser.add_option(
        "--num-bases", dest="nbases", type="int",
        help="number of bases to trim [default=%default].")

    parser.add_option(
        "--seed", dest="seed", type="int",
        help="seed for random number generator [default=%default].")

    parser.add_option(
        "--pattern-identifier", dest="renumber_pattern", type="string",
        help="rename reads in file by pattern [default=%default]")

    parser.add_option(
        "--grep-pattern", dest="grep_pattern", type="string",
        help="subset to reads matching pattern [default=%default]")

    parser.set_defaults(
        method=None,
        change_format=None,
        guess_format=None,
        sample_size=0.1,
        nbases=0,
        pair=None,
        apply=None,
        seed=None,
        renumber_pattern="read_%010i",
        grep_pattern=".*")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv, add_output_options=True)

    c = E.Counter()

    if options.method is None:
        raise ValueError("no method specified, please use --method")

    if options.method == "change-format":
        for record in Fastq.iterate_convert(options.stdin,
                                            format=options.target_format,
                                            guess=options.guess_format):
            c.input += 1
            options.stdout.write("%s\n" % record)
            c.output += 1

    elif options.method == "grep":
        for record in Fastq.iterate(options.stdin):
            if re.match(options.grep_pattern, record.seq):
                options.stdout.write("%s\n" % record)

    elif options.method == "reverse-complement":
        for record in Fastq.iterate(options.stdin):
            record.seq = Genomics.complement(record.seq)
            record.quals = record.quals[::-1]
            options.stdout.write("%s\n" % record)

    elif options.method == "sample":
        sample_threshold = min(1.0, options.sample_size)

        random.seed(options.seed)

        if options.pair:
            if not options.output_filename_pattern:
                raise ValueError(
                    "please specify output filename pattern for "
                    "second pair (--output-filename-pattern)")

            outfile1 = options.stdout
            outfile2 = IOTools.openFile(options.output_filename_pattern, "w")

            for record1, record2 in zip(
                    Fastq.iterate(options.stdin),
                    Fastq.iterate(IOTools.openFile(options.pair))):
                c.input += 1
                if random.random() <= sample_threshold:
                    c.output += 1
                    outfile1.write("%s\n" % record1)
                    outfile2.write("%s\n" % record2)
        else:
            for record in Fastq.iterate(options.stdin):
                c.input += 1
                if random.random() <= sample_threshold:
                    c.output += 1
                    options.stdout.write("%s\n" % record)

    elif options.method == "apply":
        ids = set(IOTools.readList(IOTools.openFile(options.apply)))

        for record in Fastq.iterate(options.stdin):
            c.input += 1
            if re.sub(" .*", "", record.identifier).strip() in ids:
                c.output += 1
                options.stdout.write("%s\n" % record)

    elif options.method == "trim3":
        trim3 = options.nbases
        for record in Fastq.iterate(options.stdin):
            c.input += 1
            record.trim(trim3)
            options.stdout.write("%s\n" % record)
            c.output += 1

    elif options.method == "trim5":
        trim5 = options.nbases
        for record in Fastq.iterate(options.stdin):
            c.input += 1
            record.trim5(trim5)
            options.stdout.write("%s\n" % record)
            c.output += 1

    elif options.method == "unique":
        keys = set()
        for record in Fastq.iterate(options.stdin):
            c.input += 1
            if record.identifier in keys:
                continue
            else:
                keys.add(record.identifier)
            options.stdout.write("%s\n" % record)
            c.output += 1

    # Need to change this to incorporate both pairs
    elif options.method == "sort":
        if not options.pair:
            # This is quicker for a single fastq file
            statement = "paste - - - - | sort -k1,1 -t ' ' | tr '\t' '\n'"
            os.system(statement)
        else:
            if not options.output_filename_pattern:
                raise ValueError(
                    "please specify output filename for second pair "
                    "(--output-filename-pattern)")
            E.warn(
                "consider sorting individual fastq files - "
                "this is memory intensive")
            entries1 = {}
            entries2 = {}

            for record1, record2 in zip(
                    Fastq.iterate(options.stdin),
                    Fastq.iterate(IOTools.openFile(options.pair))):
                entries1[
                    record1.identifier[:-2]] = (record1.seq, record1.quals)
                entries2[
                    record2.identifier[:-2]] = (record2.seq, record2.quals)

            outfile1 = options.stdout
            outfile2 = IOTools.openFile(options.output_filename_pattern, "w")
            assert len(set(entries1.keys()).intersection(
                set(entries2.keys()))) == len(entries1),\
                "paired files do not contain the same reads "\
                "need to reconcile files"

            for entry in sorted(entries1):
                outfile1.write("@%s/1\n%s\n+\n%s\n" %
                               (entry, entries1[entry][0], entries1[entry][1]))
                outfile2.write("@%s/2\n%s\n+\n%s\n" %
                               (entry, entries2[entry][0], entries2[entry][1]))

    elif options.method == "renumber-reads":
        id_count = 1
        for record in Fastq.iterate(options.stdin):
            record.identifier = options.renumber_pattern % id_count
            id_count += 1
            options.stdout.write("@%s\n%s\n+\n%s\n" %
                                 (record.identifier, record.seq, record.quals))

    # write footer and output benchmark information.
    E.info("%s" % str(c))
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
