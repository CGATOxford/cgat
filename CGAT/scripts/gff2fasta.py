'''
gff2fasta.py - output sequences from genomic features
=====================================================

:Tags: Genomics Intervals Sequences GFF Fasta Transformation

Purpose
-------

This script outputs the genomic sequences for intervals within
a :term:`gff` or :term: `gtf` formatted file.

The ouput can be optionally masked and filtered.

Usage
-----

If you want to convert a ``features.gff`` file with intervals information
into a :term:`fasta` file containing the sequence of each interval, use this
script as follows::

   python gff2fasta.py --genome-file=hg19 < features.gff > features.fasta

The input can also be a :term:`gtf` formatted file. In that case, use the
``--is-gtf`` option::

   python gff2fasta.py --genome-file=hg19 --is-gtf < features.gtf >\
 features.fasta

If you want to add a polyA tail onto each transcript you can use the `extend`
options:

   python gff2fasta.py --genome-file=hg19 --is-gtf
   --extend-at=3 --extend-by=125 --extend-with=A
   < features.gtf > features.fasta

If you want to merge the sequence of similar features together, please use
``--merge-overlapping``::

   python gff2fasta.py --genome-file=hg19 --merge-overlapping < features.gff >\
 features.fasta

It is possible to filter the output by selecting a minimum or maximum number
of nucleotides in the resultant fasta sequence with ``--max-length`` or
``--min-interval-length`` respectively::

   python gff2fasta.py --genome-file=hg19 --max-length=100\
 < features.gff > features.fasta

Or you can also filter the output by features name with the ``--feature``
option::

   python gff2fasta.py --genome-file=hg19 --feature=exon < features.gff\
 > features.fasta

On the other hand, low-complexity regions can be masked with the ``--masker``
option and a given :term:`gff` formatted file::

   python gff2fasta.py --genome-file=hg19 --masker=dust\
 --maskregions-bed-file=intervals.gff < features.gff > features.fasta

where ``--masker`` can take the following values: ``dust``, ``dustmasker``,
and ``softmask``.

Options
-------

``--is-gtf``
  Tells the script to expect a :term:`gtf` format file

``--genome-file``
  PATH to Fasta file of genome build to use

``--merge-overlapping``
  Merge features in :term:`gtf`/:term:`gff` file that are adjacent and share
  attributes

``--method=filter --filter-method``
  Filter on a :term:`gff` feature such as ``exon`` or ``CDS``

``--maskregions-bed-file``
  Mask sequences in intervals in :term:`gff` file

``--remove-masked-regions``
  Remove sequences in intervals in :term:`gff` file rather than masking them

``--min-interval-length``
  Minimum output sequence length

``--max-length``
  Maximum output sequence length

``--extend-at``
  Extend sequence at 3', 5' or both end.  Optionally '3only' or '5only' will
  return only the 3' or 5' extended sequence

``--extend-by``
  Used in conjunction with ``--extend-at``, the number of nucleotides to extend
  by

``--extend-with``
  Optional. Used in conjunction with ``--extend-at`` and ``--extend-by``.
  Instead of extending by the genomic sequence, extend by this string repeated
  n times, where n is --entend-by


``--masker``
  Masker type to use: dust, dustmasker, soft or none

``--fold-at``
  Fold the fasta sequence every n bases

``--naming-attribute``
  Use this attribute to name the fasta entries

Command line options
--------------------
'''

import sys
import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.Genomics as Genomics
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Intervals as Intervals
import CGAT.Masker as Masker
import bx.intervals.intersection


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option("--is-gtf", dest="is_gtf", action="store_true",
                      help="input is gtf instead of gff.")

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default].")

    parser.add_option(
        "-m", "--merge-adjacent", dest="merge", action="store_true",
        help="merge adjacent intervals with the same attributes."
        " [default=%default]")

    parser.add_option(
        "-e", "--feature", dest="feature", type="string",
        help="filter by a feature, for example 'exon', 'CDS'."
        " If set to the empty string, all entries are output "
        "[%default].")

    parser.add_option(
        "-f", "--maskregions-bed-file", dest="filename_masks",
        type="string", metavar="gff",
        help="mask sequences with regions given in gff file "
        "[%default].")

    parser.add_option(
        "--remove-masked-regions", dest="remove_masked_regions",
        action="store_true",
        help="remove regions instead of masking [%default].")

    parser.add_option(
        "--min-interval-length", dest="min_length", type="int",
        help="set minimum length for sequences output "
        "[%default]")

    parser.add_option(
        "--max-length", dest="max_length", type="int",
        help="set maximum length for sequences output "
        "[%default]")

    parser.add_option(
        "--extend-at", dest="extend_at", type="choice",
        choices=("none", "3", "5", "both", "3only", "5only"),
        help="extend at no end, 3', 5' or both ends. If "
        "3only or 5only are set, only the added sequence "
        "is returned [default=%default]")

    parser.add_option(
        "--header-attributes", dest="header_attr",
        action="store_true",
        help="add GFF entry attributes to the FASTA record"
        " header section")

    parser.add_option(
        "--extend-by", dest="extend_by", type="int",
        help="extend by # bases [default=%default]")

    parser.add_option(
        "--extend-with", dest="extend_with", type="string",
        help="extend using base [default=%default]")

    parser.add_option(
        "--masker", dest="masker", type="choice",
        choices=("dust", "dustmasker", "softmask", "none"),
        help="apply masker [%default].")

    parser.add_option(
        "--fold-at", dest="fold_at", type="int",
        help="fold sequence every n bases[%default].")

    parser.add_option(
        "--fasta-name-attribute", dest="naming_attribute", type="string",
        help="use attribute to name fasta entry. Currently only compatable"
        " with gff format [%default].")

    parser.set_defaults(
        is_gtf=False,
        genome_file=None,
        merge=False,
        feature=None,
        filename_masks=None,
        remove_masked_regions=False,
        min_length=0,
        max_length=0,
        extend_at=None,
        extend_by=100,
        extend_with=None,
        masker=None,
        fold_at=None,
        naming_attribute=False,
        header_attr=False,
    )

    (options, args) = E.Start(parser)

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta(options.genome_file)
        contigs = fasta.getContigSizes()

    if options.is_gtf:
        iterator = GTF.transcript_iterator(GTF.iterator(options.stdin))
    else:
        gffs = GTF.iterator(options.stdin)
        if options.merge:
            iterator = GTF.joined_iterator(gffs)
        else:
            iterator = GTF.chunk_iterator(gffs)

    masks = None
    if options.filename_masks:
        masks = {}
        with IOTools.openFile(options.filename_masks, "r") as infile:
            e = GTF.readAsIntervals(GTF.iterator(infile))

        # convert intervals to intersectors
        for contig in list(e.keys()):
            intersector = bx.intervals.intersection.Intersecter()
            for start, end in e[contig]:
                intersector.add_interval(bx.intervals.Interval(start, end))
            masks[contig] = intersector

    ninput, noutput, nmasked, nskipped_masked = 0, 0, 0, 0
    nskipped_length = 0
    nskipped_noexons = 0

    feature = options.feature

    # iterator is a list containing groups (lists) of features.
    # Each group of features have in common the same transcript ID, in case of
    # GTF files.
    for ichunk in iterator:

        ninput += 1

        if feature:
            chunk = [x for x in ichunk if x.feature == feature]
        else:
            chunk = ichunk

        if len(chunk) == 0:
            nskipped_noexons += 1
            E.info("no features in entry from "
                   "%s:%i..%i - %s" % (ichunk[0].contig,
                                       ichunk[0].start,
                                       ichunk[0].end,
                                       str(ichunk[0])))
            continue

        contig, strand = chunk[0].contig, chunk[0].strand
        
        if options.is_gtf:
            name = chunk[0].transcript_id
        else:
            if options.naming_attribute:
                attr_dict = {x.split("=")[0]: x.split("=")[1]
                             for x in chunk[0].attributes.split(";")}
                name = attr_dict[options.naming_attribute]
            else:
                name = str(chunk[0].attributes)

        lcontig = contigs[contig]
        positive = Genomics.IsPositiveStrand(strand)
        intervals = [(x.start, x.end) for x in chunk]
        intervals.sort()

        if masks:
            if contig in masks:
                masked_regions = []
                for start, end in intervals:
                    masked_regions += [(x.start, x.end)
                                       for x in masks[contig].find(start, end)]

                masked_regions = Intervals.combine(masked_regions)
                if len(masked_regions):
                    nmasked += 1

                if options.remove_masked_regions:
                    intervals = Intervals.truncate(intervals, masked_regions)
                else:
                    raise NotImplementedError("unimplemented")

                if len(intervals) == 0:
                    nskipped_masked += 1
                    if options.loglevel >= 1:
                        options.stdlog.write("# skipped because fully masked: "
                                             "%s: regions=%s masks=%s\n" %
                                             (name,
                                              str([(x.start,
                                                    x.end) for x in chunk]),
                                              masked_regions))
                    continue

        out = intervals

        if options.extend_at and not options.extend_with:
            if options.extend_at == "5only":
                intervals = [(max(0, intervals[0][0] - options.extend_by),
                              intervals[0][0])]
            elif options.extend_at == "3only":
                intervals = [(intervals[-1][1],
                              min(lcontig,
                                  intervals[-1][1] + options.extend_by))]
            else:
                if options.extend_at in ("5", "both"):
                    intervals[0] = (max(0,
                                        intervals[0][0] - options.extend_by),
                                    intervals[0][1])
                if options.extend_at in ("3", "both"):
                    intervals[-1] = (intervals[-1][0],
                                     min(lcontig,
                                         intervals[-1][1] + options.extend_by))

        if not positive:
            intervals = [(lcontig - x[1], lcontig - x[0])
                         for x in intervals[::-1]]
            out.reverse()

        s = [fasta.getSequence(contig, strand, start, end)
             for start, end in intervals]
        # IMS: allow for masking of sequences
        s = Masker.maskSequences(s, options.masker)
        l = sum([len(x) for x in s])
        if (l < options.min_length or
                (options.max_length and l > options.max_length)):
            nskipped_length += 1
            if options.loglevel >= 1:
                options.stdlog.write("# skipped because length out of bounds "
                                     "%s: regions=%s len=%i\n" %
                                     (name, str(intervals), l))
                continue

        if options.extend_at and options.extend_with:
            extension = "".join((options.extend_with,) * options.extend_by)

            if options.extend_at in ("5", "both"):
                s[1] = extension + s[1]
            if options.extend_at in ("3", "both"):
                s[-1] = s[-1] + extension

        if options.fold_at:
            n = options.fold_at
            s = "".join(s)
            seq = "\n".join([s[i:i+n] for i in range(0, len(s), n)])
        else:
            seq = "\n".join(s)

        if options.header_attr:
            attributes = " ".join([":".join([ax, ay]) for ax, ay in chunk[0].asDict().items()])
            options.stdout.write(">%s %s:%s:%s feature:%s %s\n%s\n" % (name,
                                                                       contig,
                                                                       strand,
                                                                       ";".join(
                                                                           ["%i-%i" %
                                                                            x for x in out]),
                                                                       chunk[0].feature,
                                                                       attributes,
                                                                       seq))
        else:
            options.stdout.write(">%s %s:%s:%s\n%s\n" % (name,
                                                         contig,
                                                         strand,
                                                         ";".join(
                                                             ["%i-%i" %
                                                              x for x in out]),
                                                         seq))

        noutput += 1

    E.info("ninput=%i, noutput=%i, nmasked=%i, nskipped_noexons=%i, "
           "nskipped_masked=%i, nskipped_length=%i" %
           (ninput, noutput, nmasked, nskipped_noexons,
            nskipped_masked, nskipped_length))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
