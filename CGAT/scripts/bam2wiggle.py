"""bam2wiggle.py - convert bam to wig/bigwig file
==============================================

:Tags: Genomics NGS Intervals Conversion BAM WIGGLE BIGWIG BEDGRAPH

Purpose
-------

convert a bam file to a bigwig or bedgraph file.

Depending on options chosen, this script either computes the densities
itself or makes use of faster solutions if possible. The script
requires the executables :file:`wigToBigWig` and :file:`bedToBigBed`
to be in the user's PATH.

If no --shift-size or --extend option are given, the coverage is computed
directly on reads.  Counting can be performed at a certain resolution.

The counting currently is not aware of spliced reads, i.e., an
inserted intron will be included in the coverage.

If --shift-size or --extend are given, the coverage is computed by shifting
read alignment positions upstream for positive strand reads or
downstream for negative strand reads and extend them by a fixed
amount.

For RNASEQ data it might be best to run genomeCoverageBed directly on
the bam file.

Usage
-----

Type::

   cgat bam2wiggle \
          --output-format=bigwig \
          --output-filename-pattern=out.bigwig in.bam

to convert the :term:`bam` file file:`in.bam` to :term:`bigwig` format
and save the result in :file:`out.bigwig`.

Command line options
--------------------

"""

import os
import sys
import tempfile
import shutil
import subprocess
import CGAT.Experiment as E
import pysam
import CGAT.IOTools as IOTools
import CGAT.scripts._bam2bed as _bam2bed


class SpanWriter(object):

    '''output values within spans.
    values are collected according to a span and an average is
    output.
    '''

    def __init__(self, span):
        self.span = span
        self.laststart = 0
        self.lastend = None
        self.val = 0
        self.lastout = None

    def __call__(self, outfile, contig, start, end, val):

        # deal with previous window
        if self.lastend:
            last_window_start = self.lastend - self.lastend % self.span
            last_window_end = last_window_start + self.span
        else:
            last_window_start = start - start % self.span
            last_window_end = start + self.span

        # print start, end, val, last_window_start, last_window_end

        if self.lastend and start > last_window_end:
            # no overlap, output previous span
            assert self.lastout != last_window_start, \
                ("start=%i, end=%i, laststart=%i, "
                 "lastend=%i, last_window_start=%i") % \
                (start, end, self.laststart,
                 self.lastend, last_window_start)
            self.lastout = last_window_start
            v = self.val / float(self.span)
            outfile.write("%i\t%f\n" % (last_window_start, v))
            self.val = 0

            last_window_start = start - start % self.span
            last_window_end = start + self.span

        if end < last_window_end:
            # window too small to output, simply add values
            self.val += val * (end - start)
            self.lastend = max(end, self.lastend)
        else:
            # output first window
            v = self.val + val * \
                (self.span - start % self.span) / float(self.span)

            s = last_window_start
            assert self.lastout != s, \
                "start=%i, end=%i, laststart=%i, lastend=%i, s=%i" % \
                (start, end, self.laststart, self.lastend, s)
            outfile.write("%i\t%f\n" % (s, v))
            self.lastout = s
            self.val = 0

            # Output middle windows
            for x in range(start + self.span - start % self.span,
                           end - self.span, self.span):
                assert self.lastout != x
                outfile.write("%i\t%f\n" % (x, val))
                self.lastout = x

            if end % self.span:
                # save rest
                self.lastend = end
                self.val = val * end % self.span
            elif end - self.span != last_window_start:
                # special case, end ends on window
                assert self.lastout != end - self.span, \
                    "start=%i, end=%i, laststart=%i, lastend=%i" % \
                    (start, end, self.laststart, self.lastend)
                self.lastout = end - self.span
                outfile.write("%i\t%f\n" % (end - self.span, val))
                self.lastend = None
            else:
                # special case, end ends on window and only single
                # window - already output as start
                self.lastend = None

    def flush(self, outfile):
        if self.lastend:
            outfile.write("%i\t%f\n" % (
                self.lastend - self.lastend % self.span,
                self.val / (self.lastend % self.span)))


def main(argv=None):
    """script main.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option("-o", "--output-format", dest="output_format",
                      type="choice",
                      choices=(
                          "bedgraph", "wiggle", "bigbed",
                          "bigwig", "bed"),
                      help="output format [default=%default]")

    parser.add_option("-s", "--shift-size", dest="shift", type="int",
                      help="shift reads by a certain amount (ChIP-Seq) "
                      "[%default]")

    parser.add_option("-e", "--extend", dest="extend", type="int",
                      help="extend reads by a certain amount "
                      "(ChIP-Seq) [%default]")

    parser.add_option("-p", "--wiggle-span", dest="span", type="int",
                      help="span of a window in wiggle tracks "
                      "[%default]")

    parser.add_option("-m", "--merge-pairs", dest="merge_pairs",
                      action="store_true",
                      help="merge paired-ended reads into a single "
                      "bed interval [default=%default].")

    parser.add_option("--scale-base", dest="scale_base", type="float",
                      help="number of reads/pairs to scale bigwig file to. "
                      "The default is to scale to 1M reads "
                      "[default=%default]")

    parser.add_option("--scale-method", dest="scale_method", type="choice",
                      choices=("none", "reads",),
                      help="scale bigwig output. 'reads' will normalize by "
                      "the total number reads in the bam file that are used "
                      "to construct the bigwig file. If --merge-pairs is used "
                      "the number of pairs output will be used for "
                      "normalization. 'none' will not scale the bigwig file"
                      "[default=%default]")

    parser.add_option("--max-insert-size", dest="max_insert_size",
                      type="int",
                      help="only merge if insert size less that "
                      "# bases. 0 turns of this filter "
                      "[default=%default].")

    parser.add_option("--min-insert-size", dest="min_insert_size",
                      type="int",
                      help="only merge paired-end reads if they are "
                      "at least # bases apart. "
                      "0 turns of this filter. [default=%default]")

    parser.set_defaults(
        samfile=None,
        output_format="wiggle",
        shift=0,
        extend=0,
        span=1,
        merge_pairs=None,
        min_insert_size=0,
        max_insert_size=0,
        scale_method='none',
        scale_base=1000000,
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv, add_output_options=True)

    if len(args) >= 1:
        options.samfile = args[0]
    if len(args) == 2:
        options.output_filename_pattern = args[1]
    if not options.samfile:
        raise ValueError("please provide a bam file")

    # Read BAM file using Pysam
    samfile = pysam.AlignmentFile(options.samfile, "rb")

    # Create temporary files / folders
    tmpdir = tempfile.mkdtemp()
    E.debug("temporary files are in %s" % tmpdir)
    tmpfile_wig = os.path.join(tmpdir, "wig")
    tmpfile_sizes = os.path.join(tmpdir, "sizes")

    # Create dictionary of contig sizes
    contig_sizes = dict(list(zip(samfile.references, samfile.lengths)))
    # write contig sizes
    outfile_size = IOTools.openFile(tmpfile_sizes, "w")
    for contig, size in sorted(contig_sizes.items()):
        outfile_size.write("%s\t%s\n" % (contig, size))
    outfile_size.close()

    # Shift and extend only available for bigwig format
    if options.shift or options.extend:
        if options.output_format != "bigwig":
            raise ValueError(
                "shift and extend only available for bigwig output")

    # Output filename required for bigwig / bigbed computation
    if options.output_format == "bigwig":
        if not options.output_filename_pattern:
            raise ValueError(
                "please specify an output file for bigwig computation.")

        # Define executable to use for binary conversion
        if options.output_format == "bigwig":
            executable_name = "wigToBigWig"
        else:
            raise ValueError("unknown output format `%s`" %
                             options.output_format)

        # check required executable file is in the path
        executable = IOTools.which(executable_name)
        if not executable:
            raise OSError("could not find %s in path." % executable_name)

        # Open outout file
        outfile = IOTools.openFile(tmpfile_wig, "w")
        E.info("starting output to %s" % tmpfile_wig)
    else:
        outfile = IOTools.openFile(tmpfile_wig, "w")
        E.info("starting output to stdout")

    # Set up output write functions
    if options.output_format in ("wiggle", "bigwig"):
        # wiggle is one-based, so add 1, also step-size is 1, so need
        # to output all bases
        if options.span == 1:
            outf = lambda outfile, contig, start, end, val: \
                outfile.write(
                    "".join(["%i\t%i\n" % (x, val)
                             for x in range(start + 1, end + 1)]))
        else:
            outf = SpanWriter(options.span)
    elif options.output_format == "bedgraph":
        # bed is 0-based, open-closed
        outf = lambda outfile, contig, start, end, val: \
            outfile.write("%s\t%i\t%i\t%i\n" % (contig, start, end, val))

    # initialise counters
    ninput, nskipped, ncontigs = 0, 0, 0

    # set output file name
    output_filename_pattern = options.output_filename_pattern
    if output_filename_pattern:
        output_filename = os.path.abspath(output_filename_pattern)

    # shift and extend or merge pairs. Output temporay bed file
    if options.shift > 0 or options.extend > 0 or options.merge_pairs:
        # Workflow 1: convert to bed intervals and use bedtools
        # genomecov to build a coverage file.
        # Convert to bigwig with UCSC tools bedGraph2BigWig

        if options.merge_pairs:
            # merge pairs using bam2bed
            E.info("merging pairs to temporary file")
            counter = _bam2bed.merge_pairs(
                samfile,
                outfile,
                min_insert_size=options.min_insert_size,
                max_insert_size=options.max_insert_size,
                bed_format=3)
            E.info("merging results: {}".format(counter))
            if counter.output == 0:
                raise ValueError("no pairs output after merging")
        else:
            # create bed file with shifted/extended tags
            shift, extend = options.shift, options.extend
            shift_extend = shift + extend
            counter = E.Counter()

            for contig in samfile.references:
                E.debug("output for %s" % contig)
                lcontig = contig_sizes[contig]

                for read in samfile.fetch(contig):
                    pos = read.pos
                    if read.is_reverse:
                        start = max(0, read.pos + read.alen - shift_extend)
                    else:
                        start = max(0, read.pos + shift)

                    # intervals extending beyond contig are removed
                    if start >= lcontig:
                        continue

                    end = min(lcontig, start + extend)
                    outfile.write("%s\t%i\t%i\n" % (contig, start, end))
                    counter.output += 1

        outfile.close()

        if options.scale_method == "reads":
            scale_factor = float(options.scale_base) / counter.output

            E.info("scaling: method=%s scale_quantity=%i scale_factor=%f" %
                   (options.scale_method,
                    counter.output,
                    scale_factor))
            scale = "-scale %f" % scale_factor
        else:
            scale = ""

        # Convert bed file to coverage file (bedgraph)
        tmpfile_bed = os.path.join(tmpdir, "bed")
        E.info("computing coverage")
        # calculate coverage - format is bedgraph
        statement = """bedtools genomecov -bg -i %(tmpfile_wig)s %(scale)s
        -g %(tmpfile_sizes)s > %(tmpfile_bed)s""" % locals()
        E.run(statement)

        # Convert bedgraph to bigwig
        E.info("converting to bigwig")
        tmpfile_sorted = os.path.join(tmpdir, "sorted")
        statement = ("sort -k 1,1 -k2,2n %(tmpfile_bed)s > %(tmpfile_sorted)s;"
                     "bedGraphToBigWig %(tmpfile_sorted)s %(tmpfile_sizes)s "
                     "%(output_filename_pattern)s" % locals())
        E.run(statement)

    else:

        # Workflow 2: use pysam column iterator to build a
        # wig file. Then convert to bigwig of bedgraph file
        # with UCSC tools.
        def column_iter(iterator):
            start = None
            end = 0
            n = None
            for t in iterator:
                if t.pos - end > 1 or n != t.n:
                    if start is not None:
                        yield start, end, n
                    start = t.pos
                    end = t.pos
                    n = t.n
                end = t.pos
            yield start, end, n

        if options.scale_method != "none":
            raise NotImplementedError(
                "scaling not implemented for pileup method")

        # Bedgraph track definition
        if options.output_format == "bedgraph":
            outfile.write("track type=bedGraph\n")

        for contig in samfile.references:
            # if contig != "chrX": continue
            E.debug("output for %s" % contig)
            lcontig = contig_sizes[contig]

            # Write wiggle header
            if options.output_format in ("wiggle", "bigwig"):
                outfile.write("variableStep chrom=%s span=%i\n" %
                              (contig, options.span))

            # Generate pileup per contig using pysam and iterate over columns
            for start, end, val in column_iter(samfile.pileup(contig)):
                # patch: there was a problem with bam files and reads
                # overextending at the end. These are usually Ns, but
                # need to check as otherwise wigToBigWig fails.
                if lcontig <= end:
                    E.warn("read extending beyond contig: %s: %i > %i" %
                           (contig, end, lcontig))
                    end = lcontig
                    if start >= end:
                        continue

                if val > 0:
                    outf(outfile, contig, start, end, val)
            ncontigs += 1

        # Close output file
        if type(outf) == type(SpanWriter):
            outf.flush(outfile)
        else:
            outfile.flush()

        E.info("finished output")

        # Report counters
        E.info("ninput=%i, ncontigs=%i, nskipped=%i" %
               (ninput, ncontigs, nskipped))

        # Convert to binary formats
        if options.output_format == "bigwig":
            outfile.close()

            E.info("starting %s conversion" % executable)
            try:
                retcode = subprocess.call(
                    " ".join((executable,
                              tmpfile_wig,
                              tmpfile_sizes,
                              output_filename_pattern)),
                    shell=True)
                if retcode != 0:
                    E.warn("%s terminated with signal: %i" %
                           (executable, -retcode))
                    return -retcode
            except OSError as msg:
                E.warn("Error while executing bigwig: %s" % msg)
                return 1
            E.info("finished bigwig conversion")
        else:
            with open(tmpfile_wig) as inf:
                sys.stdout.write(inf.read())

    # Cleanup temp files
    shutil.rmtree(tmpdir)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
