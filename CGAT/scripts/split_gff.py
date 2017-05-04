'''split_gff - split a gff file into chunks
=============================================

:Tags: Genomics Intervals Genesets GFF Manipulation

Purpose
-------

Split gff file into chunks. Overlapping entries will always be output
in the same chunk. Input is read from stdin unless otherwise
specified. The input needs to be contig/start position sorted.

Options
-------

-i --min-chunk-size

    This option specifies how big each chunck should
    be, in terms of the number of gff lines to be
    included. Because overlapping lines are always
    output to the same file, this should be considered
    a minimum size.

-n, --dry-run

    This options tells the script not to actaully write
    any files, but it will output a list of the files
    that would be output.

Example
-------

    cgat splitgff -i 1 < in.gff

where in.gff looks like:

    chr1	.	exon	1	10	.	+	.
    chr1	.	exon	8	100	.	+	.
    chr1	.	exon	102	150	.	+	.

will produce two files that look like:

    000001.chunk:
    chr1	.	exon	1	10	.	+	.
    chr1	.	exon	8	100	.	+	.

    000002.chunk:
    chr1	.	exon	102	150	.	+	.

Usage
-----

   cgat splitgff [OPTIONS]

Will read a gff file from stdin and split into multiple gff files.

   cgat split_gff -I GFF [OPTIONS]

Will read the gff file GFF and split into multiple gff files.

Command line options
--------------------

'''

import sys
import os
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.Experiment as E


class OutputChunk:

    def __init__(self, output_filename_pattern, dry_run=False):
        self.nchunk = 0
        self.output_filename_pattern = output_filename_pattern
        self.dry_run = dry_run

    def createOpen(self, mode="w", header=None):
        """open file. Check first, if directory exists.
        """

        self.nchunk += 1
        filename = self.output_filename_pattern % self.nchunk

        if self.dry_run:
            E.info("opening file %s" % filename)
            returnIOTools.openFile("/dev/null", mode)

        if mode in ("w", "a"):
            dirname = os.path.dirname(filename)
            if dirname and not os.path.exists(dirname):
                os.makedirs(dirname)

        if os.path.exists(filename):
            existed = True
        else:
            existed = False

        f = IOTools.openFile(filename, mode)

        if header and not existed:
            f.write(header + "\n")

        return f

    def __call__(self, chunk):
        """output a chunk into a new file."""
        outfile = self.createOpen()
        for c in chunk:
            outfile.write(str(c) + "\n")
        outfile.close()
        return len(chunk)


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option(
        "-i", "--min-chunk-size", dest="min_chunk_size", type="int",
        help="minimum chunk size [default=%default].")

    parser.add_option(
        "-n", "--dry-run", dest="dry_run", action="store_true",
        help="do not create any files [default=%default].")

    parser.set_defaults(
        method="overlap",
        dry_run=False,
        min_chunk_size=2,
        output_filename_pattern="%06i.chunk",
    )

    (options, args) = E.Start(parser, add_output_options=True)

    gffs = GTF.iterator(options.stdin)

    ninput, noutput, nchunks = 0, 0, 0

    outputChunk = OutputChunk(options.output_filename_pattern,
                              dry_run=options.dry_run)

    if options.method == "overlap":

        last_contig, last_to = None, 0
        chunk = []
        for gff in gffs:
            ninput += 1
            if len(chunk) >= options.min_chunk_size and \
                    (gff.contig != last_contig or
                     gff.start > last_to):
                noutput += outputChunk(chunk)
                nchunks += 1
                chunk = []
                last_contig, last_to = gff.contig, gff.end

            chunk.append(gff)
            last_to = max(gff.end, last_to)

        noutput += outputChunk(chunk)
        nchunks += 1

    E.info("ninput=%i, noutput=%i, nchunks=%i" % (ninput, noutput, nchunks))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
