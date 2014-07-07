'''
split_gff - split a gff file into chunks
=============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics Intervals Genesets GFF Manipulation

Purpose
-------

Split gff file into chunks. Overlapping entries will 
always be output in the same chunk. Input is read from
stdin unless otherwise specified. The input needs to
be contig/start position sorted.

Options
-------

-p --output-pattern

    Specifies how to name the output files. The placeholder
    should be included and will specify where the chunk
    number will go. 

    E.g. -p output_chunk_%i.gff 

    will produce output files:
    output_chunck_1.gff
    output_chunck_2.gff
    output_chunck_3.gff
    ...


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

-I, -L 

    Use these options to redirect the stdin, stdout and
    stderror to files. Gziped files will be handled
    transparently

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

Will read a gff file from stdin and split into multiple
gff files. 

   cgat split_gff -I GFF [OPTIONS]

Will read the gff file GFF and split into multiple gff
files. 

Command line options
--------------------

'''

import sys
import re
import os
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.Experiment as E


class OutputChunk:

    def __init__(self, options):
        self.nchunk = 0
        self.options = options

    def createOpen(self, mode="w", header=None):
        """open file. Check first, if directory exists.
        """

        self.nchunk += 1
        filename = self.options.output_pattern % self.nchunk

        if self.options.dry_run:
            E.info("opening file %s" % filename)
            return open("/dev/null", mode)

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
        version="%prog version: $Id: gff2chunks.py 2781 2009-09-10 11:33:14Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-p", "--output-pattern", dest="output_pattern", type="string",
                      help="output pattern for filenames. Should contain a '%i' [default=%default].")

    parser.add_option("-i", "--min-chunk-size", dest="min_chunk_size", type="int",
                      help="minimum chunk size [default=%default].")

    parser.add_option("-n", "--dry-run", dest="dry_run", action="store_true",
                      help="do not create any files [default=%default].")

    parser.set_defaults(
        method="overlap",
        dry_run=False,
        output_pattern="%06i.chunk",
        min_chunk_size=2,
    )

    (options, args) = E.Start(parser)

    gffs = GTF.iterator(options.stdin)

    ninput, noutput, nchunks = 0, 0, 0

    outputChunk = OutputChunk(options)

    if options.method == "overlap":

        last_contig, last_to = None, None
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
