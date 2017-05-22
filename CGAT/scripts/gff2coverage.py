'''
gff2coverage.py - compute genomic coverage of gff intervals
===========================================================

:Tags: Genomics Intervals Summary GFF

Purpose
-------

This script computes the genomic coverage of intervals in
a :term:`gff` formatted file. The coverage is computed per feature.

Usage
-----

You can use two methods to compute the coverage: **genomic** and **histogram**.

Let us explain their usage with this ``small.gtf`` file::

   19 processed_transcript exon 16 16 . - . gene_id
   19 processed_transcript exon 27 27 . - . gene_id
   19 processed_transcript exon 8  8  . - . gene_id
   19 processed_transcript exon 19 19 . - . gene_id
   19 processed_transcript exon 5  5  . - . gene_id

and this toy example (``small.fasta``) of an indexed fasta file::

   >chr19
   GCCGGCCTCTACCTGCAGCAGATGCCCTAT

Both files (``small.gtf`` and ``small.fasta``) are included
in the `GitHub <https://github.com/CGATOxford/cgat>`_ repository.

genomic method
++++++++++++++

The **genomic** method computes the coverage of intervals
accross the genome file given as input. Let us see how to
apply the genomic method to the small examples above::

   python gff2coverage.py --method=genomic --genome-file=small < small.gtf

The output (wrapped to fit here) will be::

   contig  source  feature  intervals  bases  p_coverage  total_p_coverage
   19      trans.  exon     5          5      16.666667   16.666667

As you can see the information displayed is the following: contig name,
source, feature name, number of intervals within the contig, number of
bases, percentage of coverage in the contig, and percentage of coverage
in the genome file.

histogram method
++++++++++++++++

On the contrary, if you want to compute the coverage of intervals
within the :term:`gff` file itself summarized as an histogram and
grouped by contig name, please use the histogram method.

To use the histogram method with the input files above, please type::

 python gff2coverage.py\
 --method=histogram\
 --window=5\
 --features=exon\
 --output-filename-pattern=%s.hist < small.gtf

In this case the output (written to file ``19.hist``) is::

   abs_pos  rel_pos  abs_exon  rel_exon
   0        0.0000   1         0.2000
   5        0.1852   2         0.4000
   10       0.3704   2         0.4000
   15       0.5556   4         0.8000
   20       0.7407   4         0.8000
   25       0.9259   5         1.0000

The output is given as a pair of columns. The first pair of columns always
appears and lists the cumulative numbers of nucleotides in each window or
bin --absolute and relative values in the former and latter columns,
respectively. The subsequent pair of columns depends on the values given to
the ``--features`` option. In this example there is an extra column for the
``exon`` feature but you could especify as many of them as you wanted among
those features listed in your :term:`gff` file.

On the other hand, the ``--num-bins`` option can be used instead of
``--window`` along with ``--genome-file`` to define the number of bins for the
resultant histogram. This parameter is used by default (with value: 1000)
when using the histogram method.

Please note the following:

- you need to specify the feature name explicitly (with the ``--feature`` \
option) to compute the genomic coverage of that feature. You can also use\
a comma-separated list of feature names.

- the output of the histogram method goes to a file (in the current working\
directory) which is named as the contig name by default. To change this\
behaviour, please use the ``--output-filename-pattern`` option where \
``%s`` will be substituted by the contig name.

Command line options
--------------------
'''

import sys
import math
import collections

import CGAT.Experiment as E
import CGAT.IndexedFasta as IndexedFasta
import CGAT.GTF as GTF


def printValues(contig, max_size, window_size, values, options):
    """output values."""

    outfile = E.openOutputFile(contig, "w")

    outfile.write("abs_pos\trel_pos")

    for feature in options.features:
        outfile.write("\tabs_%s\trel_%s" % (feature, feature))
    outfile.write("\n")

    max_vv = []

    for f in range(len(options.features)):
        max_vv.append(float(max([x[f] for x in values])))

    bin = 0
    for vv in values:
        outfile.write("%i\t" % bin)
        outfile.write(options.value_format % (float(bin) / max_size))

        for x in range(len(options.features)):
            outfile.write("\t%i\t%s" % (
                vv[x],
                options.value_format % (vv[x] / max_vv[x])))
        outfile.write("\n")
        bin += window_size

    outfile.close()


def processChunk(contig, chunk, options, fasta=None):
    """
    This function requires segments to be non-overlapping.
    """

    if len(chunk) == 0:
        return

    # check whether there are overlapping features or not
    checked = []
    for feature in chunk:
        checked.append(feature)
        others = [x for x in chunk if x not in checked]
        for otherFeature in others:
            if GTF.Overlap(feature, otherFeature):
                raise ValueError(" Histogram could not be created"
                                 " since the file contains overlapping "
                                 "features! \n%s\n%s  "
                                 % (feature, otherFeature))
    # clear auxiliary list
    del checked[:]

    # compute max_coordinate for the histogram
    max_coordinate = max([x.end for x in chunk])
    # compute window size
    if options.window_size:
        window_size = options.window_size
        num_bins = int(math.ceil((float(max_coordinate) / window_size)))
    elif options.num_bins and fasta:
        contig_length = fasta.getLength(contig)
        assert max_coordinate <= contig_length, ("maximum coordinate (%i) "
                                                 "larger than contig size (%i)"
                                                 " for contig %s"
                                                 % (max_coordinate,
                                                    contig_length,
                                                    contig))
        max_coordinate = contig_length
        window_size = int(math.floor(float(contig_length) / options.num_bins))
        num_bins = options.num_bins
    else:
        raise ValueError("please specify a window size of provide "
                         "genomic sequence with number of bins.")

    values = [[] for x in range(num_bins)]

    # do several parses for each feature, slow, but easier to code
    # alternatively: sort by feature and location.
    for feature in options.features:
        total = 0
        bin = 0
        end = window_size
        for entry in chunk:
            if entry.feature != feature:
                continue

            while end < entry.start:
                values[bin].append(total)
                bin += 1
                end += window_size

            while entry.end > end:
                seg_start = max(entry.start, end - window_size)
                seg_end = min(entry.end, end)
                total += seg_end - seg_start
                values[bin].append(total)
                end += window_size
                bin += 1
            else:
                seg_start = max(entry.start, end - window_size)
                seg_end = min(entry.end, end)
                total += seg_end - seg_start

        while bin < num_bins:
            values[bin].append(total)
            bin += 1

    printValues(contig, max_coordinate, window_size, values, options)


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: "
                            "$Id: gff2coverage.py 2781 2009-09-10 11:33:14Z "
                            "andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default]")

    parser.add_option("-f", "--features", dest="features", type="string",
                      action="append", help="features to collect "
                      "[default=%default]")

    parser.add_option("-w", "--window-size", dest="window_size", type="int",
                      help="window size in bp for histogram computation. "
                      "Determines the bin size.  "
                      "[default=%default]")

    parser.add_option("-b", "--num-bins", dest="num_bins", type="int",
                      help="number of bins for histogram computation "
                      "if window size is not given. "
                      "[default=%default]")

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=("genomic", "histogram", ),
                      help="methods to apply. "
                      "[default=%default]")

    parser.set_defaults(
        genome_file=None,
        window_size=None,
        num_bins=1000,
        value_format="%6.4f",
        features=[],
        method="genomic",
    )

    (options, args) = E.Start(parser, add_output_options=True)

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta(options.genome_file)
    else:
        fasta = None

    if options.method == "histogram":

        gff = GTF.readFromFile(options.stdin)

        gff.sort(key=lambda x: (x.contig, x.start))

        chunk = []
        last_contig = None

        for entry in gff:

            if last_contig != entry.contig:
                processChunk(last_contig, chunk, options, fasta)
                last_contig = entry.contig
                chunk = []

            chunk.append(entry)

        processChunk(last_contig, chunk, options, fasta)

    elif options.method == "genomic":
        intervals = collections.defaultdict(int)
        bases = collections.defaultdict(int)
        total = 0
        for entry in GTF.iterator(options.stdin):
            intervals[(entry.contig, entry.source, entry.feature)] += 1
            bases[(entry.contig, entry.source, entry.feature)
                  ] += entry.end - entry.start
            total += entry.end - entry.start

        options.stdout.write("contig\tsource\tfeature\tintervals\tbases")
        if fasta:
            options.stdout.write(
                "\tpercent_coverage\ttotal_percent_coverage\n")
        else:
            options.stdout.write("\n")

        total_genome_size = sum(
            fasta.getContigSizes(with_synonyms=False).values())

        for key in sorted(intervals.keys()):
            nbases = bases[key]
            nintervals = intervals[key]
            contig, source, feature = key
            options.stdout.write("\t".join(("\t".join(key),
                                            str(nintervals),
                                            str(nbases))))
            if fasta:
                options.stdout.write(
                    "\t%f" % (100.0 * float(nbases) / fasta.getLength(contig)))
                options.stdout.write(
                    "\t%f\n" % (100.0 * float(nbases) / total_genome_size))
            else:
                options.stdout.write("\n")

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
