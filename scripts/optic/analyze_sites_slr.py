##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################
'''
optic/analyze_sites_slr.py -
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python optic/analyze_sites_slr.py --help

Type::

   python optic/analyze_sites_slr.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import re
from types import *
import CGAT.Genomics as Genomics
import CGAT.Experiment as E
import scipy
import scipy.stats
import CGAT.WrapperSlr as WrapperSlr
import CGAT.Mali as Mali
import CGAT.CSV as CSV

USAGE = """analyze results from site-specific codeml run.

The input is either:

   * filenames with a set of files of related codeml runs

"""


class Result:

    def __init__(self, id=0, nfiltered=0, ntotal=0, nsynonymous=0, nnegative=0, npositive=0, significance=0):
        self.mId = id
        self.mNPositive = npositive
        self.mNNegative = nnegative
        self.mNSynonymous = nsynonymous
        self.mNFiltered = nfiltered
        self.mNTotal = ntotal
        self.mSignificance = significance

    def __str__(self):
        return "%s\t%i\t%i\t%i\t%i\t%i\t%5.2e" % (self.mId, self.mNFiltered, self.mNTotal, self.mNSynonymous, self.mNNegative, self.mNPositive, self.mSignificance)

    def getHeader(self):
        return "\t".join(("prefix", "nfiltered", "ntotal", "nsyn", "nneg", "npos", "p"))


class PositionInformation:

    def __init__(self, aa, seq_pos, mali_pos, context):
        self.mAA = aa
        self.mSequencePosition = seq_pos
        self.mContext = context
        self.mMaliPosition = mali_pos


def ProcessResult(result, options, mali=None, prefix=None, p_value=None):

    counts = None

    if options.method == "summary-slr":

        thresholds = "95%", "99%", "95% corrected", "99% corrected"

        if prefix:
            options.stdout.write("%s\t" % prefix)

        options.stdout.write("%5.2f\t%5.2f\t%5.2f\t%6.4f\t%i\t%i\t%i\t" % (result.mTreeLength,
                                                                           result.mOmega,
                                                                           result.mKappa,
                                                                           result.mLogLikelihood,
                                                                           len(result.mSites),
                                                                           result.mNSitesSynonymous,
                                                                           result.mNSitesGaps +
                                                                           result.mNSitesSingleChar,
                                                                           ))
        options.stdout.write(
            "\t".join(map(lambda x: "%i" % result.mNPositiveSites[x][0], thresholds)))
        options.stdout.write("\t")
        options.stdout.write(
            "\t".join(map(lambda x: "%i" % result.mNNegativeSites[x], thresholds)))
        options.stdout.write("\n")

    elif options.method in ("summary-filtered",
                            "positive-site-table", "negative-site-table", "neutral-site-table",
                            "positive-site-list", "negative-site-list", "neutral-site-list"):

        mali_length = mali.getLength()
        mali_width = mali.getWidth()
        column_data = map(
            lambda x: Mali.MaliData(x, gap_chars="Nn", mask_chars="-."), mali.getColumns())

        # sanity check: do lengths of mali and # of sites correspond
        if len(result.mSites) * 3 != mali_width:
            raise "mali (%i) and # of sites (%i) do not correspond." % (
                mali_width, len(result.mSites))

        if options.method == "summary-filtered":
            # count sites, but filter with multiple alignment
            ntotal = 0
            npositive = 0
            nnegative = 0
            nneutral = 0
            nfiltered = 0
            nsynonymous = 0

            if prefix:
                options.stdout.write("%s\t" % prefix)

            for x in range(len(result.mSites)):
                site = result.mSites[x]
                column = column_data[x * 3]

                if column.mNChars != mali_length:
                    nfiltered += 1
                    continue

                if site.isPositive(options.significance_threshold, options.use_adjusted):
                    npositive += 1
                elif site.isNegative(options.significance_threshold, options.use_adjusted):
                    nnegative += 1

                if site.isSynonymous():
                    nsynonymous += 1

                ntotal += 1

            options.stdout.write("%5.2f\t%5.2f\t%5.2f\t%6.4f\t%i\t%i\t%i\t%i\t%i\t%i\n" % (result.mTreeLength,
                                                                                           result.mOmega,
                                                                                           result.mKappa,
                                                                                           result.mLogLikelihood,
                                                                                           len(result.mSites),
                                                                                           nfiltered, ntotal,
                                                                                           nsynonymous, nnegative, npositive))
            counts = Result(
                nfiltered, ntotal, nsynonymous, nnegative, npositive)

        elif options.method in ("positive-site-table", "negative-site-table", "neutral-site-table",
                                "positive-site-list", "negative-site-list", "neutral-site-list",
                                ):

            select_positive_sites = options.method in (
                "positive-site-table", "positive-site-list")
            select_negative_sites = options.method in (
                "negative-site-table", "negative-site-list")

            # iterate over sites and output those under xxx selection
            identifiers = mali.getIdentifiers()
            chars_per_row = [[] for x in range(mali_length)]

            sites = []

            for col in range(len(result.mSites)):

                site = result.mSites[col]
                column = column_data[col * 3]

                if column.mNChars != mali_length:
                    continue

                keep = False

                if select_positive_sites and site.isPositive(options.significance_threshold, options.use_adjusted):
                    keep = True

                elif select_negative_sites and site.isNegative(options.significance_threshold, options.use_adjusted):
                    keep = True

                if not keep:
                    continue

                sites.append((col, site))

            nsites = len(sites)

            if options.truncate_sites_list:
                # truncate sites list, sort by significance
                sites.sort(lambda x, y: cmp(x[1].mPValue, y[1].mPValue))
                sites = sites[:options.truncate_sites_list]

            for col, site in sites:

                site = result.mSites[col]
                xcol = col * 3

                for row in range(mali_length):
                    id = identifiers[row]
                    x = max(xcol - options.context_size * 3, 0)
                    y = min(xcol + 3 + options.context_size * 3, mali_width)
                    segment = mali[id][x:y]
                    codon = mali[id][xcol:xcol + 3]
                    pos = mali.getResidueNumber(id, xcol)
                    pos /= 3

                    # save as real-world coordinates
                    chars_per_row[row].append(PositionInformation(Genomics.MapCodon2AA(codon),
                                                                  pos + 1,
                                                                  xcol,
                                                                  Genomics.TranslateDNA2Protein(segment).upper()))

            if p_value is not None:
                pp_value = p_value
            else:
                pp_value = "na"

            if options.method in ("positive-site-table", "negative-site-table", "neutral-site-table"):

                if options.context_size:
                    for row in range(mali_length):
                        if prefix:
                            options.stdout.write("%s\t" % prefix)

                        options.stdout.write("%s\t%i\t%s\t%s\n" % (
                            identifiers[row],
                            nsites,
                            pp_value,
                            ";".join(["%s%i in %s" % (x.mAA, x.mSequencePosition, x.mContext) for x in chars_per_row[row]])))
                else:
                    for row in range(mali_length):
                        if prefix:
                            options.stdout.write("%s\t" % prefix)

                        options.stdout.write("%s\t%i\t%s\t%s\n" % (
                            identifiers[row],
                            nsites,
                            pp_value,
                            ";".join(["%s%i" % (x.mAA, x.mSequencePosition) for x in chars_per_row[row]])))

            elif options.method in ("positive-site-list", "negative-site-list", "neutral-site-list"):

                for row in range(mali_length):

                    if prefix:
                        xprefix = "%s\t%s" % (prefix, identifiers[row])
                    else:
                        xprefix = "%s" % (identifiers[row])
                    x = 0
                    for chars in chars_per_row[row]:
                        x += 1
                        options.stdout.write("%s\t%i\t%s\t%i\t%i\t%s\n" % (xprefix, x,
                                                                           chars.mAA, chars.mSequencePosition,
                                                                           chars.mMaliPosition, chars.mContext))

    options.stdout.flush()

    return counts


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/analyze_sites_slr.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("--method", dest="method", type="choice",
                      choices=("summary-slr", "summary-filtered", "over-representation",
                               "positive-site-table", "negative-site-table", "neutral-site-table",
                               "positive-site-list", "negative-site-list", "neutral-site-list"),
                      help="method to apply.")

    parser.add_option("--column-prefix", dest="prefix", type="string",
                      help="prefix for rows.")

    parser.add_option("-s", "--filename-sites", dest="filename_sites", type="string",
                      help="filename with sites information.")

    parser.add_option("-l", "--filename-log", dest="filename_log", type="string",
                      help="filename with logging information.")

    parser.add_option("-m", "--filename-mali", dest="filename_mali", type="string",
                      help="filename of multiple alignment, that was input to SLR. If given, is used to filter indels.")

    parser.add_option("--filter-probability", dest="filter_probability", type="float",
                      help="threshold for probability above which to include positive sites.")

    parser.add_option("--no-header", dest="write_header", action="store_false",
                      help="only output header.")

    parser.add_option("--only-header", dest="only_header", action="store_true",
                      help="only output header.")

    parser.add_option("--significance-threshold", dest="significance_threshold", type="float",
                      help="threshold for significance tests [%default].")

    parser.add_option("--use-adjusted", dest="use_adjusted", action="store_true",
                      help="use SLR adjusted probability values.")

    parser.add_option("--truncate-sites-list", dest="truncate_sites_list", type="int",
                      help="truncate sites list after ## entries (0 for all).")

    parser.add_option("--context-size", dest="context_size", type="int",
                      help="size of left/right context around a selected residue.")

    parser.set_defaults(
        prefix=None,
        filter_probability=0,
        filter_omega=0,
        filename_sites="-",
        filename_log=None,
        filename_mali=None,
        significance_threshold=0.05,
        write_header=True,
        only_header=False,
        use_adjusted=False,
        context_size=0,
        truncate_sites_list=0,
    )

    (options, args) = E.Start(parser)

    slr = WrapperSlr.Slr()

    # write headers
    if "%s" in options.filename_sites:
        options.prefix = True

    if options.method == "summary-slr":

        # write header
        if options.write_header or options.only_header:

            if options.loglevel >= 1:
                options.stdlog.write(
                    """# Numbers of positive/neutral/negative sites according to SLR
#
# This uses the thresholds as set in SLR. Use "counts" for filtering
# residues based on your own thresholds
""")
            thresholds = "95%", "99%", "95% corrected", "99% corrected"

            if options.prefix:
                options.stdout.write("prefix\t")
            options.stdout.write(
                "ltree\tomega\tkappa\tlnL\tnsites\tnsyn\tngap\t")
            options.stdout.write(
                "\t".join(map(lambda x: "npos_" + x.replace(" ", "_"), thresholds)))
            options.stdout.write("\t")
            options.stdout.write(
                "\t".join(map(lambda x: "nneg_" + x.replace(" ", "_"), thresholds)))
            options.stdout.write("\n")

    elif options.method == "summary-filtered":

        # write header
        if options.write_header or options.only_header:
            if options.loglevel >= 1:
                options.stdlog.write(
                    """# Numbers of positive/neutral/negative sites according to SLR
#
# This method uses the supplied threshold and the multiple alignment to filter.
# All positions that are above the threshold (P-Value) and which are located in
# indels: >= 1 sequence missing from column, are removed.
""")

            if options.prefix:
                options.stdout.write("prefix\t")
            options.stdout.write(
                "ltree\tomega\tkappa\tlnL\tnsites\tnfiltered\tntotal\tnsyn\tnneg\tnpos\n")

    elif options.method in ("positive-site-table", "negative-site-table", "neutral-site-table"):

        # write header
        if options.write_header or options.only_header:
            if options.loglevel >= 1:
                options.stdlog.write(
                    """# Numbers of positive/neutral/negative sites according to SLR
#
# Note: sequence positions are 1-based, but mali positions are 0-based.
# Residues in indel positions have been removed and signifnicance was determined according
# with a threshold of %5.2e
""" % options.significance_threshold)

            if options.prefix:
                options.stdout.write("prefix\t")
            options.stdout.write("cluster\tnsites\tp-value\tsites\n")

    elif options.method in ("positive-site-list", "negative-site-list", "neutral-site-list"):

        # write header
        if options.write_header or options.only_header:
            if options.loglevel >= 1:
                options.stdlog.write(
                    """# Sites under positive/neutral/negative selection according to SLR
#
# Note: sequence positions are 1-based, but mali positions are 0-based.
# Residues in indel positions have been removed and signifnicance was determined according
# with a threshold of %5.2e
""" % options.significance_threshold)

            if options.prefix:
                options.stdout.write("prefix\t")

            options.stdout.write(
                "sequence\tn\taa\tseq_pos\tmali_pos\tcontext\n")

    elif options.method == "over-representation":

        # write header
        if options.write_header or options.only_header:
            if options.loglevel >= 1:
                options.stdlog.write(
                    """# Genes with over-represented sites.
#
# This method uses as input the output of summary-filtered.
""")

    if options.only_header:
        sys.exit(0)

    if options.method in ("summary-slr", "summary-filtered",
                          "positive-site-table", "negative-site-table", "neutral-site-table",
                          "positive-site-list", "negative-site-list", "neutral-site-list"):

        ninput, noutput, nskipped = 0, 0, 0

        if "%s" in options.filename_sites:

            headers, table = CSV.ReadTable(sys.stdin)

            fprefix = headers.index("prefix")

            try:
                fsignificance = headers.index("p")
            except ValueError:
                fsignificance = None

            for row in table:

                id = row[fprefix]
                if fsignificance is not None:
                    p_value = row[fsignificance]
                else:
                    p_value = None

                ninput += 1

                fn = re.sub("%s", id, options.filename_sites)
                if not os.path.exists(fn):
                    nskipped += 1
                    continue

                lines_sites = open(fn, "r").readlines()
                if options.filename_log:
                    lines_log = open(
                        re.sub("%s", id, options.filename_log), "r").readlines()

                result = slr.parseOutput(lines_sites, lines_log)

                if options.method in ("summary-filtered", "positive-site-table", "negative-site-table", "neutral-site-table"):
                    mali = Mali.Mali()
                    mali.readFromFile(
                        open(re.sub("%s", id, options.filename_mali), "r"))
                else:
                    mali = None

                ProcessResult(result, options, mali,
                              prefix=id,
                              p_value=p_value)
                noutput += 1
        else:
            if options.filename_sites == "-":
                lines_sites = sys.stdin.readlines()
            else:
                lines_sites = open(options.filename_sites, "r").readlines()

            ninput += 1
            if options.filename_log:
                lines_log = open(options.filename_log, "r").readlines()

            result = slr.parseOutput(lines_sites, lines_log)

            if options.filename_mali:
                mali = Mali.Mali()
                mali.readFromFile(open(options.filename_mali, "r"))
            else:
                if options.method == "summary-filtered":
                    raise "please supply a multiple alignment for filtering."

                mali = None

            ProcessResult(result, options, mali, prefix=options.prefix)
            noutput += 1

        if options.loglevel >= 1:
            options.stdlog.write(
                "# ninput=%i, noutput=%i, nskipped=%i.\n" % (ninput, noutput, nskipped))

    else:
        if options.method == "over-representation":

            results = []
            for line in sys.stdin:
                if line[0] == "#":
                    continue
                data = line[:-1].split("\t")
                if data[0] == "prefix":
                    continue

                results.append(Result(
                    data[0], int(data[6]), int(data[7]), int(data[8]), int(data[9]), int(data[10])))

            # probability of a single site being positive
            ntotal = sum(map(lambda x: x.mNTotal, results))
            npositives = sum(map(lambda x: x.mNPositive, results))
            p = float(npositives) / float(ntotal)

            if options.loglevel >= 1:
                options.stdlog.write(
                    "# sites: total=%i, positive=%i, p=%f\n" % (ntotal, npositives, p))

            new_results = []
            for result in results:
                if result.mNTotal == 0:
                    continue

                # use -1, because I need P( x >= X)
                # sf = 1 - cdf and cdf = P( x <= X ), thus sf = 1 - P( x <= X )
                # = P (x > X ).
                r = scipy.stats.binom.sf(
                    result.mNPositive - 1, result.mNTotal, p)

                result.mSignificance = r

                if r < options.significance_threshold:
                    new_results.append(result)

            new_results.sort(
                lambda x, y: cmp(x.mSignificance, y.mSignificance))

            options.stdlog.write(Result().getHeader() + "\n")

            for result in new_results:
                options.stdout.write(str(result) + "\n")

            if options.loglevel >= 1:
                options.stdlog.write(
                    "# ntotal=%i, npos=%i\n" % (len(results), len(new_results)))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
