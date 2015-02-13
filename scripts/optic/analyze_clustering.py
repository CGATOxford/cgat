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
optic/analyze_clustering.py - 
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

   python optic/analyze_clustering.py --help

Type::

   python optic/analyze_clustering.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import re
import random
import CGAT.Experiment as E
import numpy
import scipy
import scipy.stats


"""analyse values based on a clustering.

input: a graph of similarities/distances between transcripts

output: statistics of the edge weights for all species pairs. Depending
on user choice (--method), summary statistics and/or histograms are
output.

optional input:
        a map of transcripts to clusters. This allows to filter:
        
        --between-clusters: only take values between clusters (default is within)

        --best-per-compoment: only take the best (minimum) value per cluster / cluster pair

"""

parser = E.OptionParser(
    version="%prog version: $Id: optic/analyze_clustering.py 2781 2009-09-10 11:33:14Z andreas $")


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser.add_option("-m", "--map", dest="map", type="string",
                      help="filename of map with ids to components.")
    parser.add_option("-e", "--method", dest="methods", type="choice", action="append",
                      choices=("hists", "stats"),
                      help="methods to apply to dataset.")
    parser.add_option("-s", "--pattern-species", dest="species_pattern", type="string",
                      help="regular expression to extract species [default='%default%'].")
    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to take for the analysis.")
    parser.add_option("-p", "--pattern-identifier", dest="pattern", type="string",
                      help="output pattern. Should contain two %s wildcards, one for the secion, the second for the method [default='%default'].")
    parser.add_option("-b", "--best-per-component", dest="best_per_component", action="store_true",
                      help="only take smallest value for each species pair in a component.")
    parser.add_option("-t", "--add-total", dest="add_total", action="store_true",
                      help="compute everything and add it.")
    parser.add_option("-o", "--only-total", dest="only_total", action="store_true",
                      help="only compute total (not per species).")
    parser.add_option("-a", "--maxbin", dest="maxbin", type="float",
                      help="maximum bin.")
    parser.add_option("-i", "--minbin", dest="minbin", type="float",
                      help="minimum bin.")
    parser.add_option("-n", "--numbin", dest="numbins", type="float",
                      help="number of bins.")
    parser.add_option("-z", "--binsize", dest="binsize", type="float",
                      help="size of a bin.")
    parser.add_option("--max-pairs", dest="max_pairs", type="int",
                      help="limit size to xxx data points per pair.")
    parser.add_option("--no-titles", dest="notitles", action="store_true",
                      help="no column titles.")

    parser.set_defaults(
        map=None,
        species_pattern="^([^@:|]+)[@:|]",
        columns="3",
        pattern="analysis_%s.%s",
        methods=[],
        numbins=20,
        minbin=0.0,
        maxbin=10.0,
        binsize=0.1,
        best_per_component=False,
        between_clusters=False,
        add_total=False,
        only_total=False,
        max_pairs=1000000,
        notitles=False,
    )

    (options, args) = E.Start(parser)

    options.columns = map(lambda x: int(x) - 1, options.columns.split(","))

    if not options.methods:
        raise "please supply at least one method to apply."

    if options.map:
        map_id2component = {}
        lines = filter(
            lambda x: x[0] != "#", open(options.map, "r").readlines())
        for i, c in map(lambda x: x[:-1].split("\t"), lines):
            map_id2component[i] = c

    arrays = {}

    rx_species = re.compile(options.species_pattern)
    species = {}
    component = "d"
    nskipped = 0
    ninput = 0
    ntaken = 0
    nmissed = 0

    # collect data and save it in arrays for all pairs of species.
    if options.loglevel >= 1:
        options.stdlog.write("# reading data (|=1000000, .=100000): ")
        sys.stdout.flush()

    for line in sys.stdin:

        if line[0] == "#":
            continue

        ninput += 1

        if ninput == 1 and not options.notitles:
            continue

        if options.loglevel >= 1:
            if not (ninput % 1000000):
                options.stdlog.write("|")
                options.stdlog.flush()
            elif not (ninput % 100000):
                options.stdlog.write(".")
                options.stdlog.flush()

        data = line[:-1].split("\t")
        try:
            species1 = rx_species.search(data[0]).groups()[0]
            species2 = rx_species.search(data[1]).groups()[0]
        except AttributeError:
            raise ValueError, "could not extract species from line: %s" % line[
                :-1]

        if options.map:
            if data[0] not in map_id2component or data[1] not in map_id2component:
                nmissed += 1
                continue

            component1 = map_id2component[data[0]]
            component2 = map_id2component[data[1]]
            same = component1 == component2
            if not same and not options.between_clusters:
                nskipped += 1
                continue
            elif same and options.between_clusters:
                nskipped += 1
                continue

            if component1 < component2:
                component1, component2 = component2, component1

            component = "%s-%s" % (component1, component2)

        species[species1] = 1
        species[species2] = 1
        if species1 > species2:
            species1, species2 = species2, species1

        key = "%s-%s" % (species1, species2)

        if key not in arrays:
            arrays[key] = []

        if len(arrays[key]) > options.max_pairs:
            x = random.randrange(len(arrays[key]))
            arrays[key][x] = ((component, float(data[options.columns[0]])))
        else:
            arrays[key].append((component, float(data[options.columns[0]])))

        ntaken += 1

    if options.loglevel >= 1:
        options.stdlog.write("done\n")
        options.stdlog.flush()

    # filter, if only best per component is to be taken
    if options.best_per_component:
        for k, vals in arrays.items():
            vals.sort()
            new_array = []
            last_c = None
            for c, v in vals:
                if last_c == c:
                    continue
                last_c = c
                new_array.append(v)

            arrays[k] = numpy.array(new_array)
    else:
        for k, vals in arrays.items():
            arrays[k] = numpy.array(map(lambda x: x[1], vals))

    # calculate histograms for all pairs of species
    r = numpy.arange(options.minbin, options.maxbin, options.binsize)
    numbins = len(r)

    for method in options.methods:

        if method == "hists":
            saved_hists = {}

            grand_total = scipy.stats.histogram2((), r)

            for s1 in species.keys():

                if options.loglevel >= 2:
                    options.stdlog.write(
                        "# processing histograms for species %s\n" % s1)
                    options.stdlog.flush()

                histograms = []
                titles = []
                total = scipy.stats.histogram2((), r)

                for s2 in species.keys():
                    if s1 < s2:
                        key = "%s-%s" % (s1, s2)
                    else:
                        key = "%s-%s" % (s2, s1)

                    if key not in arrays:
                        continue
                    if key in saved_hists:
                        h = saved_hists[key]
                    else:
                        h = scipy.stats.histogram2(arrays[key], r)
                        saved_hists[key] = h

                    for x in range(numbins):
                        total[x] += h[x]

                    titles.append(s2)
                    histograms.append(h)

                if options.only_total:
                    titles = []
                    histograms = []

                if options.only_total or options.add_total:
                    titles.append("total")
                    histograms.append(total)

                outfile = open(options.pattern % (s1, method), "w")
                outfile.write("bin\t%s\n" % "\t".join(titles))

                for x in range(numbins):
                    outfile.write("%f\t%s\n" % (
                        r[x], "\t".join([str(histograms[y][x]) for y in range(len(histograms))])))

                outfile.close()

                for x in range(numbins):
                    grand_total[x] += total[x]

            titles = ["total"]
            outfile = open(options.pattern % ("total", method), "w")
            outfile.write("bin\t%s\n" % "\t".join(titles))
            for x in range(numbins):
                outfile.write("%f\t%i\n" % (r[x], grand_total[x]))
            outfile.close()

        elif method == "stats":
            outfile = open(options.pattern % ("all", method), "w")
            outfile.write(
                "species1\tspecies2\tcounts\tmin\tmax\tmean\tmedian\tstd\n")

            for s1 in species.keys():
                histograms = []
                titles = []
                for s2 in species.keys():
                    if s2 < s1:
                        continue
                    key = "%s-%s" % (s1, s2)
                    if key not in arrays:
                        continue
                    v = arrays[key]
                    outfile.write("\t".join(("%s" % s1,
                                             "%s" % s2,
                                             "%i" % len(v),
                                             "%5.2f" % min(v),
                                             "%5.2f" % max(v),
                                             "%5.2f" % scipy.mean(v),
                                             "%5.2f" % scipy.median(v),
                                             "%5.2f" % numpy.std(v))) + "\n")
            v = numpy.concatenate(arrays.values())
            outfile.write("\t".join(("%s" % "all",
                                     "%s" % "all",
                                     "%i" % len(v),
                                     "%5.2f" % min(v),
                                     "%5.2f" % max(v),
                                     "%5.2f" % scipy.mean(v),
                                     "%5.2f" % scipy.median(v),
                                     "%5.2f" % numpy.std(v))) + "\n")

            outfile.close()

    options.stdlog.write("# ninput=%i, ntaken=%i, nskipped=%i, nmissed=%i\n" % (
        ninput, ntaken, nskipped, nmissed))
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
