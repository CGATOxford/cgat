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
optic/analyze_duplications.py -
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Analyse duplications in various species

The file expects a file of duplication events:

species1 species2 location function height

is for a duplication of species2 in species.

Usage
-----

Example::

   python optic/analyze_duplications.py --help

Type::

   python optic/analyze_duplications.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import CGAT.Experiment as E
import pgdb
import numpy
import CGAT.MatlabTools as MatlabTools
import CGAT.Histogram as Histogram
import CGAT.TreeTools as TreeTools
import scipy
import scipy.optimize


def GetSubset(input_data, l, f):
    """get subset of data based on locations/functions."""

    data = []
    for x in input_data:

        location, function = x[3], x[4]

        if l != "all":
            if l == "nojunk":
                if location not in ("local", "muller", "trans", "nonlocal"):
                    continue
            elif l == "cis":
                if location not in ("local", "nonlocal"):
                    continue
            else:
                if location != l:
                    continue

        if f != "all" and \
                function != f:
            continue

        data.append(x)

    return data


def residuals_decay(p, y, x):
    """decay function

    y = A0 + A1 * exp( k * x )
    """
    A0, A1 = p

    err = y - (A0 * numpy.exp(A1 * x))

    return err


def fit(histogram, parameters, fresiduals=residuals_decay):
    """fit a function to a histogram.

    (function of the form x = a + b * exp(-c x )
    """
    in_x = numpy.array(map(lambda x: x[0], histogram))
    in_y = numpy.array(map(lambda x: x[1], histogram))

    try:
        result = scipy.optimize.leastsq(
            fresiduals, parameters, args=(in_y, in_x),
            full_output=True)
    except scipy.optimize.minpack.error:
        return None

    params, info, ier, msg = result

    if ier != 1:
        return None

    return params


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/analyze_duplications.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-s", "--species", dest="species", type="string",
                      help="species to use.")

    parser.add_option("-p", "--column-prefix", dest="prefix", type="string",
                      help="prefix to use for temporary files.")

    parser.add_option("-m", "--method", dest="method", type="string",
                      help="method to use [counts|lists|hists|links].")

    parser.add_option("-o", "--filename-output", dest="filename_output", type="string",
                      help="output filename.")

    parser.add_option("-f", "--functions", dest="functions", type="string",
                      help="functions to grep [functional|pseudo|all].")

    parser.add_option("-l", "--locations", dest="locations", type="string",
                      help="locations to grep [local|nojunk|all|...].")

    parser.add_option("-b", "--bin-size", dest="bin_size", type="string",
                      help="bin size.")

    parser.add_option("-i", "--fit", dest="fit", type="string",
                      help="fitting method [decay|power]")

    parser.add_option("--min-value", dest="min_value", type="float",
                      help="minimum value for histogram.")

    parser.add_option("--max-value", dest="max_value", type="float",
                      help="maximum value for histogram.")

    parser.add_option("--use-relative-height", dest="use_relative_height", action="store_true",
                      help="use relative height values.")

    parser.add_option("--reverse", dest="reverse", action="store_true",
                      help="""reverse species. Histograms will show the age of duplications for
                      duplicates in other genomes.""")

    parser.set_defaults(species="",
                        functions="functional,pseudo,all",
                        locations="local,nojunk,all",
                        filename_output=None,
                        bin_size=1.0,
                        min_value=None,
                        max_value=None,
                        nonnull=None,
                        use_relative_height=False,
                        header=True,
                        fit=None,
                        reverse=False,
                        method="counts")

    (options, args) = E.Start(parser, add_database_options=True)

    options.species = options.species.split(",")
    options.locations = options.locations.split(",")
    options.functions = options.functions.split(",")

    if len(options.species) == 0:
        raise "please supply list of species."

    dbhandle = pgdb.connect(options.psql_connection)

    input_data = map(lambda x: x[:-1].split("\t"),
                     filter(lambda x: x[0] != "#", sys.stdin.readlines()))

    # remove header
    if options.header:
        del input_data[0]

    # decide which columns to take
    # 1st column: species1: this is the species in which duplications have occured.
    # 2nd column: species2: this is the species with respect to which duplications occured.
    # 3rd column: clusterid
    # 4th column: chromosomes
    # 5th column: function
    # 6th column: height
    # 7th column: relative height
    # 8th column: locations
    # 9th column: tree
    if options.use_relative_height:
        take = (0, 1, 2, 3, 4, 6, 7, 8)
    else:
        take = (0, 1, 2, 3, 4, 5, 7, 8)

    for x in range(len(input_data)):
        input_data[x] = tuple([input_data[x][y] for y in take])

    map_pos2species = []
    map_species2pos = {}
    for x in range(len(options.species)):
        map_species2pos[options.species[x]] = x
        map_pos2species.append(options.species[x])

    outfile = None

    if options.method in ("counts", "medians"):

        if options.method == "counts":
            func = len
        elif options.method == "medians":
            func = numpy.median

        for location in options.locations:

            for function in options.functions:
                matrix = numpy.zeros(
                    (len(options.species), len(options.species)), numpy.Float)

                data = GetSubset(input_data, location, function)

                # sort by species1 and species2
                data.sort()

                last_species1, last_species2 = None, None
                values = []
                for species1, species2, cluster_id, l, f, height, locations, tree in data:

                    if last_species1 != species1 or last_species2 != species2:

                        if len(values) > 0:
                            matrix[map_species2pos[last_species1], map_species2pos[last_species2]] = func(
                                values)

                        values = []
                        last_species1 = species1
                        last_species2 = species2

                    values.append(float(height))

                if len(values) > 0:
                    matrix[map_species2pos[last_species1],
                           map_species2pos[last_species2]] = func(values)

                if options.filename_output:
                    dict = {"f": function,  "l": location}
                    outfile = open(options.filename_output % dict, "w")
                else:
                    outfile = sys.stdout
                    outfile.write("matrix for method %s: location: %s, function: %s\n" % (
                        options.method, location, function))

                if options.method == "medians":
                    format = "%6.4f"
                elif options.method == "counts":
                    format = "%i"
                MatlabTools.WriteMatrix(matrix,
                                        outfile=outfile,
                                        format=format,
                                        row_headers=options.species,
                                        col_headers=options.species)

                if options.filename_output:
                    outfile.close()

    elif options.method in ("lists", "lists-union"):
        # write lists of duplicated genes in species1 as compared to species2
        # according to location/function
        # First field : gene name
        # Second field: cluster id
        # Third field : number of other genes in cluster
        # Fourth field: location of gene
        written = {}
        for location in options.locations:

            for function in options.functions:

                values = [[[] for y in range(len(options.species))]
                          for x in range(len(options.species))]

                data = GetSubset(input_data, location, function)

                # sort by species1 and species2
                data.sort()

                last_species1, last_species2 = None, None

                for species1, species2, cluster_id, l, f, height, locations, tree in data:

                    if last_species1 != species1 or last_species2 != species2:

                        # write trees per cluster
                        if options.filename_output:
                            if options.method == "lists":
                                if outfile:
                                    outfile.close()
                                dict = {
                                    "f": function,  "l": location, "s": species1, "o": species2}
                                written = {}
                                outfile = open(
                                    options.filename_output % dict, "w")
                            elif options.method == "lists-union":
                                if last_species1 != species1:
                                    if outfile:
                                        outfile.close()
                                    dict = {
                                        "f": function,  "l": location, "s": species1}
                                    written = {}
                                    outfile = open(
                                        options.filename_output % dict, "w")
                        else:
                            outfile = sys.stdout
                            if options.method == "lists":
                                outfile.write("location: %s, function: %s, species1: %s, species2: %s\n" % (
                                    location, function, species1, species2))
                                written = {}
                            elif options.method == "lists-union":
                                if last_species1 != species1:
                                    outfile.write(
                                        "location: %s, function: %s, species1: %s\n" % (location, function, species1))
                                    written = {}

                        last_species1 = species1
                        last_species2 = species2

                    # get tree
                    tt = TreeTools.Newick2Tree(tree)
                    taxa = TreeTools.GetTaxa(tt)
                    for t in taxa:
                        if t in written:
                            continue
                        outfile.write("%s\t%s\t%i\n" %
                                      (t, cluster_id, len(taxa)))
                        written[t] = 1

    elif options.method in ("hists", "fit-decay"):

        for location in options.locations:

            for function in options.functions:

                values = [[[] for y in range(len(options.species))]
                          for x in range(len(options.species))]

                data = GetSubset(input_data, location, function)

                data.sort()

                ###############################################################
                # convert to matrix of list
                # values[x][y] contains heights of duplications in species x
                # with reference to y

                for species1, species2, cluster_id, l, f, height, locations, tree in data:
                    try:
                        values[map_species2pos[species1]][
                            map_species2pos[species2]].append(float(height))
                    except KeyError:
                        continue

                ###############################################################
                ###############################################################
                ###############################################################
                # calculate histograms per species
                ###############################################################
                for s in options.species:
                    histograms = []
                    headers = []

                    if options.filename_output:
                        dict = {"f": function,  "l": location, "s": s}
                        outfile = open(options.filename_output % dict, "w")
                    else:
                        outfile = sys.stdout
                        outfile.write(
                            "location: %s, function: %s\n" % (location, function))

                    for x in range(len(options.species)):

                        if options.reverse:
                            # duplications in species x
                            vv = values[x][map_species2pos[s]]
                        else:
                            # duplications in species s
                            vv = values[map_species2pos[s]][x]

                        if len(vv) == 0:
                            pass
                        else:
                            headers.append(options.species[x])
                            h = Histogram.Calculate(vv,
                                                    increment=options.bin_size,
                                                    min_value=options.min_value,
                                                    max_value=options.max_value,
                                                    no_empty_bins=True)

                            if options.method == "fit-decay":
                                result = fit(h, [2.0, -1.0])
                                if result:
                                    outfile.write("%s\t%s\t%s\t%i\t%f\t%f\ty = %f * exp ( %f * x )\n" % (
                                        "function",
                                        s,
                                        options.species[x],
                                        h[0][1],
                                        result[0],
                                        result[1],
                                        result[0],
                                        result[1],
                                    ))
                            elif options.method == "hists":
                                histograms.append(h)

                    if options.method == "hists":
                        combined_histogram = Histogram.Combine(
                            histograms, missing_value="-")

                        outfile.write("bin\t" + "\t".join(headers) + "\n")
                        Histogram.Write(outfile, combined_histogram)

                    if options.filename_output:
                        outfile.close()
                    else:
                        outfile.flush()

    elif options.method == "pairs":

        # get branches with 0 branchlength

        for location in options.locations:

            if options.loglevel >= 2:
                options.stdlog.write("# processing location %s\n" % location)

            for function in options.functions:

                if options.loglevel >= 2:
                    options.stdlog.write(
                        "#   processing function %s " % function)
                    options.stdlog.flush()

                data = GetSubset(input_data, location, function)

                if options.loglevel >= 2:
                    options.stdlog.write("with %i data points\n" % len(data))
                    options.stdlog.flush()

                data.sort()
                last_species1, last_species2, last_cluster_id = None, None, None

                values = []
                for species1, species2, cluster_id, l, f, height, locations, tree in data:

                    if last_species1 != species1 or last_species2 != species2:

                        # write trees per cluster
                        if options.filename_output:
                            if outfile:
                                outfile.close()
                            dict = {
                                "f": function,  "l": location, "s": species1, "o": species2}
                            outfile = open(options.filename_output % dict, "w")
                        else:
                            outfile = sys.stdout
                            outfile.write("location: %s, function: %s, species1: %s, species2: %s\n" % (
                                location, function, species1, species2))

                        last_species1 = species1
                        last_species2 = species2
                        last_cluster_id = None

                    if last_cluster_id != cluster_id:
                        if last_cluster_id is not None:
                            pass

                        last_cluster_id = cluster_id

                    outfile.write("%s\t%s\t%s\t%s\n" %
                                  (cluster_id, height, locations, tree))

    elif options.method == "links":

        # write a tree for each species pair:
        # each node is a gene+location, the weight of the vertex is the height
        # further info added: cluster_id for the duplication

        for location in options.locations:

            if options.loglevel >= 2:
                options.stdlog.write("# processing location %s\n" % location)

            for function in options.functions:

                if options.loglevel >= 2:
                    options.stdlog.write(
                        "#   processing function %s " % function)
                    options.stdlog.flush()

                data = GetSubset(input_data, location, function)

                if options.loglevel >= 2:
                    options.stdlog.write("with %i data points\n" % len(data))
                    options.stdlog.flush()

                # stores duplications within first species as compared to
                # second species
                values = [[[] for y in range(len(options.species))]
                          for x in range(len(options.species))]

                for species1, species2, cluster_id, l, f, height, locations, tree in data:
                    values[map_species2pos[species1]][map_species2pos[species2]].append(
                        (cluster_id, -len(locations), locations, tree))

                # get links per species
                for s in options.species:
                    if options.loglevel >= 2:
                        options.stdlog.write(
                            "#     processing species %s\n" % s)

                    headers = []
                    for x in range(len(options.species)):

                        if map_pos2species[x] == s:
                            continue

                        vv = values[map_species2pos[s]][x]
                        vv.sort()

                        # write trees per cluster
                        if options.filename_output:
                            dict = {
                                "f": function,  "l": location, "s": s, "o": map_pos2species[x]}
                            outfile = open(options.filename_output % dict, "w")
                        else:
                            outfile = sys.stdout
                            outfile.write("location: %s, function: %s, species1: %s, species2: %s\n" % (
                                location, function, s, map_pos2species[x]))

                        # only print out largest tree
                        last_cluster_id = None
                        for cluster_id, n, locations, tree in vv:
                            if cluster_id != last_cluster_id:
                                outfile.write("%s\t%s\t%s\n" %
                                              (cluster_id, locations, tree))
                                last_cluster_id = cluster_id

                        if options.filename_output:
                            outfile.close()

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
