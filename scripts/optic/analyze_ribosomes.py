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
optic/analyze_ribosomes.py - 
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

   python optic/analyze_ribosomes.py --help

Type::

   python optic/analyze_ribosomes.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import CGAT.Experiment as E
import CGAT.CSV as CSV
import CGAT.IOTools as IOTools
import CGAT.Stats as Stats

USAGE = """ program $Id: optic/analyze_ribosomes.py 2781 2009-09-10 11:33:14Z andreas $

analyse ribosomal proteins.

Input is a file of orthologs, a list of schemas and a list
of sequence identifiers.
"""


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/analyze_ribosomes.py 2781 2009-09-10 11:33:14Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-s", "--schemas", dest="schemas", type="string",
                      help="schemas in the set.")

    parser.add_option("-e", "--field-extract", dest="field_extract", type="string",
                      help="pattern for the field to extract.")

    parser.add_option("-c", "--field-compare", dest="field_compare", type="string",
                      help="pattern for the field to compare.")

    parser.add_option("-i", "--identifiers-tsv-file", dest="filename_identifiers", type="string",
                      help="identifiers in the positive set.")

    parser.add_option("-u", "--filename-subset", dest="filename_subset", type="string",
                      help="subset in the positive set.")

    parser.add_option("--filter-min-ratio", dest="filter_min_ratio", type="float",
                      help="minimum boundary for filter.")

    parser.add_option("--filter-max-ratio", dest="filter_max_ratio", type="float",
                      help="maximum boundary for filter.")

    parser.add_option("-o", "--output-fields", dest="output_fields", type="string",
                      help="output fields, choices are: zscore, val, nvals, sum, min, max, stddev, mean, median.")

    parser.add_option("--output-filename-pattern", dest="output_pattern", type="string",
                      help="pattern for table headers, should contain %s for schema and %s for field anme.")

    parser.add_option("-f", "--output-format", dest="output_format", type="choice",
                      choices=("table", "list", "values"),
                      help="output format. Tabular form (one row per ortholog) or list form.")

    parser.add_option("--format", dest="format", type="string",
                      help="output format for numbers.")

    parser.add_option("--remove-na", dest="remove_na", action="store_true",
                      help="remove entries with any na values.")

    parser.set_defaults(
        field_extract="%s_length",
        field_compare="%s_length",
        filename_identifiers=None,
        filename_subset=None,
        filter_min_ratio=0.00,
        filter_max_ratio=0.00,
        schemas="",
        output_fields="",
        output_pattern="%s_%s",
        output_format="table",
        format="%6.4f",
        remove_na=False,
    )

    (options, args) = E.Start(parser, add_csv_options=True)

    options.schemas = options.schemas.split(",")
    if not options.schemas:
        raise "please supply schemas."

    if options.output_fields:
        options.output_fields = options.output_fields.split(",")
    else:
        options.output_fields = ()

    fields, table = CSV.ReadTable(sys.stdin)

    map_fields2column = {}
    for x in fields:
        map_fields2column[x] = len(map_fields2column)

    if options.loglevel >= 1:
        options.stdlog.write("# read a %i x %i table.\n" %
                             (len(table), len(fields)))

    if options.filename_subset:
        subset, nerrors = IOTools.ReadList(open(options.filename_subset, "r"))
        subset = set(subset)

        table = filter(lambda x: x[0] in subset, table)

        if options.loglevel >= 1:
            options.stdlog.write("# subset of %i entries reduced table to a %i x %i table.\n" % (
                len(subset), len(table), len(fields)))

    if options.filename_identifiers:
        identifiers, nerrors = IOTools.ReadList(
            open(options.filename_identifiers, "r"))
    else:
        identifiers = []

    identifiers = set(identifiers)

    # extract rows with positive identifiers
    positive_rows = filter(lambda x: x[0] in identifiers, table)

    if options.loglevel >= 1:
        options.stdlog.write("# subset of %i identifiers gives %i positive entries.\n" % (
            len(identifiers), len(positive_rows)))

    if options.output_format == "table":
        options.stdout.write("id")
        for schema in options.schemas:
            if options.output_fields:
                for field in options.output_fields:
                    options.stdout.write(
                        "\t" + options.output_pattern % (schema, field))
            else:
                options.stdout.write("\t%s" % (schema))

        options.stdout.write("\n")
    else:
        options.stdout.write("schema\tvalue\n")

    if identifiers:
        for row in positive_rows:

            if options.output_format == "table":
                options.stdout.write(row[0])

            for schema in options.schemas:

                # set fields for extraction
                f_extract = map_fields2column[options.field_extract % schema]
                f_compare = map_fields2column[options.field_compare % schema]

                # get region for extraction
                if row[f_compare] != "na":
                    r = float(row[f_compare])
                    if options.filter_min_ratio or options.filter_max_ratio:
                        mi = r * options.filter_min_ratio
                        ma = r * options.filter_max_ratio
                        f = lambda x: x[f_compare] != "na" and float(x[f_compare]) >= mi and float(
                            x[f_compare]) <= ma and x[0] not in identifiers and x[f_extract] != "na"
                    else:
                        f = lambda x: x[0] not in identifiers and x[
                            f_extract] != "na"
                    # extract values: filter by minimum and maximum range and remove
                    # positive identifiers.
                    v = float(row[f_extract])
                    values = map(
                        lambda x: float(x[f_extract]), filter(f, table))

                    stats = Stats.DistributionalParameters(values)
                else:
                    v = None

                for field in options.output_fields:

                    if v is not None:
                        if field == "zscore":
                            f = options.format % stats.getZScore(v)
                        elif field == "diff":
                            f = options.format % (v - stats["mean"])
                        elif field == "reldiff":
                            f = options.format % (
                                (v - stats["mean"]) / stats["mean"])
                        elif field == "val":
                            f = options.format % v
                        else:
                            f = options.format % stats[field]
                    else:
                        f = "na"

                    if options.output_format == "table":
                        options.stdout.write("\t%s" % f)
                    elif options.output_format == "list":
                        options.stdout.write("%s\t%s\n" % (schema, f))
                    elif options.output_format == "values":
                        options.stdout.write("%s\t%s\t%5.2f\t%s\n" % (
                            row[0], schema, v, ",".join(map(lambda x: options.format % x, values))))

            if options.output_format == "table":
                options.stdout.write("\n")

    else:

        extract_columns = []

        for schema in options.schemas:
            extract_columns.append(
                map_fields2column[options.field_extract % schema])

        # simply dump a subset of values
        for row in table:

            skip = False

            if options.filter_min_ratio or options.filter_max_ratio:

                master = options.schemas[0]

                v = row[map_fields2column[options.field_compare % master]]

                if v == "na":
                    continue

                v = float(v)

                mi = v * options.filter_min_ratio
                ma = v * options.filter_max_ratio

                for schema in options.schemas[1:]:

                    r = row[map_fields2column[options.field_compare % schema]]

                    if r == "na":
                        if options.remove_na:
                            skip = True
                        continue

                    r = float(r)

                    if r < mi or r > ma:
                        skip = True
                        if options.loglevel >= 3:
                            if options.format == "table":
                                options.stdout.write("* ")
                                options.stdout.write("%s\t" % row[0])
                                options.stdout.write(
                                    "\t".join([row[y] for y in extract_columns]))
                                options.stdout.write("\n")
                        break

            if skip:
                continue

            if options.output_format == "table":
                options.stdout.write("%s\t" % row[0])
                options.stdout.write(
                    "\t".join([row[y] for y in extract_columns]))
                options.stdout.write("\n")

            elif options.output_format == "list":
                has_na = False
                for x in range(len(options.schemas)):
                    v = row[extract_columns[x]]
                    if v == "na":
                        has_na = True

                if has_na and options.remove_na:
                    continue

                for x in range(len(options.schemas)):
                    options.stdout.write(
                        "%s\t%s\n" % (options.schemas[x], row[extract_columns[x]]))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
