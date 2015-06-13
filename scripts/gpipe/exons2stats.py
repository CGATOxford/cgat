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
gpipe/exons2stats.py - 
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

   python gpipe/exons2stats.py --help

Type::

   python gpipe/exons2stats.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import csv
import pgdb
import CGAT.Experiment as E
import CGAT.Exons as Exons

parser = E.OptionParser(
    version="%prog version: $Id: gpipe/exons2stats.py 2781 2009-09-10 11:33:14Z andreas $")


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser.add_option("-q", "--quality", dest="quality", type="string",
                      help="quality categories to take into account.")
    parser.add_option("-f", "--format=", dest="format", type="string",
                      help="input format [exons|gff|table]")

    parser.add_option("-e", "--exons-file=", dest="tablename_exons", type="string",
                      help="table name with exons.")
    parser.add_option("-p", "--predictions=", dest="tablename_predictions", type="string",
                      help="table name with predictions.")
    parser.add_option("-n", "--non-redundant", dest="non_redundant", action="store_true",
                      help="only non-redundant predictions.")
    parser.add_option("-s", "--schema", dest="schema", type="string",
                      help="schema to use.")

    parser.set_defaults(
        fields=["Id", "NumExons", "GeneLength", "MinExonLength",
                "MaxExonLength", "MinIntronLength", "MaxIntronLength"],
        tablename_exons="exons",
        tablename_predictions="predictions",
        quality=None,
        non_redundant=False,
        schema=None,
        tablename_redundant="redundant",
        tablename_quality="quality",
        format="exons",
    )

    (options, args) = E.Start(
        parser, add_csv_options=True, add_database_options=True)

    if options.quality:
        options.quality = options.quality.split(",")

    if options.format == "table":
        dbhandle = pgdb.connect(options.psql_connection)
        exons = Exons.GetExonBoundariesFromTable(dbhandle,
                                                 options.tablename_predictions,
                                                 options.tablename_exons,
                                                 non_redundant_filter=options.non_redundant,
                                                 quality_filter=options.quality,
                                                 table_name_quality=options.tablename_quality,
                                                 table_name_redundant=options.tablename_redundant,
                                                 schema=options.schema)
    else:
        exons = Exons.ReadExonBoundaries(sys.stdin)

    stats = Exons.CalculateStats(exons)

    print "\t".join(options.fields)

    writer = csv.DictWriter(sys.stdout,
                            options.fields,
                            dialect=options.csv_dialect,
                            lineterminator=options.csv_lineterminator,
                            extrasaction='ignore')

    for k, v in stats.items():
        v["Id"] = k
        writer.writerow(v)

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
