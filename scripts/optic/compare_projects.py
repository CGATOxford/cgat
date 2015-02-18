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
optic/compare_projects.py -
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

   python optic/compare_projects.py --help

Type::

   python optic/compare_projects.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import string
import getopt
import CGAT.Experiment as E
import pgdb
import numpy

USAGE = """python %s [OPTIONS] schema1 schema2 [...]

Version: $Id: optic/compare_projects.py 2387 2009-01-07 16:33:07Z andreas $

Dump out various overview for comparisions between genomes

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-C, --connection=               postgres connection string
-r, --output-report=                   type of report to create
-f, --fields=                   fields to include in report
-a, --associated-schemas        fields will be grouped by schema
-s, --separator=                separator to use [default=tab]
-u, --summary=                  calculate a summary report
""" % sys.argv[0]

param_loglevel = 1

param_long_options = ["verbose=", "help", "connection=",
                      "report=", "fields=",
                      "associate-schemas",
                      "separator=", "summary=",
                      "version"]

param_short_options = "v:hC:r:f:as:"
param_schemas = None
param_report = None
param_fields = None

param_connection = "db:andreas"

param_associate_fields = 1

param_separator = "\t"

param_summary = None


def GetExpandedList(l, s):
    """return a list where every elment in l is duplicated s times.
    """
    nl = []
    for x in l:
        nl += [x] * s
    return nl


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    try:
        optlist, args = getopt.getopt(
            sys.argv[1:], param_short_options, param_long_options)
    except getopt.error, msg:
        print USAGE, msg
        sys.exit(2)

    for o, a in optlist:
        if o in ("-v", "--verbose"):
            param_loglevel = int(a)
        elif o in ("--version", ):
            print "version="
            sys.exit(0)
        elif o in ("-h", "--help"):
            print USAGE
            sys.exit(0)
        elif o in ("-C", "--connection"):
            param_connection = a
        elif o in ("-r", "--output-report"):
            param_report = a
        elif o in ("-f", "--fields"):
            param_fields = string.split(a, ",")
        elif o in ("-a", "--associate-schemas"):
            param_associate_fields = 0
        elif o in ("-s", "--separator"):
            param_separator = a
        elif o in ("-u", "--summary"):
            param_summary = string.split(a, ",")

    if len(args) < 2:
        print USAGE
        print "please supply at least two schemas to compare."
        sys.exit(1)
    else:
        param_schemas = args

    print E.GetHeader()
    print E.GetParams()

    dbhandle = pgdb.connect(param_connection)
    headers1 = []
    headers2 = []

    ##########################################################################
    # setup the main statement
    if param_report == "queries":
        table_name = "queries"
        statement = """SELECT %s.%s.query_token, %%s
        FROM %s
        WHERE %s
        """ % (param_schemas[0], table_name,
               string.join(map(lambda x: "%s.%s" % (x, table_name),
                               param_schemas), ", "),
               string.join(map(lambda x: "%s.%s.query_token = %s.%s.query_token" %
                               (param_schemas[0], table_name, x, table_name),
                               param_schemas[1:]), " AND "))
        headers1 = [""]
        headers2 = ["queries"]

    elif param_report == "domains":
        table_name = "domains_summary"
        statement = """SELECT %s.%s.domain_id, %%s
        FROM %s
        WHERE %s
        """ % (param_schemas[0], table_name,
               string.join(map(lambda x: "%s.%s" % (x, table_name),
                               param_schemas), ", "),
               string.join(map(lambda x: "%s.%s.domain_id = %s.%s.domain_id" %
                               (param_schemas[0], table_name, x, table_name),
                               param_schemas[1:]), " AND "))
        headers1 = [""]
        headers2 = ["domains"]

    ##########################################################################
    # associated_fields == 1, if fields shall be grouped together
    # build the SQL Statement
    # summary_fields: starting coordinates of field separators for summary
    # calculation
    if param_associate_fields:
        fields = string.join(map(lambda x, y: "%s.%s.%s" % (x, table_name, y),
                                 param_schemas * len(param_fields),
                                 GetExpandedList(param_fields, len(param_schemas))), ",")
        summary_fields = range(
            1, len(param_schemas) * (len(param_fields) + 1), len(param_schemas))
        if param_summary:
            for x in param_summary:
                headers1 += [x] + [""] * (len(param_fields) - 1)
                headers2 += param_fields
        else:
            for x in param_fields:
                headers1 += [x] + [""] * (len(param_schemas) - 1)
                for y in param_schemas:
                    headers2 += [y]

    else:
        fields = string.join(map(lambda x, y: "%s.%s.%s" % (x, table_name, y),
                                 GetExpandedList(
                                     param_schemas, len(param_fields)),
                                 param_fields * len(param_schemas)
                                 ), ",")
        summary_fields = range(
            1, (len(param_schemas) + 1) * len(param_fields), len(param_fields))
        if param_summary:
            for x in param_schemas:
                headers1 += [x] + [""] * (len(param_summary) - 1)
                headers2 += param_summary
        else:
            for x in param_schemas:
                headers1 += [x] + [""] * (len(param_fields) - 1)
                for y in param_fields:
                    headers2 += [y]

    if param_loglevel >= 2:
        print "#", statement % fields

    cc = dbhandle.cursor()
    try:
        cc.execute(statement % fields)
        result = cc.fetchall()
    except pgdb.DatabaseError, msg:
        print "# query failed with message", msg
        result = None

    print string.join(headers1, param_separator)
    print string.join(headers2, param_separator)

    if param_summary:
        funcs = []
        for method in param_summary:
            if method == "min":
                funcs.append(min)
            elif method == "max":
                funcs.append(max)
            elif method == "sum":
                funcs.append(sum)
            elif method == "mean":
                funcs.append(numpy.mean)
            elif method == "samplestd":
                funcs.append(numpy.std)
            elif method == "count":
                funcs.append(lambda x: len(filter(lambda y: y > 0, x)))
            else:
                raise "unknown method"

        for r in result:
            sys.stdout.write(r[0])

            for f in funcs:
                for x in range(0, len(summary_fields) - 1):
                    sys.stdout.write(
                        "\t" + str(f(r[summary_fields[x]:summary_fields[x + 1]])))
            sys.stdout.write("\n")
    else:
        for r in result:
            print string.join(map(str, r), param_separator)

    print E.GetFooter()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
