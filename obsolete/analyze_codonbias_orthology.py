'''
analyze_codonbias_orthology.py - 
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

   python analyze_codonbias_orthology.py --help

Type::

   python analyze_codonbias_orthology.py --help

for command line help.

Command line options
--------------------

'''
import sys
import psycopg2

"""analyze codonbias of orthologs.

The input looks like:

schema1 gene1 cai1 schema2, gene2, cai2
"""

import CGAT.Experiment as E


def GetThreshold(dbhandle, schema):
    """get dominant set threshold for schema."""

    statement = "SELECT cai FROM %s.codonbias WHERE is_selected = TRUE ORDER BY cai DESC"

    # get thresholds
    cc = dbhandle.cursor()
    cc.execute(statement)
    r = cc.fetchall()
    cc.close()

    return threshold

# ------------------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: analyze_codonbias_orthology.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("--filters", dest="filters", type="string",
                      help="Filters to use for filtering sequences [all|codon1|codon2|codon3|44].")
    parser.add_option("--fields", dest="fields", type="string",
                      help="Fields to output [aligned|nunaligned1|nunaligned2|identical|transitions|transversions|jc69|t92].")

    parser.set_defaults(
        filename_map=None,
        gap_char="-",
    )

    (options, args) = E.Start(parser, add_pipe_options=True,
                              add_database_options=True)

    data = map(lambda x: x[:-1].split("\t"),
               filter(lambda x: x[0] != "#", sys.stdin.readlines()))

    schema1 = data[0][0]
    schema2 = data[0][3]

    dbhandle = psycopg2.connect(options.psql_connection)

    threshold1 = GetThreshold(dbhandle, schema1)
    threshold2 = GetThreshold(dbhandle, schema2)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
