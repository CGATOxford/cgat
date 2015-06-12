'''
analyze_go.py - 
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

   python analyze_go.py --help

Type::

   python analyze_go.py --help

for command line help.

Command line options
--------------------

'''
import os
import sys
import string
import re
import tempfile
import subprocess
import optparse

import CGAT.Experiment as E
import CGAT.Database as Database


def WriteBackground(go_type, options, suffix):

    statement = """SELECT DISTINCTROW
    gene_stable_id, glook_%s_id, description, olook_evidence_code
    FROM %s.%s_go_%s__go_%s__main
    WHERE glook_%s_id IS NOT NULL
    GROUP BY gene_stable_id, glook_%s_id, description
    ORDER BY gene_stable_id
    """ % (go_type,
           options.database, options.species, go_type, go_type,
           go_type, go_type)

    result = dbhandle.Execute(statement).fetchall()

    outfile = open("%s%s" % (options.prefix, suffix), "w")

    for r in result:
        outfile.write(" : ".join(map(str, r)) + "\n")

    outfile.close()


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: analyze_go.py 309 2005-12-01 15:50:26Z andreas $")

    dbhandle = Database.Database()

    parser.add_option("-s", "--species", dest="species", type="string",
                      help="species to use.")

    parser.add_option("-p", "--column-prefix", dest="prefix", type="string",
                      help="prefix to use for temporary files.")

    parser.set_defaults(species="dmelanogaster")
    parser.set_defaults(database="ensembl_mart_31")
    parser.set_defaults(prefix="dm_go_")

    (options, args) = E.Start(parser, add_database_options=True)

    dbhandle.Connect(options)

    WriteBackground("biol_process", options, "bp")
    WriteBackground("cell_location", options, "lm")
    WriteBackground("mol_function", options, "fm")


if __name__ == "__main__":
    sys.exit(main(sys.argv))
