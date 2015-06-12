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
optic/analyze_multiple_orthologs.py - 
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

   python optic/analyze_multiple_orthologs.py --help

Type::

   python optic/analyze_multiple_orthologs.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import pgdb

"""get ortholog counts for a set of predictions in a genome.
"""

parser = E.OptionParser(
    version="%prog version: $Id: optic/analyze_multiple_orthologs.py 2781 2009-09-10 11:33:14Z andreas $")


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser.add_option("-s", "--species", dest="species", type="string",
                      help="schema of master species.")

    parser.set_defaults(
        tablename_orthologs="orthology_pairwise1v5.orthologlinks_first",
        filename_ids="-",
        schemas=None,
        species=None,
    )

    (options, args) = E.Start(parser, add_database_options=True)

    dbhandle = pgdb.connect(options.psql_connection)

    if options.filename_ids == "-":
        ids, errors = IOTools.ReadList(sys.stdin)

    extra_options = ["schema1 = '%s'" % options.species,
                     "prediction_id1 IN ('%s')" % "','".join(ids)]

    if options.schemas:
        extra_options.append("schema2 IN ('%s')" % "','".join(options.schemas))

    statement = """SELECT prediction_id1, schema2, prediction_id2, gene_id2, gd1, gd2, td1, td2
    FROM %s
    WHERE schema1 != schema2 AND %s
    ORDER BY prediction_id1""" % (options.tablename_orthologs,
                                  " AND ".join(extra_options))

    cc = dbhandle.cursor()
    cc.execute(statement)
    result = cc.fetchall()
    cc.close()

    if options.schemas:
        schemas = options.schemas
    else:
        schemas = set(map(lambda x: x[1], result))

    # compute counts
    degeneracies = {}
    for x in ids:
        degeneracies[x] = {}
        for s in schemas:
            degeneracies[x][s] = (0, 0, 0, 0)

    for prediction_id1, schema2, prediction_id2, gene_id2, gd1, gd2, td1, td2 in result:
        degeneracies[prediction_id1][schema2] = (gd1, gd2, td1, td2)

    # output
    options.stdout.write("%s\t%s\n" % ("prediction_id", "\t".join(schemas)))
    for x in ids:
        options.stdout.write("%s" % x)
        for s in schemas:
            options.stdout.write("\t%s:%s:%s:%s" % degeneracies[x][s])
        options.stdout.write("\n")

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
