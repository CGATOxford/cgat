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
gpipe/get_predictions.py - 
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

   python gpipe/get_predictions.py --help

Type::

   python gpipe/get_predictions.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import pgdb
import CGAT.Experiment as E

""" program $Id: gpipe/get_predictions.py 2781 2009-09-10 11:33:14Z andreas $

gpipe/get_predictions.py

retrieve predictions/exons

"""


def GetResult(dbhandle, keys, options):

    if "exon_from" in keys:

        statement = """
        SELECT *
        FROM %(schema)s.exons AS e,
        %(schema)s.predictions AS p
        WHERE p.prediction_id = e.prediction_id AND
        p.sbjct_token = '%(sbjct_token)s' AND
        p.sbjct_strand = '%(sbjct_strand)s' AND
        (NOT 
        ( e.genome_exon_to < '%(exon_from)i' OR e.genome_exon_from > '%(exon_to)i' ) )
        """ % keys
    else:
        raise "don't know what to do."

    cc = dbhandle.cursor()
    cc.execute(statement)
    results = cc.fetchall()
    cc.close()

    return results


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gpipe/get_predictions.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-s", "--schema", dest="schema", type="string",
                      help="database schema to use.")

    parser.set_defaults(
        schema=None,
        method="ranges",
        result="exons",
    )

    (options, args) = E.Start(parser, add_database_options=True)

    dbhandle = pgdb.connect(options.psql_connection)

    keys = {'schema': options.schema}

    ntested, nmissed = 0, 0
    if options.method == "ranges":

        lines = filter(lambda x: x[0] != "#", sys.stdin.readlines())

        for line in lines:
            sbjct_token, sbjct_strand, range_from, range_to = line.split(
                "\t")[:4]
            range_from, range_to = map(int, (range_from, range_to))

            keys["sbjct_token"] = sbjct_token
            keys["sbjct_strand"] = sbjct_strand
            keys["exon_from"] = range_from
            keys["exon_to"] = range_to

            ntested += 1
            result = GetResult(dbhandle, keys, options)

            if len(result) == 0:
                nmissed += 1
                continue

            print "########################"
            print "# matches for:", line[:-1]
            for r in result:
                print "\t".join(map(str, r))

    print "# ntested=%i, nmissed=%i" % (ntested, nmissed)
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
