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
optic/analyze_genes.py - 
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

   python optic/analyze_genes.py --help

Type::

   python optic/analyze_genes.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import CGAT.Experiment as E
import pgdb

""" program $Id: optic/analyze_genes.py 2781 2009-09-10 11:33:14Z andreas $

count number of matches found in gene prediction pipeline.

Count by category and by gene.
"""

parser = E.OptionParser(
    version="%prog version: $Id: optic/analyze_genes.py 2781 2009-09-10 11:33:14Z andreas $")


def GetCounts(r, priority):
    """count transcript by genes sorted by priority.
    """
    last_g = None

    qq = {}
    dd = {}

    for x in priority:
        dd[x] = 0

    for g, q, p in r:
        if last_g != g:
            if last_g:
                for x in priority:
                    if x in qq:
                        dd[x] += 1
                        break
            qq = {}
            last_g = g

        qq[q] = 1

    for x in priority:
        if x in qq:
            dd[x] += 1
            break

    return dd


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser.add_option("-g", "--genomes", dest="genomes", type="string",
                      help="genomes to analyse.")
    parser.add_option("-i", "--priority", dest="priority", type="string",
                      help="quality priority.")
    parser.add_option("-s", "--method=sort --sort-order", dest="sort", type="string",
                      help="sort order.")
    parser.add_option("-t", "--table-orthology", dest="table_orthology", type="string",
                      help="tablename with orthology relationships.")
    parser.add_option("-e", "--schema", dest="orthology_reference_schema", type="string",
                      help="schema name of reference.")
    parser.add_option("-f", "--method=filter --filter-method", dest="filename_filter", type="string",
                      help="filename with schema|prediction_id|gene to use as filter. The prediction_ids are used for filtering.")

    parser.set_defaults(
        genomes="",
        priority="CG,PG,SG,RG,CP,PP,SP,RP,CF,PF,SF,UG,UP,UF,BF,UK",
        sort="CG,PG,SG,RG,CP,PP,SP,RP,CF,PF,SF,UG,UP,UF,BF,UK",
        table_orthology=None,
        orthology_reference_schema=None,
        filename_filter=None,
        separator="|",
    )

    (options, args) = E.Start(parser, add_database_options=True)

    if options.genomes:
        options.genomes = options.genomes.split(",")
    if options.priority:
        options.priority = options.priority.split(",")
    if options.sort:
        options.sort = options.sort.split(",")

    if len(options.sort) != len(options.priority):
        raise "different number of classes in sort order and priority order"

    subset = {}
    if options.filename_filter:
        data = map(lambda x: x[:-1].split(options.separator)[:3],
                   filter(lambda x: x[0] != "#", open(options.filename_filter, "r").readlines()))
        for s, p, g in data:
            if s not in subset:
                subset[s] = {}
            subset[s][p] = 1

    dbhandle = pgdb.connect(options.psql_connection)

    data = []

    # get data for queries
    statement = """
    SELECT rep_token,
    CASE WHEN nexons=1 THEN 'SG' ELSE 'CG' END, 0
    FROM %s.queries
    ORDER BY rep_token
    """ % options.genomes[0]

    cc = dbhandle.cursor()
    cc.execute(statement)
    r = cc.fetchall()
    cc.close()

    if options.loglevel >= 1:
        print "# retrieved %i lines from queries table in %s" % (len(r), options.genomes[0])
        sys.stdout.flush()

    data.append(GetCounts(r, options.priority))

    # get data from genomes
    for genome in options.genomes:

        if options.table_orthology:

            # check which direction the data has been entered
            statement = "SELECT COUNT(*) FROM %s WHERE schema1 = '%s' AND schema2 = '%s'" % (options.table_orthology,
                                                                                             options.orthology_reference_schema,
                                                                                             genome)

            cc = dbhandle.cursor()
            cc.execute(statement)
            n = cc.fetchone()[0]
            cc.close()

            if n > 0:
                a = options.orthology_reference_schema
                b = genome
                c = 2
            else:
                a = genome
                b = options.orthology_reference_schema
                c = 1

            statement = """
            SELECT DISTINCT o.gene_id, o.class, o.prediction_id
            FROM %s.overview AS o,
            %s AS l
            WHERE
            l.schema1 = '%s' AND l.schema2 = '%s' AND
            l.prediction_id%i = o.prediction_id AND
            gene_id > 0 ORDER BY gene_id
            """ % (genome,
                   options.table_orthology,
                   a, b, c)

        else:
            statement = """
            SELECT gene_id, class, prediction_id FROM %s.overview WHERE gene_id != '0' ORDER BY gene_id, class
            """ % genome

        cc = dbhandle.cursor()
        cc.execute(statement)
        r = cc.fetchall()
        cc.close()

        if options.loglevel >= 1:
            print "# retrieved %i lines from %s" % (len(r), genome)
            sys.stdout.flush()

        # remove unwanted predictions
        if subset:
            r = filter(lambda x: str(x[2]) in subset[genome], r)
            if options.loglevel >= 1:
                print "# after filtering: %i lines from %s" % (len(r), genome)
                sys.stdout.flush()

        data.append(GetCounts(r, options.priority))

    print "class\tquery\t%s" % "\t".join(options.genomes)

    totals = [0] * (len(options.genomes) + 1)
    for q in options.sort:
        sys.stdout.write(q)
        for x in range(len(options.genomes) + 1):
            sys.stdout.write("\t%i" % data[x][q])
            totals[x] += data[x][q]
        sys.stdout.write("\n")

    sys.stdout.write("all\t%s\n" % "\t".join(map(lambda x: "%i" % x, totals)))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
