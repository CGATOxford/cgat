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
gpipe/select_transcripts_per_gene.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
Remove transcripts from a gene that are likely predictions by paralogs.

Usage
-----

Example::

   python gpipe/select_transcripts_per_gene.py --help

Type::

   python gpipe/select_transcripts_per_gene.py --help

for command line help.


Documentation
-------------

Code
----

'''
import sys
import CGAT.Experiment as E
import pgdb


def Write(outfile, gene_id, qualities, priority):
    """write best quality for gene.
    """

    qq = map(lambda x: x[0], qualities)

    best = None
    for x in priority:
        if x in qq:
            best = x
            break

    if best:
        # select best transcript (=highest score):
        s = []
        for quality, id, index in qualities:
            if x == best:
                s.append((index, id))

        s.sort()
        best_id = s[-1][1]

        outfile.write("%s\t%i\t%s\t%s\n" %
                      (str(gene_id), len(qualities), best_id, best))
    else:
        outfile.write("%s\t%i\t%s\t%s\n" %
                      (str(gene_id), len(qualities), "", ""))


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gpipe/select_transcripts_per_gene.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-i", "--priority", dest="priority", type="string",
                      help="quality priority.")
    parser.add_option("-q", "--table-quality", dest="tablename_quality", type="string",
                      help="table with quality information.")
    parser.add_option("-g", "--table-genes", dest="tablename_genes", type="string",
                      help="table with gene information.")
    parser.add_option("-p", "--table-predictions", dest="tablename_predictions", type="string",
                      help="table with prediction information.")
    parser.add_option("-k", "--table-kaks", dest="tablename_kaks", type="string",
                      help="table with kaks information.")
    parser.add_option("-f", "--infile-transcripts", dest="filename_input", type="string",
                      help="table with transcript information.")

    parser.set_defaults(
        priority="CG,PG,SG,RG,CP,PP,SP,RP,CF,PF,SF,UG,UP,UF,BF,UK",
        sort="CG,PG,SG,RG,CP,PP,SP,RP,CF,PF,SF,UG,UP,UF,BF,UK",
        table_genes=None,
        table_quality=None,
        table_predictions=None,
        table_kaks=None,
        filename_input=None,
        separator="|",
    )

    (options, args) = E.Start(parser, add_database_options=True)

    if options.priority:
        options.priority = options.priority.split(",")

    if options.filename_input:
        infile = open(options.filename_input, "r")

        rows = []
        for line in infile:
            if line[0] == "#":
                continue
            data = line[:-1].split("\t")
            s, t, g, q = data[0].split(options.separator)
            rows.append((g, t, q, 1.0))

        infile.close()

    else:
        if not options.tablename_genes:
            raise "please specify a table with gene information."
        if not options.tablename_quality:
            raise "please specify a table with quality information."
        if not options.tablename_predictions:
            raise "please specify a table with prediction information."
        if not options.tablename_kaks:
            raise "please specify a table with kaks information."

        dbhandle = pgdb.connect(options.psql_connection)

        data = []

        # get data for queries
        statement = """
        SELECT g.gene_id, g.prediction_id, q.class, k.ds
        FROM %s AS g, %s as q, %s AS p, %s AS k
        WHERE g.prediction_id = q.prediction_id AND
        p.prediction_id = q.prediction_id AND
        k.prediction_id = p .predction_id
        g.gene_id > 0
        ORDER BY g.gene_id
        """ % (options.tablename_genes, options.tablename_quality, options.tablename_predictions)

        cc = dbhandle.cursor()
        cc.execute(statement)
        rows = cc.fetchall()
        cc.close()

    if options.loglevel >= 1:
        print "# retrieved %i lines." % len(rows)
        sys.stdout.flush()

    last_gene, last_id, last_quality, last_score = rows[0]
    qualities = [(last_quality, last_id, last_score), ]

    for this_gene, this_id, this_quality, this_score in rows[1:]:

        if this_gene != last_gene:
            Write(sys.stdout, last_gene, qualities, options.priority)
            qualities = []

        qualities.append((this_quality, this_id, this_score))

        last_gene = this_gene

    Write(sys.stdout, last_gene, qualities, options.priority)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
