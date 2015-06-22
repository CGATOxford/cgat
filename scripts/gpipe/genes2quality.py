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
gpipe/genes2quality.py - 
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

   python gpipe/genes2quality.py --help

Type::

   python gpipe/genes2quality.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import re
import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import pgdb


""" program $Id: gpipe/genes2quality.py 2781 2009-09-10 11:33:14Z andreas $

return a table mapping genes to quality codes.

As there can be several transcripts per gene, this script selects the highest quality
transcript as representative of the gene quality.

If the method "filter" is selected, the script outputs transcripts that
should be removed.

The following filtering options are available:

--filter-nbest: keeps n best transcripts per gene
--filter-ds: keeps only transcripts with a ds of less than #
"""

parser = E.OptionParser(
    version="%prog version: $Id: gpipe/genes2quality.py 2781 2009-09-10 11:33:14Z andreas $")


def ProcessGene(outfile, gene_id, qualities, options):
    """output best quality transcript for gene.

    The best quality transcript is decided in two steps:

    1. According to the quality index, i.e., genes better than pseudogenesx
    2. According to a quality score. The higher, the better
    """
    noutput = 0
    if options.method == "geneinfo":

        qq = map(lambda x: x[0], qualities)

        best = None
        for x in options.priority:
            if x in qq:
                best = x
                break

        if best:
            # select best transcript (=highest score):
            s = []
            for quality, id, index in qualities:
                if quality == best:
                    s.append((index, id))

            s.sort()

            best_id = s[-1][1]

            if options.filter_quality and best not in options.filter_quality:
                return 0
            outfile.write("%s\t%i\t%s\t%s\n" %
                          (str(gene_id), len(qualities), best_id, best))
        else:
            outfile.write("%s\t%i\t%s\t%s\n" %
                          (str(gene_id), len(qualities), "", ""))
        noutput += 1

    elif options.method == "filter":

        # select nbest transcripts + all CGs and where gene_id == prediction_id
        # write a list of predictions to be removed from genes
        found = 0

        qualities.sort(lambda x, y: cmp(x[2], y[2]))
        qualities.reverse()

        for quality, id, index in qualities[options.filter_nbest:]:
            if quality not in options.filter_keep and \
                    gene_id != id:
                outfile.write("%s\t%s\t%s\t%5.2f\n" %
                              (str(id), str(gene_id), quality, index))
                noutput += 1
    else:
        raise "Unknown method %s" % options.method

    return noutput


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

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
                      help="filename with transcript information.")
    parser.add_option("-a", "--infile-fasta", dest="filename_fasta", type="string",
                      help="filename with transcript information. The sequence lengths will be used to select the longest transcript in case of ties.")
    parser.add_option("-m", "--method", dest="method", type="string",
                      help="method to apply [geneinfo|filter]")
    parser.add_option("-n", "--filter-nbest", dest="filter_nbest", type="int",
                      help="when filtering, take ## best predictions.")
    parser.add_option("-u", "--filter-quality", dest="filter_quality", type="string",
                      help="only print genes with quality codes.")

    parser.add_option("-s", "--filter-ds", dest="filter_ds", type="float",
                      help="filtering predictions whose ds is larger than the threashould")

    parser.set_defaults(
        priority="CG,PG,SG,RG,CP,PP,SP,RP,CF,PF,SF,UG,UP,UF,BF,UK",
        sort="CG,PG,SG,RG,CP,PP,SP,RP,CF,PF,SF,UG,UP,UF,BF,UK",
        table_genes=None,
        table_quality=None,
        table_predictions=None,
        filename_input=None,
        filename_fasta=None,
        separator="|",
        method="geneinfo",
        filter_keep="CG",
        filter_nbest=0,
        filter_ds=None,
        filter_quality=None
    )

    (options, args) = E.Start(parser, add_database_options=True)

    if options.priority:
        options.priority = options.priority.split(",")
    if options.filter_keep:
        options.filter_keep = options.filter_keep.split(",")
    if options.filter_quality:
        options.filter_quality = options.filter_quality.split(",")

    if options.filename_fasta:

        infile = open(options.filename_fasta, "r")
        sequences = Genomics.ReadPeptideSequences(infile)
        rows = []
        for key, sequence in sequences:
            s, t, g, q = key.split(options.separator)
            # use sequence length as score
            rows.append((g, t, q, len(sequence)))

        rows.sort()
        infile.close()

    elif options.filename_input:

        infile = open(options.filename_input, "r")

        rows = []
        for line in infile:
            if line[0] == "#":
                continue
            data = re.split("\s+", line[:-1])
            s, t, g, q = data[0].split(options.separator)
            rows.append((g, t, q, 1.0))

        rows.sort()
        infile.close()

    else:

        if not options.tablename_genes:
            raise "please specify a table with gene information."
        if not options.tablename_quality:
            raise "please specify a table with quality information."
        if not options.tablename_predictions:
            raise "please specify a table with prediction information."

        dbhandle = pgdb.connect(options.psql_connection)

        data = []

        # get data for queries: Use score or kaks (if table for kaks is given)
        if options.tablename_kaks:
            statement = """
            SELECT g.gene_id, g.prediction_id, q.class, -k.ds
            FROM %s AS g, %s as q, %s AS p, %s AS k
            WHERE g.prediction_id = q.prediction_id AND
            p.prediction_id = q.prediction_id AND
            k.prediction_id = p.prediction_id
            ORDER BY g.gene_id
            """ % (options.tablename_genes,
                   options.tablename_quality,
                   options.tablename_predictions,
                   options.tablename_kaks)
        else:
            statement = """
            SELECT g.gene_id, g.prediction_id, q.class, p.score
            FROM %s AS g, %s as q, %s AS p
            WHERE g.prediction_id = q.prediction_id AND
            p.prediction_id = q.prediction_id
            ORDER BY g.gene_id
            """ \
            % (options.tablename_genes,
               options.tablename_quality,
               options.tablename_predictions)

        cc = dbhandle.cursor()
        cc.execute(statement)
        rows = cc.fetchall()
        cc.close()

    if options.loglevel >= 1:
        print "# retrieved %i lines." % len(rows)
        sys.stdout.flush()

    if len(rows) == 0:
        E.Stop()
        sys.exit(0)

    last_gene, last_id, last_quality, last_score = rows[0]
    qualities = [(last_quality, last_id, last_score), ]

    if options.method == "geneinfo":
        options.stdout.write("gene_id\tntranscripts\tbest_id\tbest_quality\n")

    elif options.method == "filter":
        options.stdout.write("prediction_id\tgene_id\tquality\tindex\n")

    ninput, noutput, nfiltered_ds = 0, 0, 0

    do_filter_ds = options.tablename_kaks and options.method == "filter" and options.filter_ds is not None

    if do_filter_ds and options.loglevel >= 1:
        options.stdlog.write(
            "# filtering values by ds <= %f" % options.filter_ds)

    for this_gene, this_id, this_quality, this_score in rows[1:]:

        if this_gene == "0":
            continue

        ninput += 1

        # ds is negative ds in the select statement
        if do_filter_ds and -this_score >= options.filter_ds:
            nfiltered_ds += 1
            options.stdout.write("%s\t%s\t%s\t%5.2f\n" % (
                str(this_id), str(this_gene), this_quality, this_score))
            continue

        if this_gene != last_gene:

            n = ProcessGene(options.stdout, last_gene, qualities, options)
            noutput += n
            qualities = []

        qualities.append((this_quality, this_id, this_score))

        last_gene = this_gene

    n = ProcessGene(options.stdout, last_gene, qualities, options)
    noutput += n

    if options.loglevel >= 1:
        options.stdlog.write(
            "# ninput=%i, noutput=%i, nfiltered_ds=%i\n" % (ninput, noutput, nfiltered_ds))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
