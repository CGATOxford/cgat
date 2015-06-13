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
gpipe/export_aaa.py - 
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

   python gpipe/export_aaa.py --help

Type::

   python gpipe/export_aaa.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import time
import CGAT.Experiment as E
import pgdb


USAGE = """python %s [OPTIONS] < in > out

output gene predictions for AAA submission.

Version: $Id: gpipe/export_aaa.py 2781 2009-09-10 11:33:14Z andreas $

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
""" % sys.argv[0]

HEADER = """##gff-version   3
#species: %(species)s
#assembly-id: %(assembly)s
#annotation-group-id: %(group_id)s_%(method_id)s
#algorithm: %(group_id)s_%(method_id)s = Oxford gene prediction pipeline based on the program exonerate.
#templates: Gene models from FlyBase 4.2.1 obtained from ENSEMBL version 37
#authors: Andreas Heger and Chris Ponting at firstname.lastname@anat.ox.ac.uk
#date: %(date)s
"""

# ------------------------------------------------------------------------


def WriteGeneTrack(outfile, dbhandle,
                   current_id,
                   filter_options,
                   options, parameters):

    statement = """
    SELECT p.sbjct_token, p.sbjct_strand, g.gene_id,
    MIN(CASE WHEN p.sbjct_strand = '+' THEN e.genome_exon_from+c.start 
    WHEN p.sbjct_strand = '-' THEN c.size-e.genome_exon_to+c.start END)+1 AS start, 
    MAX(CASE WHEN p.sbjct_strand = '+' THEN e.genome_exon_to+c.start 
    WHEN p.sbjct_strand = '-' THEN c.size-e.genome_exon_from+c.start END) AS end
    FROM
    %(tablename_predictions)s AS p,
    %(tablename_genes)s AS g,
    %(tablename_exons)s AS e,
    %(tablename_quality)s AS q,
    %(tablename_redundant)s AS r,
    %(tablename_contigs)s AS c    
    %(extra_tables)s
    WHERE g.prediction_id = p.prediction_id AND
    e.prediction_id = p.prediction_id AND
    q.prediction_id = p.prediction_id AND
    r.mem_prediction_id = p.prediction_id AND
    c.sbjct_token = p.sbjct_token AND
    e.genome_exon_from > 0 AND e.genome_exon_to > 0 AND
    g.gene_id > 0
    %(extra_where)s
    """ % (parameters)

    statement += " " + filter_options

    statement += " GROUP BY p.sbjct_token, p.sbjct_strand, g.gene_id"
    statement += " ORDER BY p.sbjct_token, start"

    if options.loglevel >= 2:
        options.stdlog.write(
            "statement for selecting genes:\n%s\n%s\n%s\n" % ("#" * 20, statement, "#" * 20))

    cc = dbhandle.cursor()
    cc.execute(statement)
    rr = cc.fetchall()
    cc.close()

    map_gene_id = {}

    for sbjct_token, sbjct_strand, gene_id, feature_from, feature_to in rr:

        if gene_id not in map_gene_id:
            map_gene_id[gene_id] = "%s_%s_%s_%s" % (
                options.species_id, options.group_id, options.method_id, current_id)
            current_id += 1

        outfile.write("\t".join((sbjct_token,
                                 options.group_id + "_" + options.method_id,
                                 "gene",
                                 str(feature_from),
                                 str(feature_to),
                                 ".",
                                 sbjct_strand,
                                 ".",
                                 "ID=%s" % map_gene_id[gene_id])) + "\n")

    return current_id, map_gene_id

# ------------------------------------------------------------------------


def WriteMRNATrack(outfile, dbhandle,
                   current_id, map_gene_id,
                   filter_options,
                   options, parameters):

    statement = """
    SELECT p.sbjct_token, p.sbjct_strand, g.gene_id, p.prediction_id,
    MIN(CASE WHEN p.sbjct_strand = '+' THEN e.genome_exon_from+c.start 
    WHEN p.sbjct_strand = '-' THEN c.size-e.genome_exon_to+c.start END)+1 AS start, 
    MAX(CASE WHEN p.sbjct_strand = '+' THEN e.genome_exon_to+c.start 
    WHEN p.sbjct_strand = '-' THEN c.size-e.genome_exon_from+c.start END) AS end
    FROM
    %(tablename_predictions)s AS p,
    %(tablename_genes)s AS g,
    %(tablename_exons)s AS e,
    %(tablename_quality)s AS q,
    %(tablename_redundant)s AS r,
    %(tablename_contigs)s AS c
    %(extra_tables)s
    WHERE g.prediction_id = p.prediction_id AND
    e.prediction_id = p.prediction_id AND
    q.prediction_id = p.prediction_id AND
    r.mem_prediction_id = p.prediction_id AND
    c.sbjct_token = p.sbjct_token AND
    e.genome_exon_from > 0 AND e.genome_exon_to > 0 AND
    g.gene_id > 0
    %(extra_where)s
    """ % (parameters)

    statement += " " + filter_options

    statement += " GROUP BY p.sbjct_token, p.sbjct_strand, g.gene_id, p.prediction_id"
    statement += " ORDER BY p.sbjct_token, start"

    if options.loglevel >= 2:
        options.stdlog.write(
            "statement for selecting mRNAs:\n%s\n%s\n%s\n" % ("#" * 20, statement, "#" * 20))

    cc = dbhandle.cursor()
    cc.execute(statement)
    rr = cc.fetchall()
    cc.close()

    map_prediction_id = {}
    for sbjct_token, sbjct_strand, gene_id, prediction_id, feature_from, feature_to in rr:

        if gene_id not in map_gene_id:
            options.stderr.write(
                "# warning: gene_id %i not in set\n" % gene_id)
            continue

        if prediction_id not in map_prediction_id:
            map_prediction_id[prediction_id] = "%s_%s_%s_%s" % (
                options.species_id, options.group_id, options.method_id, current_id)
            current_id += 1

        outfile.write("\t".join((sbjct_token,
                                 options.group_id + "_" + options.method_id,
                                 "mRNA",
                                 str(feature_from),
                                 str(feature_to),
                                 ".",
                                 sbjct_strand,
                                 ".",
                                 ";".join(("ID=%s" % map_prediction_id[prediction_id],
                                           "PARENT=%s" % map_gene_id[gene_id])),
                                 )) + "\n")

    return current_id, map_prediction_id

# ------------------------------------------------------------------------


def WriteCDSTrack(outfile, dbhandle,
                  current_id,
                  map_gene_id, map_prediction_id,
                  filter_options,
                  options, parameters):

    statement = """
    SELECT DISTINCT p.sbjct_token, p.sbjct_strand, g.gene_id, p.prediction_id,
    (CASE WHEN p.sbjct_strand = '+' THEN e.genome_exon_from+c.start 
    WHEN p.sbjct_strand = '-' THEN c.size-e.genome_exon_to+c.start END)+1 AS start, 
    CASE WHEN p.sbjct_strand = '+' THEN e.genome_exon_to+c.start 
    WHEN p.sbjct_strand = '-' THEN c.size-e.genome_exon_from+c.start END AS end
    FROM
    %(tablename_predictions)s AS p,
    %(tablename_genes)s AS g,
    %(tablename_exons)s AS e,
    %(tablename_quality)s AS q,
    %(tablename_redundant)s AS r,
    %(tablename_contigs)s AS c
    %(extra_tables)s
    WHERE g.prediction_id = p.prediction_id AND
    e.prediction_id = p.prediction_id AND
    q.prediction_id = p.prediction_id AND
    r.mem_prediction_id = p.prediction_id AND
    c.sbjct_token = p.sbjct_token AND
    e.genome_exon_from > 0 AND e.genome_exon_to > 0 AND
    g.gene_id > 0
    %(extra_where)s
    """ % (parameters)

    statement += " " + filter_options

    statement += " ORDER BY p.sbjct_token, start"

    if options.loglevel >= 2:
        options.stdlog.write(
            "statement for selecting cds:\n%s\n%s\n%s\n" % ("#" * 20, statement, "#" * 20))

    cc = dbhandle.cursor()
    cc.execute(statement)
    rr = cc.fetchall()
    cc.close()

    for sbjct_token, sbjct_strand, gene_id, prediction_id, feature_from, feature_to in rr:

        if gene_id not in map_gene_id:
            options.stderr.write(
                "# warning: gene_id %i not in set\n" % gene_id)
            continue

        if prediction_id not in map_prediction_id:
            options.stderr.write(
                "# warning: prediction_id %i not in set\n" % prediction_id)
            continue

        outfile.write("\t".join((sbjct_token,
                                 options.group_id + "_" + options.method_id,
                                 "CDS",
                                 str(feature_from),
                                 str(feature_to),
                                 ".",
                                 sbjct_strand,
                                 ".",
                                 ";".join(
                                     ("PARENT=%s" % map_gene_id[gene_id],)),
                                 )) + "\n")

    return current_id

# ------------------------------------------------------------------------


def WriteSequenceRegions(outfile, dbhandle,
                         options, parameters):

    statement = """
    SELECT c.sbjct_token, 1, c.size
    FROM
    %(tablename_contigs)s AS c    
    """ % (parameters)

    cc = dbhandle.cursor()
    cc.execute(statement)
    rr = cc.fetchall()
    cc.close()

    for sbjct_token, feature_from, feature_to in rr:
        outfile.write(" ".join(("##sequence-region",
                                sbjct_token,
                                str(feature_from),
                                str(feature_to),
                                )) + "\n")

# ------------------------------------------------------------------------


def WriteGFF(outfile, dbhandle, options):

    # pre-processing: set table names
    parameters = {}
    parameters['tablename_predictions'] = options.schema + \
        "." + options.tablename_predictions
    parameters['tablename_genes'] = options.schema + \
        "." + options.tablename_genes
    parameters['tablename_redundant'] = options.schema + \
        "." + options.tablename_redundant
    parameters['tablename_contigs'] = options.schema + \
        "." + options.tablename_contigs
    parameters['tablename_exons'] = options.schema + \
        "." + options.tablename_exons
    parameters['tablename_quality'] = options.schema + \
        "." + options.tablename_quality

    if options.tablename_orthologs:
        parameters['extra_tables'] = ",%s AS oo" % options.tablename_orthologs
        parameters['extra_where'] = "AND oo.prediction_id1 = p.prediction_id AND oo.schema1='%s' " % (
            options.schema)

    # pre-processing: set filtering options
    if options.set == "full":
        filter_options = ""
    elif options.set == "filtered":
        filter_options = "AND r.rep_prediction_id = r.mem_prediction_id"
    elif options.set == "clean":
        filter_options = """AND r.rep_prediction_id = r.mem_prediction_id AND 
        q.class IN ('%s')
        """ % "','".join(options.quality_clean_set)

    # start numbering from the first id
    current_id = options.first_id

    # output
    outfile.write(HEADER % {'species': options.species,
                            'assembly': options.assembly_id,
                            'group_id': options.group_id,
                            'method_id': options.method_id,
                            'date': time.strftime("%Y%m%d")})

    WriteSequenceRegions(outfile, dbhandle, options, parameters)

    current_id, map_gene_id = WriteGeneTrack(outfile, dbhandle, current_id,
                                             filter_options, options, parameters)

    current_id, map_prediction_id = WriteMRNATrack(outfile, dbhandle, current_id,
                                                   map_gene_id,
                                                   filter_options, options, parameters)

    WriteCDSTrack(outfile, dbhandle, current_id,
                  map_gene_id, map_prediction_id,
                  filter_options, options, parameters)

    options.stdlog.write("# genome=%s, genes=%i, transcripts=%i\n" % (
        options.schema, len(map_gene_id), len(map_prediction_id)))
    options.stdlog.flush()

    return map_gene_id, map_prediction_id


def WriteMap(options, category, map_id):
    """write map between genes/predictions to new ids."""
    outfile = open(options.filename_pattern_map %
                   (options.assembly_id, category), "w")
    ids = map_id.keys()
    ids.sort()
    for id in ids:
        outfile.write("%i\t%s\n" % (id, map_id[id]))
    outfile.close()


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gpipe/export_aaa.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-s", "--set", dest="set",
                      help="set to use.", type="choice", choices=("full", "filtered", "clean"))
    parser.add_option("-b", "--filename-batch", dest="filename_batch", type="string",
                      help="filename with batch information.")
    parser.add_option("-o", "--filename-pattern-output", dest="filename_pattern_output", type="string",
                      help="filename with one %s for batch output.")
    parser.add_option("-m", "--filename-pattern-map", dest="filename_pattern_map", type="string",
                      help="filename with one %s for maps.")
    parser.add_option("-f", "--first-id", dest="first_id", type="int",
                      help="first number to be used for identifiers.")
    parser.add_option("--tablename-orthologs", dest="tablename_orthologs", type="string",
                      help="tablename with orthology information. Only orthologs are output.")

    parser.set_defaults(
        set="clean",
        filename_batch=None,
        filename_pattern_output="%s.gff",
        filename_pattern_map=None,
        species="Drosophila pseudoobscura",
        assembly_id="dpse_caf1",
        abbreviation="dpse",
        species_id="GA",
        method_id="GPI",
        group_id="OXFD",
        schema="dpse_vs_dmel8",
        tablename_predictions="predictions",
        tablename_genes="genes",
        tablename_redundant="redundant",
        tablename_exons="exons",
        tablename_contigs="contigs",
        tablename_quality="quality",
        tablename_orthologs=None,
        first_id=2000000,
        quality_clean_set='CG,PG,SG',
    )

    (options, args) = E.Start(parser,
                              add_database_options=True,
                              add_pipe_options=True)

    options.quality_clean_set = options.quality_clean_set.split(",")

    dbhandle = pgdb.connect(options.psql_connection)

    if options.filename_batch:
        infile = open(options.filename_batch)
        for line in infile:
            if line[0] == "#":
                continue
            options.species, options.abbreviation, options.species_id, options.assembly_id, options.schema = line[
                :-1].split("\t")[:5]
            filename = options.filename_pattern_output % options.assembly_id
            outfile = open(filename, "w")

            if options.loglevel >= 1:
                options.stdlog.write(
                    "# writing data from %s to %s\n" % (options.schema, filename))
                options.stdlog.flush()

            map_gene_id, map_prediction_id = WriteGFF(
                outfile, dbhandle, options)

            outfile.close()

            if options.filename_pattern_map:
                WriteMap(options, "genes", map_gene_id)
                WriteMap(options, "predictions", map_prediction_id)

        infile.close()

    else:
        outfile = options.stdout
        WriteGFF(outfile, dbhandle, options)
        if options.filename_pattern_map:
            WriteMap(options, "genes", map_gene_id)
            WriteMap(options, "predictions", map_prediction_id)

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
