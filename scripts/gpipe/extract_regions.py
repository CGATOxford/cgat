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
gpipe/extract_regions.py - 
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

   python gpipe/extract_regions.py --help

Type::

   python gpipe/extract_regions.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import getopt
import tempfile
import time
import optparse
import math

import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IOTools as IOTools
import pgdb
import CGAT.Regions as Regions

USAGE = """python %s [OPTIONS] 

Version: $Id: gpipe/extract_regions.py 2781 2009-09-10 11:33:14Z andreas $

get sequences/sequence regions for a list of predictions and output a fasta-formatted
file.

""" % sys.argv[0]


def getMRNAs(dbhandle, schema, options, prediction_ids):

    xfrom = ""
    xwhere = ""
    if prediction_ids == "all":
        where = ""
    elif prediction_ids == "nr":
        xfrom = "%s.redundant AS r," % schema
        where = "AND p.prediction_id = r.rep_prediction_id AND r.rep_prediction_id = r.mem_prediction_id"
    else:
        where = "AND p.prediction_id  IN ('%s')" % "','".join(
            map(str, prediction_ids))

    # select genomic location for all predictions
    # retrieve: prediction_id, sbjct_token, sbjct_strand, sbjct_genome_from,
    # sbjct_genome_to
    statement = """
    SELECT prediction_id, sbjct_token, sbjct_strand, sbjct_genome_from, sbjct_genome_to
    FROM
    %s
    %s.%s AS p
    WHERE True %s
    """ % (xfrom, schema, options.tablename_predictions,
           where)

    cc = dbhandle.cursor()
    cc.execute(statement)
    locations = cc.fetchall()
    cc.close()

    return locations


def getExons(dbhandle, schema, options, prediction_ids):

    xfrom = ""
    where = ""
    if prediction_ids == "all":
        where = ""
    elif prediction_ids == "nr":
        xfrom = "%s.redundant AS r," % schema
        where = "AND p.prediction_id = r.rep_prediction_id AND r.rep_prediction_id = r.mem_prediction_id"
    else:
        where = "AND p.prediction_id  IN ('%s')" % "','".join(
            map(str, prediction_ids))

    # select genomic location for all introns
    # filter on genome_exon_to > 0 in order to eliminate
    # all non-matching exons to query.
    # retrieve: prediction_id, sbjct_token, sbjct_strand, sbjct_genome_from,
    # sbjct_genome_to
    statement = """
    SELECT DISTINCT
    p.prediction_id, sbjct_token, sbjct_strand, genome_exon_from, genome_exon_to
    FROM
    %s.%s AS p,
    %s
    %s.%s AS e
    WHERE p.prediction_id = e.prediction_id 
    %s AND
    genome_exon_to > 0
    GROUP BY p.prediction_id, sbjct_token, sbjct_strand, genome_exon_from, genome_exon_to
    ORDER BY p.prediction_id, genome_exon_from
    """ % (schema, options.tablename_predictions,
           xfrom,
           schema, options.tablename_exons,
           where)

    cc = dbhandle.cursor()
    cc.execute(statement)
    locations = cc.fetchall()
    cc.close()

    return locations


def getIntrons(dbhandle, schema, options, prediction_ids):

    result = getExons(dbhandle, schema, options, prediction_ids)

    locations = []
    last_to = None
    last_id = None
    for id, sbjct_token, sbjct_strand, sbjct_from, sbjct_to in result:
        if last_id != id:
            last_to = None
            last_id = id
        if last_to:
            locations.append(
                (id, sbjct_token, sbjct_strand, last_to, sbjct_from))
        last_to = sbjct_to

    return locations

##########################################################################


def getOrthologs(dbhandle, options, prediction_ids):

    where = "AND l.prediction_id1 IN ('%s')" % "','".join(
        map(str, prediction_ids))

    # select genomic location for all predictions
    # retrieve: prediction_id, sbjct_token, sbjct_strand, sbjct_genome_from,
    # sbjct_genome_to
    statement = """
    SELECT schema2, prediction_id2
    FROM
    %s AS l
    WHERE l.schema1 = '%s' AND l.schema2 != l.schema1
    %s
    """ % (options.tablename_orthologs, options.schema, where )

    cc = dbhandle.cursor()
    cc.execute(statement)
    result = cc.fetchall()
    cc.close()

    orthologs = {}
    for schema, prediction_id in result:
        if schema not in orthologs:
            orthologs[schema] = []
        orthologs[schema].append(prediction_id)

    return orthologs

##########################################################################


def GetIdentifierInfo(dbhandle, schema, options, prediction_ids):
    """get additional identifier info."""

    xfrom = ""
    where = ""

    if prediction_ids == "all":
        where = ""
    elif prediction_ids == "nr":
        xfrom = "%s.redundant AS r," % schema
        where = "AND p.prediction_id = r.rep_prediction_id AND r.rep_prediction_id = r.mem_prediction_id"
    else:
        where = "AND p.prediction_id  IN ('%s')" % "','".join(
            map(str, prediction_ids))

    statement = """
    SELECT p.prediction_id, g.gene_id, q.class
    FROM
    %s.%s AS p,
    %s.%s AS g,
    %s
    %s.%s AS q
    WHERE p.prediction_id = g.prediction_id AND 
    p.prediction_id = q.prediction_id
    %s
    GROUP BY p.prediction_id, g.gene_id, q.class
    """ % (schema, options.tablename_predictions,
           schema, options.tablename_genes,
           xfrom,
           schema, options.tablename_quality,
           where)

    cc = dbhandle.cursor()
    cc.execute(statement)
    result = cc.fetchall()
    cc.close()

    info = {}
    for prediction_id, gene_id, quality in result:
        info[prediction_id] = (gene_id, quality)

    return info


##########################################################################
def ProcessPredictions(dbhandle, schema, options, prediction_ids, taboo_regions=None):
    """print required regions for a set of prediction_ids from a given schema."""
    # Step 2 :

    if "%s" in options.genome_file:
        fasta = IndexedFasta.IndexedFasta(options.genome_file % schema)
    else:
        fasta = IndexedFasta.IndexedFasta(options.genome_file)

    if options.type == "mrnas":
        locations = getMRNAs(dbhandle, schema, options, prediction_ids)
    elif options.type == "introns":
        locations = getIntrons(dbhandle, schema, options, prediction_ids)
    elif options.type == "exons":
        locations = getExons(dbhandle, schema, options, prediction_ids)
    elif options.type == "cds":
        locations = getExons(dbhandle, schema, options, prediction_ids)
        options.join_regions = True

    if options.loglevel >= 1:
        options.stdlog.write(
            "# %s: retrieved %i locations.\n" % (schema, len(locations)))

    if options.loglevel >= 3:
        for location in locations:
            options.stdlog.write("# %s\n" % (str(location)))

    if len(locations) == 0:
        return

    if options.id_format == "full":
        # retrieve quality and gene information
        extra_info = GetIdentifierInfo(
            dbhandle, schema, options, prediction_ids)

    contigs = set(map(lambda x: x[1], locations))

    if options.loglevel >= 1:
        options.stdlog.write("# %s: %i locations are on %i contigs.\n" %
                             (schema, len(locations), len(contigs)))
        options.stdlog.flush()

    # retrieve contig sizes
    statement = """
    SELECT sbjct_token, size FROM %s.contigs""" % ( schema )
    cc = dbhandle.cursor()
    cc.execute(statement)
    result = cc.fetchall()
    cc.close()

    contig_size = {}
    for sbjct_token, size in result:
        contig_size[sbjct_token] = size

    # sort locations by contig, that way sequence retrieval is faster
    # sort also by identifier
    locations.sort(lambda x, y: cmp((x[2], x[0]), (y[2], y[0])))
    tokens = {}
    noutput = 0
    nskipped_length = 0
    nskipped_overlap = 0
    nskipped_doubles = 0

    output_regions = set()

    last_token = None
    this_is_last = False

    if options.join_regions:
        # add dummy to locations
        locations.append(("", "", "+", 0, 0))

    for token, sbjct_token, sbjct_strand, sbjct_genome_from, sbjct_genome_to in locations:

        if token != "":

            if token not in tokens:
                tokens[token] = 0
            tokens[token] += 1

            sbjct_genome_from -= (options.extend_region -
                                  options.shorten_region)
            sbjct_genome_from = max(0, sbjct_genome_from)
            sbjct_genome_to += (options.extend_region - options.shorten_region)
            sbjct_genome_to = min(sbjct_genome_to, contig_size[sbjct_token])

            if sbjct_genome_to - sbjct_genome_from <= 0:
                continue

            lgenome = fasta.getLength(sbjct_token)

            forward_sbjct_genome_from, forward_sbjct_genome_to = Genomics.ToForwardCoordinates(sbjct_genome_from,
                                                                                               sbjct_genome_to,
                                                                                               sbjct_strand,
                                                                                               lgenome)
            if taboo_regions:

                if options.taboo_regions == "both":
                    overlaps = taboo_regions.getOverlaps(sbjct_token, sbjct_strand,
                                                         forward_sbjct_genome_from, forward_sbjct_genome_to)
                else:
                    overlaps = taboo_regions.getOverlaps(sbjct_token, sbjct_strand,
                                                         sbjct_genome_from, sbjct_genome_to)

                if overlaps:
                    if options.loglevel >= 2:
                        options.stdlog.write(
                            "# %s eliminated due to overlap with taboo region.\n" % (token))
                    nskipped_overlap += 1
                    continue

            genomic_sequence = fasta.getSequence(sbjct_token,
                                                 sbjct_strand,
                                                 sbjct_genome_from,
                                                 sbjct_genome_to)

            if len(genomic_sequence) == 0:
                options.stderr.write("## warning: sequence for prediction %s is emtpy: %s:%s:%i:%i\n" %
                                     (str(token), sbjct_token, sbjct_strand,
                                      sbjct_genome_from, sbjct_genome_to))
                options.stderr.flush()

            # adjust sbjct_genome_to if it went outside of contig boundaries.
            sbjct_genome_to = len(genomic_sequence) + sbjct_genome_from

        do_output = False
        # decide whether segment is to be joined
        if options.join_regions:

            if last_token == token:
                do_output = False
                fragments.append(str(genomic_sequence))
                max_sbjct_genome_to = sbjct_genome_to
            else:
                if last_token is not None:
                    output_sbjct_genome_from = min_sbjct_genome_from
                    output_sbjct_genome_to = max_sbjct_genome_to
                    output_genomic_sequence = "".join(fragments)
                    output_token = last_token
                    do_output = True

                min_sbjct_genome_from = sbjct_genome_from
                max_sbjct_genome_to = sbjct_genome_to
                fragments = [str(genomic_sequence)]
                last_token = token
        else:
            output_sbjct_genome_from = sbjct_genome_from
            output_sbjct_genome_to = sbjct_genome_to
            output_genomic_sequence = str(genomic_sequence)
            output_token = token
            do_output = True

        if not do_output:
            continue

        # build output token
        if options.id_format == "id":
            output_token = str(output_token)
        elif options.id_format == "schema-id":
            output_token = "%s%s%s" % (
                schema, options.separator, str(output_token))
        elif options.id_format == "full":
            gene_id, quality = extra_info[output_token]
            output_token = options.separator.join(
                map(str, (schema, output_token, gene_id, quality)))

        if options.forward_coordinates:
            sbjct_genome_from, sbjct_genome_to = forward_sbjct_genome_from, forward_sbjct_genome_to

        # Filter output by length
        if len(output_genomic_sequence) < options.min_length:
            nskipped_length += 1
            return

        # Filter output for redundant entries
        k = "%s-%s-%i-%i" % (sbjct_token, sbjct_strand,
                             output_sbjct_genome_from, output_sbjct_genome_to)
        if k in output_regions:
            nskipped_doubles += 1
            return

        output_regions.add(k)

        if options.output_coordinate_format == "full":
            coordinates = "%s:%s:%i:%i" % (sbjct_token,
                                           sbjct_strand,
                                           output_sbjct_genome_from,
                                           output_sbjct_genome_to)

        elif options.output_coordinate_format == "long":
            coordinates = "%s:%s:%i:%i:%i" % (sbjct_token,
                                              sbjct_strand,
                                              sbjct_genome_from,
                                              sbjct_genome_to,
                                              lgenome)

        if options.output_format == "fasta":
            # print fasta formatted output
            if options.fasta_format == "id-coordinates":
                options.stdout.write(">%s %s\n%s\n" % (output_token, coordinates,
                                                       str(output_genomic_sequence)))
            elif options.fasta_format == "coordinates":
                options.stdout.write(">%s:%s:%i:%i\n%s\n" % (coordinates,
                                                             str(output_genomic_sequence)))
            elif options.fasta_format == "schema-coordinates":
                options.stdout.write(">%s|%s:%s:%i:%i\n%s\n" % (schema,
                                                                coordinates,
                                                                str(output_genomic_sequence)))

        elif options.output_format == "table":
            # print identifier and sequence in tab separated table
            options.stdout.write("%s\t%s\t%s\t%i\t%i\ts\n" % (output_token,
                                                              sbjct_token, sbjct_strand,
                                                              output_sbjct_genome_from, output_sbjct_genome_to,
                                                              str(output_genomic_sequence)))

        elif options.output_format == "region":
            options.stdout.write("%s\t%s\t%s\t%i\t%i\n" % (output_token,
                                                           sbjct_token, sbjct_strand,
                                                           output_sbjct_genome_from, output_sbjct_genome_to))

        noutput += 1

        options.stdout.flush()

    if options.loglevel >= 1:
        if not prediction_ids:
            prediction_ids = []
        options.stdlog.write("# %s: ninput=%i, npredictions=%i, nlocations=%i, noutput=%i, nskipped_length=%i, nskipped_overlap=%i, nskipped_doubles=%i\n" % (
            schema, len(prediction_ids), len(tokens), len(locations), noutput, nskipped_length, nskipped_overlap, nskipped_doubles))


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gpipe/extract_regions.py 2781 2009-09-10 11:33:14Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="pattern to look for sequence filename.")

    parser.add_option("-i", "--ids", dest="ids", type="string",
                      help="comma separated list of prediction ids. Use 'all' to use all predictions.")

    parser.add_option("-f", "--filename-ids", dest="filename_ids", type="string",
                      help="filename with prediction ids.")

    parser.add_option("-t", "--sequence-type", dest="type", type="choice",
                      choices=("mrnas", "introns", "exons", "cds"),
                      help="type to output.")

    parser.add_option("-e", "--extend-region", dest="extend_region", type="int",
                      help="regions are extended by this margin at either end.")

    parser.add_option("-r", "--shorten-region", dest="shorten_region", type="int",
                      help="regions are shortened by this margin at either end.")

    parser.add_option("-m", "--min-interval-length", dest="min_length", type="int",
                      help="minimum length of segment.")

    parser.add_option("-s", "--schema", dest="schema", type="string",
                      help="schema to take data from.")

    parser.add_option("-o", "--output-format", dest="output_format", type="choice",
                      choices=("fasta", "table", "region"),
                      help="output formats.")

    parser.add_option("--fasta-format", dest="fasta_format", type="choice",
                      choices=(
                          "id-coordinates", "coordinates", "schema-coordinates"),
                      help="output formats for fasta formatted headers.")

    parser.add_option("--orthologs", dest="orthologs", action="store_true",
                      help="lookup up orthologs of prediction ids.")

    parser.add_option("--multiple", dest="multiple", action="store_true",
                      help="""lookup up predictions in multiple species.
                       Identifiers should be given as schema|prediction_id[|additional_fields].
                       Note that the genome file locations have to be consistent.""")

    parser.add_option("--id-format", dest="id_format", type="choice",
                      choices=("id", "schema-id", "full"),
                      help="output format for ids.")

    parser.add_option("--taboo-regions", dest="taboo_regions", type="choice",
                      choices=("same", "both"),
                      help="check for overlap in same/both strands.")

    parser.add_option("--filename-taboo-regions", dest="filename_taboo_regions", type="string",
                      help="filename with information about taboo regions.")

    parser.add_option("--is-forward-coordinates", dest="forward_coordinates", action="store_true",
                      help="output coordinates are forward coordinates.")

    parser.add_option("--join-regions", dest="join_regions", action="store_true",
                      help="join regions with the same identifier.")

    parser.add_option("--output-coordinate-format", dest="output_coordinate_format", type="choice",
                      choices=("full", "long"),
                      help="""output format of coordinates. Output format is contig:strand:from:to in zero based
/forward/reverse strand coordinates in open/closed notation. 'long' includes the contig length as fifth field"""  )

    parser.set_defaults(
        genome_file="genome",
        identifiers=None,
        filename_ids="-",
        ids=None,
        extend_region=0,
        shorten_region=0,
        join_regions=False,
        tablename_predictions="predictions",
        tablename_exons="exons",
        tablename_genes="genes",
        tablename_quality="quality",
        tablename_orthologs="orthology_pairwise1v5.orthologlinks_first",
        schema=None,
        output_format="fasta",
        fasta_format="id-coordinates",
        type="mrnas",
        min_length=1,
        id_format="id",
        mmultiple=False,
        separator="|",
        filename_taboo_regions=False,
        forward_coordinates=False,
        output_coordinate_format="full",
    )

    (options, args) = E.Start(parser, add_database_options=True)

    if options.orthologs:
        options.id_format = "schema-id"

    # database handle for connecting to postgres
    dbhandle = pgdb.connect(options.psql_connection)

    # Step 1 : Input of predictions

    # read identifiers from file, command line arguments or stdin.

    if options.ids in ("all", "nr"):
        prediction_ids = options.ids
        if options.loglevel >= 1:
            options.stdlog.write("# using all prediction ids.\n")
            options.stdlog.flush()
    elif options.ids:
        prediction_ids = options.ids.split(",")
    elif len(args) > 0:
        prediction_ids = args

    elif options.filename_ids:
        prediction_ids = []

        if options.filename_ids == "-":
            prediction_ids += IOTools.ReadList(sys.stdin)[0]
        elif options.filename_ids:
            prediction_ids += IOTools.ReadList(
                open(options.filename_ids, "r"))[0]

        if len(prediction_ids) == 0:
            raise "no prediction identifiers given."

        if options.loglevel >= 1:
            options.stdlog.write(
                "# read %i prediction ids.\n" % len(prediction_ids))
            options.stdlog.flush()

    if options.filename_taboo_regions:
        # Note: the input has to be in forward coordinates in order for option
        # "both" to work.
        taboo_regions = Regions.RegionFilter()
        if options.taboo_regions == "both":
            ignore_strand = True
        else:
            ignore_strand = False
        taboo_regions.readFromFile(
            open(options.filename_taboo_regions, "r"), ignore_strand=ignore_strand)
    else:
        taboo_regions = None

    if options.orthologs:
        orthologs = getOrthologs(dbhandle, options, prediction_ids)
        for schema, ids in orthologs.items():
            ProcessPredictions(dbhandle, schema, options, ids, taboo_regions)
    elif options.multiple:
        # convert identifiers
        prediction_ids = map(
            lambda x: x.split(options.separator), prediction_ids)
        prediction_ids.sort()
        last_schema = None
        ids = []
        for x in prediction_ids:
            if x[0] != last_schema:
                ProcessPredictions(dbhandle, last_schema, options, ids)
                ids = []
                last_schema = x[0]
            ids.append(x[1])
        ProcessPredictions(dbhandle, last_schema, options, ids)
    else:
        ProcessPredictions(
            dbhandle, options.schema, options, prediction_ids, taboo_regions)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
