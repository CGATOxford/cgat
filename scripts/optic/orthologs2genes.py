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
optic/orthologs2genes.py - 
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

   python optic/orthologs2genes.py --help

Type::

   python optic/orthologs2genes.py --help

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
import optparse
import time
import warnings

import CGAT.Experiment as E
import CGAT.Orthologs as Orthologs
import pgdb
import CGAT.Components as Components

# ignore warnings from networkx/matplotlib that a display
# can not be found
with warnings.catch_warnings():
    warnings.simplefilter("ignore")

"""read orthology assignments and re-assign gene
lists.

Overlapping orthologous transcripts are assigned to the same gene.
"""

# ------------------------------------------------------------------------


def GetPrediction2Location(dbhandle, schema,
                           tablename_predictions="predictions"):

    map_transcript2location = {}

    statement = """
    SELECT prediction_id, sbjct_token, sbjct_strand,
    sbjct_genome_from, sbjct_genome_to
    FROM %s.%s 
    """ % \
        (schema, tablename_predictions)

    cc = dbhandle.cursor()
    cc.execute(statement)
    rr = cc.fetchall()
    cc.close()

    for prediction_id, sbjct_token, sbjct_strand, sbjct_genome_from, sbjct_genome_to in rr:
        map_transcript2location[str(prediction_id)] = (
            sbjct_token, sbjct_strand, sbjct_genome_from, sbjct_genome_to)

    return map_transcript2location

# ------------------------------------------------------------------------


def GetPrediction2LocationFromFile(infile, options):
    """read exon file to get map of transcript to location."""

    map_transcript2location = {}

    for line in infile:
        if line[0] == "#":
            continue
        id, sbjct_token, sbjct_strand, phase, n, p1, p2, sbjct_genome_from, sbjct_genome_to = line[
            :-1].split("\t")

        sbjct_genome_from, sbjct_genome_to = int(
            sbjct_genome_from), int(sbjct_genome_to)
        schema, prediction_id = id.split(options.separator)[:2]

        if prediction_id not in map_transcript2location:
            map_transcript2location[prediction_id] = (
                sbjct_token, sbjct_strand, sbjct_genome_from, sbjct_genome_to)
        else:
            o_sbjct_token, o_sbjct_strand, o_sbjct_genome_from, o_sbjct_genome_to = map_transcript2location[
                prediction_id]
            map_transcript2location[prediction_id] = (sbjct_token, sbjct_strand,
                                                      min(sbjct_genome_from,
                                                          o_sbjct_genome_from),
                                                      max(sbjct_genome_to, o_sbjct_genome_to))

    return schema, map_transcript2location

# ------------------------------------------------------------------------


def MapTranscripts2Genes(transcripts, map_transcript2location):
    """map all orthologous and overlapping transcripts into genes.

    The new gene is chosen at random.
    """

    graph = Components.SComponents()

    map_id2info = {}
    added = set()
    for transcript in transcripts:
        map_id2info[transcript.mTranscript] = (
            transcript.mSchema, transcript.mQuality)
        token1, strand1, from1, to1 = map_transcript2location[
            transcript.mTranscript]
        # add link to self as otherwise the component is empty
        graph.add(transcript.mTranscript, transcript.mTranscript)
        # add link to overlapping transcripts
        for x in added:
            token2, strand2, from2, to2 = map_transcript2location[x]

            if token1 == token2 and strand1 == strand2 and \
                    min(to1, to2) - max(from1, from2) > 0:
                graph.add(transcript.mTranscript, x)
        added.add(transcript.mTranscript)

    components = graph.getComponents()

    new_genes = {}
    new_transcripts = []
    for component in components:
        g = component[0]
        new_genes[g] = []
        for id in component:
            s, q = map_id2info[id]
            t = Orthologs.Transcript()
            t.mSchema = s
            t.mTranscript = id
            t.mGene = g
            t.mQuality = q
            new_genes[g].append(t)
            new_transcripts.append(t)

    return new_transcripts, new_genes


def Write(old_transcripts, map_transcript2location, fix, options):
    """map transcripts and write output."""

    if fix:
        new_transcripts = old_transcripts
    else:
        new_transcripts, genes = MapTranscripts2Genes(
            old_transcripts, map_transcript2location)

    if options.write_map:
        map_old = {}
        for x in old_transcripts:
            map_old[x.mTranscript] = x
        for x in new_transcripts:
            options.stdout.write("%s\t%s\n" % (map_old[x.mTranscript], x))
    else:
        options.stdout.write("\n".join(map(str, new_transcripts)) + "\n")


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/orthologs2genes.py 2889 2010-04-07 08:57:57Z andreas $")

    parser.add_option("-s", "--species-regex", dest="species_regex", type="string",
                      help="regular expression to extract species from identifier.")

    parser.add_option("-g", "--gene-regex", dest="gene_regex", type="string",
                      help="regular expression to extract gene from identifier.")

    parser.add_option("-1", "--fix1", dest="fix1", action="store_true",
                      help="do not remap genes for schema1.")

    parser.add_option("-2", "--fix2", dest="fix2", action="store_true",
                      help="do not remap genes for schema2.")

    parser.add_option("-m", "--write-map", dest="write_map", action="store_true",
                      help="write transcript mapping information.")

    parser.add_option("--filename-exons1", dest="filename_exons1", type="string",
                      help="filename with exon information for schema1. If not given, the database will be used.")

    parser.add_option("--filename-exons2", dest="filename_exons2", type="string",
                      help="filename with exon information for schema2. If not given, the database will be used.")

    parser.set_defaults(
        species_regex="^([^|]+)\|",
        gene_regex="^[^|]+\|[^|]+\|([^|]+)\|",
        separator="|",
        tablename_predictions="predictions",
        fix1=False,
        fix2=False,
        write_map=False,
        filename_exons1=None,
        filename_exons2=None,
        schema1=None,
        schema2=None,
    )

    (options, args) = E.Start(parser,
                              add_pipe_options=True,
                              add_database_options=True)

    if not options.filename_exons1 and not options.filename_exons2:
        dbhandle = pgdb.connect(options.psql_connection)

    rs = re.compile(options.species_regex)
    rg = re.compile(options.gene_regex)

    t0 = time.time()

    ninput, noutput, nmissed, nskipped = 0, 0, 0, 0

    orthologs = Orthologs.ReadInterpretation(sys.stdin,
                                             options.separator)

    t1 = time.time()

    E.info("read %i groups in %i seconds" % (len(orthologs), (t1 - t0)))

    if len(orthologs) == 0:
        raise IOError("empty input.")

    orthologs = Orthologs.ClusterOrthologsByGenes(orthologs)

    # map transcripts2genes
    transcripts1, transcripts2, genes1, genes2, weight = orthologs[0]
    schema1 = transcripts1[0].mSchema
    schema2 = transcripts2[0].mSchema

    if options.filename_exons1:
        rschema1, map_transcript2location1 = GetPrediction2LocationFromFile(
            open(options.filename_exons1, "r"), options)
    else:
        map_transcript2location1 = GetPrediction2Location(dbhandle, schema1,
                                                          tablename_predictions=options.tablename_predictions)

    E.info("collected %i for map_transcript2location1" %
           len(map_transcript2location1))

    if options.filename_exons2:
        rschema2, map_transcript2location2 = GetPrediction2LocationFromFile(
            open(options.filename_exons2, "r"), options)
    else:
        map_transcript2location2 = GetPrediction2Location(dbhandle, schema2,
                                                          tablename_predictions=options.tablename_predictions)

    E.info("collected %i for map_transcript2location2" %
           len(map_transcript2location2))

    # swop maps, if schemas have become confused (can happen if files are
    # supplied)
    if rschema2 == schema1 and rschema1 == schema2:
        map_transcript2location1, map_transcript2location2 = map_transcript2location2, map_transcript2location1
    elif schema1 != rschema1 or schema2 != rschema2:
        raise ValueError("mismatch in schemas of orthologs and supplied exons: %s != %s or %s != %s" %
                         (schema1, rschema1, schema2, rschema2))

    t2 = time.time()

    E.info("retrieved locations in %i seconds" % (t2 - t1))

    for transcripts1, transcripts2, genes1, genes2, weight in orthologs:
        ninput += 1
        Write(transcripts1, map_transcript2location1, options.fix1, options)
        Write(transcripts2, map_transcript2location2, options.fix2, options)
        noutput += 1

    E.info("ninput=%i, noutput=%i, nmissed=%i, skipped=%i" %
           (ninput, noutput, nmissed, nskipped))
    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
