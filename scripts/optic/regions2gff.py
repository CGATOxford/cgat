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
optic/regions2gff.py - 
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

   python optic/regions2gff.py --help

Type::

   python optic/regions2gff.py --help

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
import CGAT.GenomicIO as GenomicIO
import CGAT.IOTools as IOTools
import pgdb
import CGAT.Regions as Regions
import CGAT.GTF as GTF
import CGAT.Intervals as Intervals

USAGE = """python %s [OPTIONS] 

Version: $Id: optic/regions2gff.py 2781 2009-09-10 11:33:14Z andreas $

extract features from the genomic database and save as a gff file. 

Various filters are available:
        * genes without PHYOP orthologs
        * genes with PHYOP orthologs
        * genes without OPTIC orthologs
        * genes with OPTIC orthologs

Features are:
        * mRNAs = regions covered by a transcript
        * Exons
        * Introns
        * Intergenic regions
        * Exonic regions
        * Intronic regions
        * Genes = regions covered by a genes

Additional data fields:
        * Pairwise dS, dN and dN/dS (PHYOP)
        * Lineage specific dS, dN and dN/dS (OPTIC)

Note that the dS, dN, dN/dS values are not calculated
for each feature, but simply taken from the complete gene.
""" % sys.argv[0]


def getRestrictionOnPredictionIds(schema, prediction_ids="all"):

    xfrom = ""

    if prediction_ids == "all":
        xwhere = ""
    elif prediction_ids == "nr":
        xfrom = "%s.redundant AS r," % schema
        xwhere = "AND p.prediction_id = r.rep_prediction_id AND r.rep_prediction_id = r.mem_prediction_id"
    else:
        xwhere = "AND p.prediction_id  IN ('%s')" % "','".join(
            map(str, prediction_ids))

    return xfrom, xwhere


def convertCoordinates(old, contigs):
    """convert coordinates to +strand coordinates."""

    new = []
    for gene_id, prediction_id, contig, strand, start, end, frame in old:
        start, end = Genomics.ToForwardCoordinates(start,
                                                   end,
                                                   strand,
                                                   contigs[contig])

        new.append((gene_id, prediction_id, contig, strand, start, end, frame))

    return new


def getMRNAs(dbhandle, schema, options, prediction_ids, contigs):
    """get a list of mRNAs.

    mRNAs extend from the first residue of the first exon to
    the last residue of the last exon of a prediction.
    """

    xfrom, xwhere = getRestrictionOnPredictionIds(schema, prediction_ids)

    # select genomic location for all predictions
    # retrieve: gene_id, prediction_id, sbjct_token, sbjct_strand,
    # sbjct_genome_from, sbjct_genome_to
    statement = """
    SELECT g.gene_id, p.prediction_id, p.sbjct_token, p.sbjct_strand, p.sbjct_genome_from, p.sbjct_genome_to, '.'
    FROM %s
    %s.%s AS p, 
    %s.%s AS g 
    WHERE p.prediction_id = g.prediction_id %s
    """ % ( xfrom,
            schema, options.tablename_predictions,
            schema, options.tablename_genes,
            xwhere)

    cc = dbhandle.cursor()
    cc.execute(statement)
    locations = cc.fetchall()
    cc.close()

    return convertCoordinates(locations, contigs)


def getGenes(dbhandle, schema, options, prediction_ids, contigs):
    """get a list of genes

    mRNAs extend from the first residue of the first exon to
    the last residue of the last exon of any transcript in a gene.
    """

    xfrom, xwhere = getRestrictionOnPredictionIds(schema, prediction_ids)

    # select genomic location for all predictions
    # retrieve: gene_id, prediction_id, sbjct_token, sbjct_strand,
    # sbjct_genome_from, sbjct_genome_to
    statement = """
    SELECT g.gene_id, g.gene_id, p.sbjct_token, p.sbjct_strand, MIN(p.sbjct_genome_from), MAX(p.sbjct_genome_to), '.',
    FROM %s
    %s.%s AS p, 
    %s.%s AS g 
    WHERE p.prediction_id = g.prediction_id %s
    GROUP BY g.gene_id, p.sbjct_token, p.sbjct_strand
    """ % ( xfrom,
            schema, options.tablename_predictions,
            schema, options.tablename_genes,
            xwhere)

    cc = dbhandle.cursor()
    cc.execute(statement)
    locations = cc.fetchall()
    cc.close()

    return convertCoordinates(locations, contigs)

##########################################################################


def getExons(dbhandle, schema, options, prediction_ids, contigs):
    """get exons.
    """

    xfrom, xwhere = getRestrictionOnPredictionIds(schema, prediction_ids)

    # select genomic location for all introns
    # filter on genome_exon_to > 0 in order to eliminate
    # all non-matching exons to query.
    # retrieve: prediction_id, sbjct_token, sbjct_strand, sbjct_genome_from,
    # sbjct_genome_to
    statement = """
    SELECT DISTINCT
    g.gene_id, p.prediction_id, sbjct_token, sbjct_strand, genome_exon_from, genome_exon_to, e.exon_frame
    FROM %s
    %s.%s AS p,
    %s.%s AS g,
    %s.%s AS e
    WHERE p.prediction_id = e.prediction_id AND 
    g.prediction_id = p.prediction_id AND 
    genome_exon_to > 0
    %s
    GROUP BY g.gene_id, p.prediction_id, sbjct_token, sbjct_strand, genome_exon_from, genome_exon_to, e.exon_frame
    ORDER BY g.gene_id, p.prediction_id, genome_exon_from
    """ % (xfrom,
           schema, options.tablename_predictions,
           schema, options.tablename_genes,
           schema, options.tablename_exons,
           xwhere)

    cc = dbhandle.cursor()
    cc.execute(statement)
    locations = cc.fetchall()
    cc.close()

    return convertCoordinates(locations, contigs)

##########################################################################


def getIntrons(dbhandle, schema, options, prediction_ids, contigs):
    """get regions corresponding to introns

    The regions are NOT normalized. The same intron might appear twice due
    to two alternative transcripts sharing two adjacent introns. Similary, 
    the introns might contain exons, if one transcript skips an exon from
    another.
    """

    result = getExons(dbhandle, schema, options, prediction_ids, contigs)

    locations = []
    last_to = None
    last_id = None

    # sort by prediction_id and location
    result.sort(lambda x, y: cmp((x[1], x[4]), (y[1], y[4])))

    for gene_id, prediction_id, sbjct_token, sbjct_strand, sbjct_from, sbjct_to, frame in result:
        if last_id != prediction_id:
            last_to = None
            last_id = prediction_id
        if last_to:
            locations.append(
                (gene_id, prediction_id, sbjct_token, sbjct_strand, last_to, sbjct_from, frame))
        last_to = sbjct_to

    return locations

##########################################################################


def getIntronicRegions(dbhandle, schema, options, prediction_ids, contigs):
    """get intronic regions. These are regions that are between exons in
    all alternative transcripts for a gene.
    """

    result = getExons(dbhandle, schema, options, prediction_ids, contigs)

    locations = []
    last_to = None
    last_id = None
    last_contig = None
    last_strand = None
    last_frame = None

    # sort by gene_id
    result.sort()
    regions = []

    def processChunk(gene_id, contig, strand, frame, regions):
        if gene_id is None:
            return

        start = min(map(lambda x: x[0], regions))
        end = max(map(lambda x: x[0], regions))

        intervals = Intervals.complementIntervals(regions, start, end)
        for start, end in intervals:
            locations.append(
                (gene_id, gene_id, contig, strand, start, end, frame))

    for gene_id, prediction_id, contig, strand, start, end, frame in result:
        if last_id != gene_id:
            processChunk(
                last_id, last_contig, last_strand, last_frame, regions)
            last_id = gene_id
            last_contig = contig
            last_strand = strand
            regions = []
        regions.append((start, end))

    processChunk(last_id, last_contig, last_strand, regions)

    return locations

##########################################################################


def getExonicRegions(dbhandle, schema, options, prediction_ids, contigs):
    """get exonic regions. These are regions that are covered by exons.
    """

    result = getExons(dbhandle, schema, options, prediction_ids, contigs)

    locations = []
    last_to = None
    last_id = None
    last_contig = None
    last_strand = None
    last_frame = None
    # sort by gene_id
    result.sort()
    regions = []

    def processChunk(gene_id, contig, strand, frame, regions):
        if gene_id is None:
            return

        for start, end in Intervals.combineIntervals(regions):
            locations.append((gene_id, gene_id, contig, strand, start, end))

    for gene_id, prediction_id, contig, strand, start, end, frame in result:
        if last_id != gene_id:
            processChunk(
                last_id, last_contig, last_strand, last_frame, regions)
            last_id = gene_id
            last_contig = contig
            last_strand = strand
            regions = []
        regions.append((start, end))

    processChunk(last_id, last_contig, last_strand, regions)

    return locations

##########################################################################


def getExonsThirdCodons(dbhandle, schema, options, prediction_ids, contigs):
    """get third codon positions in exons."""

    result = getExons(dbhandle, schema, options, prediction_ids, contigs)

    # sort by prediction_id and location
    result.sort(lambda x, y: cmp((x[1], x[4]), (y[1], y[4])))

    locations = []

    def processChunk(prediction_id, gene_id, contig, strand, regions):

        if gene_id is None:
            return

        # re-arrange positions on negative strand
        if Genomics.IsNegativeStrand(strand):
            # convert to negative strand coordinates counting from 0
            coordinate_offset = max(map(lambda x: x[1], regions))
            regions = map(
                lambda x: (coordinate_offset - x[1], coordinate_offset - x[0]), regions)
            regions.sort()
        else:
            coordinate_offset = 0

        offset = 0
        for start, end in regions:
            start -= offset
            for x in range(start + 2, end, 3):
                if coordinate_offset:
                    # the factor -1 results from the open/closed
                    # bracket notation
                    c = coordinate_offset - x - 1
                else:
                    c = 0
                locations.append(
                    (prediction_id, gene_id, contig, strand, c, c + 1))
            offset = (end - start) % 3

        if (offset != 0):
            if options.loglevel >= 1:
                options.stdlog.write("# WARNING: prediction=%s, gene=%s on %s:%s : frame did not add up\n" % (
                    prediction_id, gene_id, contig, strand))

    regions = []
    last_to = None
    last_prediction_id = None
    last_gene_id = None
    last_contig = None
    last_strand = None

    if options.loglevel >= 2:
        options.stdlog.write("# processing %i entries.\n" % len(result))
        options.stdlog.flush()

    niterations = 0
    for gene_id, prediction_id, contig, strand, start, end in result:
        if last_prediction_id != prediction_id:
            processChunk(
                last_prediction_id, last_gene_id, last_contig, last_strand, regions)
            last_prediction_id = prediction_id
            last_gene_id = gene_id
            last_contig = contig
            last_strand = strand

        regions.append((start, end))

        niterations += 1
        if options.loglevel >= 2 and niterations % options.report_step == 0:
            options.stdlog.write("# iteration: %i (%5.2f%%)\n" % (
                niterations, 100.0 * niterations / len(result)))
            options.stdlog.flush()

    processChunk(
        last_prediction_id, gene_id, last_contig, last_strand, regions)

    return locations

##########################################################################


def getIntergenicRegions(dbhandle, schema, options, prediction_ids, contigs):
    """get intergenic regions. These are regions that are lying between genes.
    """

    result = getMRNAs(dbhandle, schema, options, prediction_ids, contigs)

    locations = []

    # sort by contig
    result.sort(lambda x, y: cmp(x[2], y[2]))
    last_contig = None
    regions = []

    def processChunk(contig, regions):
        if contig is None:
            return

        start = 0
        end = contigs[contig]

        regions = Intervals.combineIntervals(regions)
        for xstart, xend in Intervals.complementIntervals(regions, start, end):
            locations.append(
                ("intergenic", "intergenic", contig, "+", xstart, xend, "."))

    for gene_id, prediction_id, contig, strand, start, end, frame in result:
        if contig != last_contig:
            processChunk(last_contig, regions)
            regions = []
            last_contig = contig

        regions.append((start, end))

    processChunk(last_contig, regions)

    return locations

##########################################################################


def getLocations(dbhandle, schema, options, prediction_ids, contigs, taboo_regions=None):
    """retrieve locations for a set of gene ids."""

    t_start = time.time()

    if options.loglevel >= 1:
        options.stdlog.write(
            "# %s: retrieving candidate regions.\n" % (schema))
        options.stdlog.flush()

    if options.type == "mrnas":
        locations = getMRNAs(
            dbhandle, schema, options, prediction_ids, contigs)
    elif options.type == "genes":
        locations = getGenes(
            dbhandle, schema, options, prediction_ids, contigs)
    elif options.type == "introns":
        locations = getIntrons(
            dbhandle, schema, options, prediction_ids, contigs)
    elif options.type == "intronic":
        locations = getIntronicRegions(
            dbhandle, schema, options, prediction_ids, contigs)
    elif options.type == "exons":
        locations = getExons(
            dbhandle, schema, options, prediction_ids, contigs)
    elif options.type == "exonic":
        locations = getExonicRegions(
            dbhandle, schema, options, prediction_ids, contigs)
    elif options.type == "exons-third-codons":
        locations = getExonsThirdCodons(
            dbhandle, schema, options, prediction_ids, contigs)
    elif options.type == "intergenic":
        locations = getIntergenicRegions(
            dbhandle, schema, options, prediction_ids, contigs)
    else:
        raise "unknown location type %s" % options.type

    if options.loglevel >= 1:
        options.stdlog.write("# %s: retrieved %i candidate regions in %i seconds.\n" % (
            schema, len(locations), time.time() - t_start))
        options.stdlog.flush()

    return locations

##########################################################################


def getMapFeature2Property(options):
    """get a map of features to a property."""

    if options.filename_properties:
        map_feature2property = {}
        infile = open(options.filename_properties, "r")
        for line in infile:
            if line[0] == "#":
                continue
            feature, property = line[:-1].split("\t")[:2]
            map_feature2property[feature] = property
    else:
        map_feature2property = None

    return map_feature2property

##########################################################################


def processPredictions(dbhandle, schema, options, prediction_ids, taboo_regions=None, map_feature2property=None):
    """print required regions for a set of prediction_ids from a given schema."""

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

    locations = getLocations(
        dbhandle, schema, options, prediction_ids, contig_size, taboo_regions)

    contigs = set(map(lambda x: x[2], locations))

    if options.loglevel >= 1:
        options.stdlog.write("# %s: %i locations are on %i contigs.\n" %
                             (schema, len(locations), len(contigs)))
        options.stdlog.flush()

    if options.loglevel >= 3:
        for location in locations:
            options.stdlog.write("# %s\n" % (str(location)))

    if len(locations) == 0:
        return

    # sort locations by contig, that way sequence retrieval is faster
    # sort also by identifier
    locations.sort(lambda x, y: cmp((x[2], x[0]), (y[2], y[0])))
    tokens = {}
    noutput = 0
    nskipped_length = 0
    nskipped_overlap = 0
    nskipped_doubles = 0
    nskipped_property = 0

    output_regions = set()
    prediction_ids = set()
    gene_ids = set()

    for gene_id, prediction_id, token, strand, start, end, frame in locations:

        if taboo_regions:
            overlaps = taboo_regions.getOverlaps(token, strand, start, end)

            if overlaps:
                if options.loglevel >= 2:
                    options.stdlog.write(
                        "# %s eliminated due to overlap with taboo region.\n" % (token))
                nskipped_overlap += 1
                continue

        # filter output for redundant entries
        k = "%s-%s-%i-%i" % (token, strand, start, end)
        if k in output_regions:
            nskipped_doubles += 1
            continue

        output_regions.add(k)

        gff = GTF.Entry()
        gff.name = token
        gff.stand = strand
        gff.start = start
        gff.strand = strand
        gff.end = end
        gff.feature = options.type
        gff.frame = frame

        if map_feature2property:
            found = True
            if prediction_id in map_feature2property:
                value = map_feature2property[prediction_id]
            elif gene_id in map_feature2property:
                value = map_feature2property[gene_id]
            else:
                value = "."
                found = False

            # decide which ones to keep
            if (not options.invert_property and not found) or \
                    (options.invert_property and found):
                nskipped_property += 1
                continue

            gff.score = value
        else:
            # the default score is the segment length
            gff.score = gff.end - gff.start

        gff.addAttribute("gene", gene_id)
        gff.addAttribute("transcript", prediction_id)

        gene_ids.add(gene_id)
        prediction_ids.add(prediction_id)

        options.stdout.write(str(gff) + "\n")
        options.stdout.flush()
        noutput += 1

    if options.loglevel >= 1:
        if not prediction_ids:
            prediction_ids = []

        options.stdlog.write("# %s: ninput=%i, npredictions=%i, ngenes=%i, noutput=%i, nskipped_property=%i, nskipped_length=%i, nskipped_overlap=%i, nskipped_doubles=%i\n" % (
            schema, len(locations), len(prediction_ids), len(gene_ids), noutput, nskipped_property, nskipped_length, nskipped_overlap, nskipped_doubles))


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/regions2gff.py 2781 2009-09-10 11:33:14Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="pattern to look for sequence filename.")

    parser.add_option("-i", "--ids", dest="ids", type="string",
                      help="comma separated list of prediction ids. Use 'all' to use all predictions.")

    parser.add_option("-f", "--filename-ids", dest="filename_ids", type="string",
                      help="filename with prediction ids.")

    parser.add_option("-t", "--sequence-type", dest="type", type="choice",
                      choices=("genes", "mrnas", "introns", "intronic",
                               "exons", "exonic", "intergenic", "exons-third-codons"),
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

    parser.add_option("--filename-properties", dest="filename_properties", type="string",
                      help="filename with mapping information between features and properties.")

    parser.add_option("--invert-properties", dest="invert-properties", action="store_true",
                      help="instead of printing features which have properties, print those that have not.")

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
        tablename_predictions="predictions",
        tablename_exons="exons",
        tablename_genes="genes",
        tablename_quality="quality",
        schema=None,
        output_format="fasta",
        fasta_format="id-coordinates",
        type="mrnas",
        min_length=1,
        id_format="id",
        mmultiple=False,
        separator="|",
        filename_taboo_regions=False,
        output_coordinate_format="full",
        filename_properties=None,
        invert_property=False,
        report_step=10000
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

    map_feature2property = getMapFeature2Property(options)

    processPredictions(dbhandle, options.schema, options,
                       prediction_ids, taboo_regions, map_feature2property)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
