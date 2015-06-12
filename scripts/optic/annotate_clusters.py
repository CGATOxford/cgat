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
optic/annotate_clusters.py - 
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

   python optic/annotate_clusters.py --help

Type::

   python optic/annotate_clusters.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import CGAT.IOTools as IOTools
import CGAT.Experiment as E

USAGE = """python optic/annotate_clusters.py [OPTIONS] < clusters > output

annotate a selection of clusters with external information such as INTERPRO or PFAM.

input:

stdin: a list of clusters to analyze. If empty, all clustes are taken.

filename-map: a map of gene_ids to cluster

--filename-[xx]: filename with annotation. Currently implemented
are PFAM and INTERPRO.

PFAM format: ignored, gene_id, ignored, name, short_description, long_description

INTERPRO format: ignored, gene_id, ignored, name, short_description, long_description
"""


class Annotation:

    def __init__(self, identifier, short, description):
        self.mIdentifier = identifier
        self.mShort = short
        self.mDescription = description


def readAnnotationInterpro(infile):
    map_id2annotation = {}

    for line in infile:
        if line[0] == "#":
            continue
        data = line[:-1].split("\t")
        id, identifier, short, description = data[1], data[3], data[4], data[5]
        if short == "NULL":
            continue

        if id not in map_id2annotation:
            map_id2annotation[id] = []

        map_id2annotation[id].append(
            Annotation(identifier, short, description))

    return map_id2annotation


def readAnnotationPfam(infile):
    map_id2annotation = {}

    for line in infile:
        if line[0] == "#":
            continue
        data = line[:-1].split("\t")
        if data[3] == "NULL":
            continue
        if len(data) != 6:
            continue

        id, identifier, short, description = data[1], data[3], data[4], data[5]

        if id not in map_id2annotation:
            map_id2annotation[id] = []

        map_id2annotation[id].append(
            Annotation(identifier, short, description))

    return map_id2annotation


# ------------------------------------------------------------------------
def printAnnotations(outfile, annotations, options):
    """print annotations."""

    # removed redundant identifiers
    identifiers = annotations.keys()
    identifiers.sort()
    descriptions = []
    for id in identifiers:
        descriptions.append(annotations[id].mDescription)

    outfile.write("\t%s\t%s" % (options.separator_fields.join(identifiers),
                                options.separator_fields.join(descriptions)))

# ------------------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/annotate_clusters.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-r", "--species-regex", dest="species_regex", type="string",
                      help="regular expression to extractspecies from identifier.")

    parser.add_option("--map-tsv-file", dest="filename_map_id2cluster", type="string",
                      help="filename with mapping information from id to cluster.")

    parser.add_option("--filename-interpro", dest="filename_interpro", type="string",
                      help="filename with interpro domain information.")

    parser.add_option("--filename-pfam", dest="filename_pfam", type="string",
                      help="filename with pfam domain information.")

    parser.set_defaults(
        master_species="dmel_vs_dmel4",
        separator="|",
        filename_map_id2cluster="input.map",
        filename_interpro="/home/andreas/projects/flies/data_1v5/interpro.list",
        filename_pfam="/home/andreas/projects/flies/data_1v5/pfam.list",
        write_no_annotation=True,
        separator_fields=";",
    )

    (options, args) = E.Start(
        parser, add_database_options=True, add_csv_options=True)

    clusters, nerrors = IOTools.ReadList(sys.stdin)

    map_id2cluster, map_cluster2id = IOTools.ReadMap(open(options.filename_map_id2cluster, "r"),
                                                     both_directions=True)

    if len(clusters) == 0:
        clusters = map_cluster2id.keys()
        clusters.sort()

    if options.filename_interpro:
        map_id2interpro = readAnnotationInterpro(
            open(options.filename_interpro, "r"))

    if options.filename_pfam:
        map_id2pfam = readAnnotationPfam(open(options.filename_pfam, "r"))

    ninput, noutput, nnomaster, nnoannotation = 0, 0, 0, 0
    nskipped = 0

    options.stdout.write("cluster\tgenes")

    if map_id2interpro:
        options.stdout.write("\tinterpro\tidescription")
    if map_id2pfam:
        options.stdout.write("\tpfam\tpdescription")
    options.stdout.write("\n")

    for cluster in clusters:

        ninput += 1
        if cluster not in map_cluster2id:
            if options.loglevel >= 1:
                options.stdlog.write("# cluster %s not in map.\n" % cluster)
            nskipped += 1
            continue

        genes = set()

        for id in map_cluster2id[cluster]:

            s, t, g, q = id.split(options.separator)

            if s != options.master_species:
                continue

            genes.add(g)

        if not genes:
            nnomaster += 1
            continue

        annotations_interpro = {}
        if map_id2interpro:
            for gene in genes:
                if gene in map_id2interpro:
                    for annotation in map_id2interpro[gene]:
                        annotations_interpro[
                            annotation.mIdentifier] = annotation

        annotations_pfam = {}

        if map_id2pfam:
            for gene in genes:
                if gene in map_id2pfam:
                    for annotation in map_id2pfam[gene]:
                        annotations_pfam[annotation.mIdentifier] = annotation

        nannotations = max(len(annotations_pfam), len(annotations_interpro))

        if nannotations == 0 and not options.write_no_annotation:
            nnoannotation += 1
            continue

        options.stdout.write("%s\t%s" % (cluster,
                                         ";".join(genes)))

        if map_id2interpro:
            printAnnotations(options.stdout, annotations_interpro, options)

        if map_id2pfam:
            printAnnotations(options.stdout, annotations_pfam, options)

        options.stdout.write("\n")

        noutput += 1

    if options.loglevel >= 1:
        options.stdlog.write("# ninput=%i, noutput=%i, nskipped=%i, nnomaster=%i, nnoannotation=%i\n" % (
            ninput, noutput, nskipped, nnomaster, nnoannotation))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
