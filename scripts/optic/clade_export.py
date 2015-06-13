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
optic/clade_export.py - 
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

   python optic/clade_export.py --help

Type::

   python optic/clade_export.py --help

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
import time
import sets
import optparse
import math
import tempfile
import glob
import shutil

import CGAT.Experiment as E

import pgdb

""" program $Id: optic/clade_export.py 2781 2009-09-10 11:33:14Z andreas $

Export data for a clade.

This script gets its data from the postgres database and the export_clustering.dir
directory.
"""

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# -------------------------------------------------------------------


def writeReadme(filename, readme, options):
    """write a readme file."""

    dirname = os.path.dirname(filename)
    if dirname and not os.path.exists(dirname):
        os.makedirs(dirname)

    outfile = open(filename, "w")

    outfile.write("Project: %s\n" % options.project_name)
    outfile.write("Release: %s\n\n" % options.release)

    outfile.write(readme + "\n")

    outfile.write(
        "\nAndreas Heger and Chris Ponting, MRC Functional Genetics Unit, Oxford, UK.\n")
    outfile.write("\nCreated on %s.\n\n" %
                  (time.asctime(time.localtime(time.time()))))

    outfile.close()

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# -------------------------------------------------------------------


def exportToFile(dbhandle, statement, filename, options):
    """export a table to a file."""

    dirname = os.path.dirname(filename)
    if dirname and not os.path.exists(dirname):
        os.makedirs(dirname)

    outfile = open(filename, "w")

    cc = dbhandle.cursor()
    cc.execute(statement)

    # write column headers from description tuple
    outfile.write("\t".join(map(lambda x: x[0], cc.description)) + "\n")

    for line in cc.fetchall():
        outfile.write("\t".join(map(str, line)) + "\n")

    cc.close()
    outfile.close()

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# -------------------------------------------------------------------


def exportMalis(dbhandle, options):
    """export multiple alignments.

    Multiple alignments are exported as a tar/gzipped file.
    """

    # export directories

    statement_genes = \
        """SELECT g.group_id AS group_id,
    m.schema || '%s' || m.gene_id AS gene,
    alignment AS alignment
    FROM
    %s.%%s AS m, %s.%s AS g
    WHERE g.cluster_id = m.cluster_id AND 
    g.schema = m.schema AND 
    g.gene_id = m.gene_id 
    ORDER BY group_id    
    """ % (options.separator,
           options.schema,
           options.schema,
           options.table_name_groups_members)

    # do NOT join on prediction id between groups_members and malis,
    # as the prediction_id is not set correctly in the former.
    statement_transcripts = \
        """SELECT g.group_id AS group_id,
    m.schema || '%s' || m.prediction_id || '%s' || m.gene_id || '%s' || m.class AS transcript,
    alignment AS alignment
    FROM
    %s.%%s AS m, %s.%s AS g
    WHERE g.cluster_id = m.cluster_id AND
    g.schema = m.schema AND 
    g.gene_id = m.gene_id 
    ORDER BY group_id
    """ % (options.separator, options.separator, options.separator,
           options.schema,
           options.schema,
           options.table_name_groups_members)

    tables = ((statement_genes % "malis_genes_aa", "genes_aa"),
              (statement_genes % "malis_genes_na", "genes_na"),
              (statement_transcripts %
               "malis_transcripts_aa", "transcripts_aa"),
              (statement_transcripts % "malis_transcripts_na", "transcripts_na"))

    for statement, filename in tables:
        exportToFile(dbhandle,
                     statement,
                     options.build_dir + "/multiple_alignments/%s" % filename,
                     options)

    readme = """Multiple alignments of orthologous groups.

This directory contains multiple alignments of orthologous groups. Multiple alignments
were build from concatenated exons of each gene using the translated sequence.

The nucleotide alignments were built by threading the nucleotide sequence back onto
the aligned peptide sequences.

Stop-codons have been masked and frameshifts have been removed, thus there can be
differences between the cds sequence and the aligned sequence.

File            Contents
genes_aa        amino acid alignments of concatenated exons, one entry per gene
genes_na        nucelotide alignments of concatenated exons, one entry per gene
transcripts_aa  amino acid alignments of transcripts, multiple entries per gene possible
transcripts_na  nucleotide alignments of transcripts, multiple entries per gene possible
"""

    writeReadme(options.build_dir + "/multiple_alignments/readme",
                readme,
                options)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# -------------------------------------------------------------------


def exportOrthologs(dbhandle, options):
    """export orthology information.
    """

    statement = """SELECT
    group_id AS group_id,
    schema || '%s' || gene_id AS gene
    FROM %s.%s
    ORDER BY group_id
    """ % (options.separator, options.schema, options.table_name_groups_members )

    exportToFile(dbhandle,
                 statement,
                 options.build_dir + "/orthologs/groups",
                 options)

    statement = """SELECT
    m.set_id AS set_id,
    m.schema || '%s' || m.gene_id AS gene
    FROM %s.%s AS m, %s.%s AS s
    WHERE s.set_id = m.set_id AND POSITION('0' in s.pattern) = 0
    ORDER BY m.set_id
    """ % (options.separator,
           options.schema, options.table_name_sets_members,
           options.schema, options.table_name_sets)

    exportToFile(dbhandle,
                 statement,
                 options.build_dir + "/orthologs/sets",
                 options)

    def getTrees(table_name, filename):

        statement = """SELECT
        '>' || group_id || E'\n' || nh
        FROM %s.%s
        ORDER BY group_id
        """ % (options.schema, table_name)

        exportToFile(dbhandle,
                     statement,
                     filename,
                     options)

    getTrees(options.table_name_trees_njtree,
             options.build_dir + "/orthologs/nj_trees")

    getTrees(options.table_name_trees_njtree,
             options.build_dir + "/orthologs/ds_trees")

    readme = """Orthologous groups.

Pairwise orthology relationships are grouped into clusters of orthologs across all
species in the clade.

File            Contents
groups          members of orthologous groups
sets            members of strict ortholog sets that include just one gene per species
nj_trees        trees build with njtree for orthologous groups
ds_trees        trees with ds branch lengths estimated with PAML
"""

    writeReadme(options.build_dir + "/orthologs/readme",
                readme,
                options)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# -------------------------------------------------------------------


def exportTrees(dbhandle, options):
    """export tree information.
    """

    readme = """Phylogenetic trees of orthologous groups.

Pairwise orthology relationships are grouped into clusters of orthologs across all
species in the clade.

File            Contents
groups          members of orthologous groups
sets            members of strict ortholog sets that include just one gene per species
"""

    writeReadme(options.build_dir + "/orthologs/readme",
                readme,
                options)


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# -------------------------------------------------------------------
def exportGenes(dbhandle, options):
    """export gene sets used in this analysis.
    """

    files = glob.glob( options.dir_export_predictions + "/export_clustering_*.fasta" ) +\
        glob.glob(
            options.dir_export_predictions + "/export_clustering_*.exons")

    if options.loglevel >= 1:
        options.stdlog.write("# exporting genes: %i files\n" % (len(files)))

    if len(files) == 0:
        return

    nfiles = len(files)
    ndirs_ok = 0
    nfiles_ok = 0

    dirname = options.build_dir + "/genes"

    if dirname and not os.path.exists(dirname):
        os.makedirs(dirname)

    for file_src in files:

        filename = os.path.basename(file_src)
        file_dest = dirname + "/" + filename

        if options.loglevel >= 3:
            options.stdlog.write(
                "# processing %s: %s\n" % (file_src, file_dest))

        shutil.copyfile(file_src, file_dest)
        nfiles_ok += 1

    if options.loglevel >= 1:
        options.stdlog.write("section %s - summary\n" % section)
        options.stdlog.write("files:\t\t%6i out of %6i (%6.4f%%)\n" % (
            nfiles_ok, nfiles, 100 * float(nfiles_ok) / nfiles))

    readme = """Gene predictions.

Contents:

xxx_peptides.fasta:        peptide sequences of predicted genes.
xxx_cds.fasta:             nucleotide sequences of predicted genes.
xxx.exons:                 exons of predicted genes.

The exons file is a tab-separated format.

Column  Content
1       Prediction identifier
2       Contig identifier
3       Strand
4       Phase
5       exon_id
6       location of exon on cds sequence
7       location of exon on contig
"""
    writeReadme(options.build_dir + "/genes/readme",
                readme,
                options)

    return


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/clade_export.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-r", "--release", dest="release", type="string",
                      help="release.")

    parser.add_option("-p", "--project", dest="project_name", type="string",
                      help="project.")

    parser.add_option("-e", "--sections", dest="sections", type="string",
                      help="sections to export: malis,orthologs,genes")

    parser.add_option("-a", "--package-name", dest="package_name", type="string",
                      help="Package name. This also provides the build directory.")

    parser.add_option("-s", "--schema", dest="schema", type="string",
                      help="schema to export from.")

    parser.add_option("-c", "--compress", dest="compress", action="store_true",
                      help="Compress files.")

    parser.set_defaults(
        sections="malis,orthologs,genes",
        dir_export_predictions="../export/export_clustering.dir",
        build_dir=None,
        compress=False,
        schema=None,
        separator="|",
        table_name_groups="groups",
        table_name_groups_members="groups_members",
        table_name_sets="ortholog_sets",
        table_name_sets_members="ortholog_sets_members",
        table_name_trees_njtree="groups_nj_trees",
        table_name_trees_dstree="ds_trees",
        release="unknown",
        project_name="unknown",
    )

    (options, args) = E.Start(parser, add_database_options=True)

    options.build_dir = options.package_name + "_" + options.release
    options.sections = options.sections.split(",")
    dbhandle = pgdb.connect(options.psql_connection)

    if not os.path.exists(options.build_dir):
        os.mkdir(options.build_dir)

    for section in options.sections:

        if section == "malis":
            exportMalis(dbhandle, options)
        elif section == "orthologs":
            exportOrthologs(dbhandle, options)
        elif section == "genes":
            exportGenes(dbhandle, options)

    readme = """This directory contains multiple alignments, orthologs and gene sets.

Director                Contents

orthologs               clusters of orthologous groups and strict ortholog sets.
multiple_alignments     multiple alignments of orthologous groups
genes                   gene sets used in this analysis

Transcript identifiers

The format for transcript identifiers is species|transcript_id|gene_id|quality:

species:        species name (TEXT)
transcript_id:  transcript identifier
gene_id:        gene identifier 
quality:        quality code. See http://wwwfgu.anat.ox.ac.uk:8080/flies/documentation.html
                for more information.

For example, the James Bond transcript would be: mm_hsapiens|007|007|CG

Gene identifiers

Gene identifiers are similar to transript identifiers, except that they just contain species and
gene_id, for example: mm_hsapiens|007

I am aware that '|' is an unfortunate choice for a separator (I was looking for something readable,
that could be used in trees, belvu alignments, and others, which ruled out ':' and '/' and others).
I might clean up identifiers in the future.
"""

    writeReadme(options.build_dir + "/readme",
                readme,
                options)

    if options.compress:
        if options.loglevel >= 1:
            options.stdlog.write("# compressing files started.\n")
            options.stdlog.flush()
        os.system(
            'find %s -name "*" -not -name "readme" -exec gzip -q {} \;' % options.build_dir)
        if options.loglevel >= 1:
            options.stdlog.write("# compressing files finished.\n")
            options.stdlog.flush()

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
