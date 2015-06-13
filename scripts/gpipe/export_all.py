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
gpipe/export_all.py - 
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

   python gpipe/export_all.py --help

Type::

   python gpipe/export_all.py --help

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

""" program $Id: gpipe/export_all.py 2781 2009-09-10 11:33:14Z andreas $

prepare export of fly data.
"""


def ExportFilesFromSubdirs(directories,
                           build_dir,
                           file_list,
                           options,
                           section="",
                           extract_id_from_dirname=None,
                           extract_id_from_filename=None,
                           within_subdirs=None,
                           readme=None):

    missing = []
    nfiles = 0
    ndirs_ok = 0
    nfiles_ok = 0

    if options.loglevel >= 1:
        options.stdlog.write(
            "# exporting section %s: %i directories\n" % (section, len(directories)))

    if os.path.exists(build_dir):
        os.rmdir(build_dir)

    os.mkdir(build_dir)

    for src_dir in directories:
        if extract_id_from_dirname:
            dir_id = extract_id_from_dirname.search(src_dir).groups()[0]
        else:
            dir_id = None

        has_error = False

        if within_subdirs:
            dest_dir = build_dir + "/" + within_subdirs % dir_id
            os.mkdir(dest_dir)
        else:
            dest_dir = build_dir

        for file_src, file_dest, converter in file_list:

            if options.loglevel >= 3:
                options.stdlog.write(
                    "# processing %s: %s\n" % (file_src, file_dest))

            # if we have not id and there is a wildcard in file_src: use glob
            # operator to get all files
            if not dir_id and re.search("%s", file_src):
                x = re.sub("%s", "*", file_src)
                files = glob.glob(src_dir + "/" + x)
                if options.stdlog >= 4:
                    options.stdlog.write(
                        "# found %i files using pattern %s\n" % (len(files), src_dir + "/" + x))
                single_file = False
            else:
                single_file = True
                if re.search("%s", file_src):
                    file_src = file_src % dir_id
                files = (src_dir + "/" + file_src, )

            if not file_dest:
                file_dest = dest_dir + "/" + os.path.basename(file_src)
            else:
                file_dest = dest_dir + "/" + file_dest

            for f in files:
                nfiles += 1

                if not dir_id:
                    if extract_id_from_filename:
                        id = extract_id_from_filename.search(f).groups()[0]
                    else:
                        id = ""
                else:
                    id = dir_id

                if single_file and re.search("%s", file_src):
                    file_src = file_src % id
                else:
                    file_src = f

                if re.search("%s", file_dest):
                    fdest = file_dest % id
                else:
                    fdest = file_dest

                if not os.path.exists(file_src):
                    missing.append(file_src)
                    has_error = True
                    continue

                if options.loglevel >= 5:
                    options.stdlog.write(
                        "# copying %s to %s\n" % (file_src, fdest))

                if converter:
                    os.system(converter % {'src': file_src, 'dest': fdest})
                else:
                    shutil.copyfile(file_src, fdest)
                nfiles_ok += 1

                id = None

        if not has_error:
            ndirs_ok += 1

    if options.loglevel >= 1:
        options.stdlog.write("section %s - summary\n" % section)
        options.stdlog.write("directories:\t%6i out of %6i (%6.4f%%)\n" % (
            ndirs_ok, len(directories), 100 * float(ndirs_ok) / len(directories)))
        options.stdlog.write("files:\t\t%6i out of %6i (%6.4f%%)\n" % (
            nfiles_ok, nfiles, 100 * float(nfiles_ok) / nfiles))

    if options.loglevel >= 2:
        options.stderr.write(
            "# section %s - %i missing files\n" % (section, len(missing)))
        for m in missing:
            options.stderr.write("# %s\n" % (m))

    outfile = open(build_dir + "/readme", "w")

    outfile.write(
        "Protein coding gene predictions and analysis of gene predictions of 12 fly genomes.\n")
    outfile.write(
        "\nAndreas Heger and Chris Ponting, MRC Functional Genetics Unit, Oxford, UK.\n")
    outfile.write("\nThis is release %s.\n" % options.release)
    outfile.write("\nSection %s created on %s.\n\n" %
                  (section, time.asctime(time.localtime(time.time()))))

    if readme:
        outfile.write(readme)

    outfile.close()

    return


def ExportMalis(options):
    """export multiple alignments.

    Multiple alignments are exported as a tar/gzipped file.
    """

    # export directories
    directories = glob.glob(options.dir_export_malis + "/data.dir/*.dir")

    file_list = (("cluster_%s.aa_mali", None, None, ),
                 ("cluster_%s.aa_mali", None, None, ),
                 ("cluster_%s.raw_mali", None, None,),
                 ("cluster_%s_ks.tree", None, None, ),
                 ("cluster_%s.bl_kaks", "cluster_%s.kaks",
                  "cut -f  1-12, %(src)s > %(dest)s"),
                 ("cluster_%s.bs_partitions", None, None),
                 ("cluster_%s.bs_evaluate.components", None, None),
                 ("cluster_%s.bs_evaluate.inconsistencies", None, None),
                 ("cluster_%s.bs_evaluate.subtrees", None, None),
                 )

    readme = """Contents: Multiple alignments of orthologs.

Input is the cds of predicted orthologs. Clusters are given by the multiple orthology
procedure.

.raw_mali: multiple alignment created by dialign (< 50 sequences) or muscle (<500 sequences).
           Lower case characters correspond to unaligned regions.
.na_mali : multiple alignment with frame-shifts and stop-codons removed.
           Unaligned regions have been removed and exon/intron transitions are marked by
           case changes in the alignment.
.aa_mali : translated .aa_mali multiple alignment
.bl_mali : multiple alignment after Gblocks clean up (badly aligned columns are removed).
.kaks    : pairwise dn and ds values for multiply aligned sequences (using the .bl_mali as input).
.tree    : phylogenetic tree in New Hampshire format. The tree was estimated by kitsch from the
           PHYLIP package using a ks distance matrix.

All the .xx_mali files are simple fasta format.

Dn and ds is claculated using codonml from the PAML package. The format of the .kaks file is a tab-separated format with the
columns corresponding to:

seq1    seq2    dN      dS      dN/dS   N       S       dN_std_err      dS_std_err      kappa   lnL     tau     error_str
"""

    ExportFilesFromSubdirs(directories,
                           options.build_dir + "/multiple_alignments",
                           file_list,
                           options,
                           extract_id_from_dirname=re.compile(
                               "cluster_(\S+).dir"),
                           section="multiple alignments",
                           within_subdirs="cluster_%s",
                           readme=readme)


def ExportOrthologs(options):
    """export orthology information."""

    file_list = (("orthologs_consistent.components.map", "map_sequence2cluster", None, ),
                 ("orthologs_consistent.orgs_per_cluster",
                  "clusters.info", None, ),
                 ("orthologs_consistent.patterns", "patterns", None, ),
                 ("orthologs_consistent.summary", "summary", None),
                 )

    readme = """Contents: Clusters of orthologs

Pairwise orthology relationships are grouped into clusters of orthologs across all
fly genomes. The method is clustering by components and then filtering by edge and
vertex consistency (using neighbourhood overlap in the case of edges and the clustering
coefficient for vertices).

map_transcript2cluster: map of transcripts to clusters.

cluster_info:           number of genes/transcripts per cluster from each species.

patterns:               species patterns and whether the absences of species is
                        consistent with the phylogenetic tree.
                        
summary:                summary of clustering.
"""

    ExportFilesFromSubdirs((options.dir_export_orthologs, ),
                           options.build_dir + "/orthologs",
                           file_list,
                           options,
                           section="orthologs",
                           within_subdirs=None,
                           readme=readme)


def ExportCodonbias(options):
    """export codon bias information."""

    directories = glob.glob(options.dir_export_codonbias + "/bias_*")

    file_list = (("all.data", "%s.data", None, ),
                 )

    readme = """Sequence properties for predicted sequences.

For each predicted transcript in each species, various sequence features are listed. The
data is organized in tab separated tables. Some of the information in the file:

Header          Content

GENENAME        name of the transcripts
is_selected     whether or not this sequence was used in computing the codon usage table
                (only "good quality" sequenes were used)
is_dominant     whether this sequence is part of the dominant set as calculate the method Carbone et al. (2005)
CAI             Codon adaptive index, using species specific codon weights calculated by the method of Carbone et al. (2005).
CAIREF          Codon adaptive index, calculated using codon weights from D. melanogaster.
CAIBIAS         Codon adaptive index using species specific codon weights derived ribosomal protein sequences.
ENC             Effective number of codons, calculated with codonW.
length          sequence length in terms of nucleotides
codons          sequence length in terms of codons
nA-nY           absolute frequency of amino acid X in sequence
pA-pY           relative frequence of amino acid X in sequence
nstops          number of stop codons
nsites1d        number of non-degenerate sites
nsites2d        number of two-fold degenerate sites
nsites3d        number of three-fold degenerate sites
nsites4d        number of four-fold degenerate sites
ngc             number of G+Cs
ngc3            number of G+Cs in third codon position
n4gc3           number of G+Cs in four-fold degenerate third codon positions
pgc             percent G+C of sequence
pgc3            percent G+C in third codon positions
p4gc3           percent G+Cs in four-fold degenerate third codon positionsp4gc3

+ various information theoretic measures comparing codon usage versus different
kinds of reference codon usage:

0: dominant set (dominant set as derived from the method of Carbone et al (2005).
1: non-dominant set (all sequences not in the dominant set)
2: bulk set (dominant set + non-dominant set)
3: uniform frequencies
4: ribosomal protein set (codon usage in ribosomal proteins).

ml:     message length
rel_ml: relative message length
kl:     symmetric Kullback-Leibler divergence
"""

    ExportFilesFromSubdirs(directories,
                           options.build_dir + "/codonbias",
                           file_list,
                           options,
                           extract_id_from_dirname=re.compile("bias_(\S+)"),
                           section="codon bias",
                           within_subdirs=None,
                           readme=readme)

    return


def ExportAAA(options):
    """export gene predictions for AAA."""

    directories = glob.glob(options.dir_export_aaa, )

    file_list = (("%s_caf1_genes.map", None, None, ),
                 ("%s_caf1_predictions.map", None, None, ),
                 ("%s_caf1.gff", None, None, ),
                 )

    readme = """Gene predictions as submitted to AAA. The files ending in .map
map transcript/gene identifiers to the identfiers required by AAA. See
http://rana.lbl.gov/drosophila/wiki/index.php/Formats_and_Naming_Schemes
for more information.
"""

    ExportFilesFromSubdirs(directories,
                           options.build_dir + "/predictions_aaa",
                           file_list,
                           options,
                           extract_id_from_dirname=None,
                           extract_id_from_filename=re.compile(
                               "export/export_aaa.dir/(\S+)_caf1"),
                           section="gene predictions for AAA",
                           within_subdirs=None,
                           readme=readme)

    return


def ExportPredictions(options):
    """export gene predictions for AAA."""

    directories = glob.glob(options.dir_export_predictions, )

    file_list = (("export_clustering_%s.exons", "%s.exons", None, ),
                 ("export_clustering_%s_cds.fasta", "%s_cds.fasta", None, ),
                 ("export_clustering_%s_peptides.fasta",
                  "%s_peptides.fasta", None, ),
                 )

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

    ExportFilesFromSubdirs(directories,
                           options.build_dir + "/predictions",
                           file_list,
                           options,
                           extract_id_from_dirname=None,
                           extract_id_from_filename=re.compile(
                               "export/export_clustering.dir/(\S+_vs_dmel\d+)"),
                           section="gene predictions",
                           within_subdirs=None,
                           readme=readme)
    return


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gpipe/export_all.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-r", "--release", dest="release", type="string",
                      help="release.")

    parser.add_option("-s", "--sections", dest="sections", type="string",
                      help="sections to export: malis,codonbias,aaa,predictions,orthologs")

    parser.add_option("-p", "--package-name", dest="package_name", type="string",
                      help="Package name. This also provides the build directory.")

    parser.add_option("-c", "--compress", dest="compress", action="store_true",
                      help="Compress files.")

    parser.set_defaults(
        sections="malis,codonbias,orthologs,aaa,predictions",
        dir_export_malis="../orthology_malis",
        dir_export_orthologs="../orthology_multiple",
        dir_export_codonbias="../codonbias",
        dir_export_aaa="../export/export_aaa.dir",
        dir_export_predictions="../export/export_clustering.dir",
        package_name="gpipe",
        build_dir=None,
        compress=False,
        release="unknown",
    )

    (options, args) = E.Start(parser, add_database_options=True)

    options.build_dir = options.package_name + "_" + options.release
    options.sections = options.sections.split(",")
    dbhandle = pgdb.connect(options.psql_connection)

    if not os.path.exists(options.build_dir):
        os.mkdir(options.build_dir)

    for section in options.sections:

        if section == "malis":
            ExportMalis(options)
        elif section == "codonbias":
            ExportCodonbias(options)
        elif section == "orthologs":
            ExportOrthologs(options)
        elif section == "aaa":
            ExportAAA(options)
        elif section == "predictions":
            ExportPredictions(options)

    outfile = open(options.build_dir + "/readme", "w")

    outfile.write(
        "Protein coding gene predictions and analysis of gene predictions of 12 fly genomes.\n")
    outfile.write(
        "\nAndreas Heger and Chris Ponting, MRC Functional Genetics Unit, Oxford, UK.\n")
    outfile.write("\ncreated on %s.\n\n" %
                  (time.asctime(time.localtime(time.time()))))

    readme = """This directory contains gene predictions and annotatios in the twelve fruit fly
genomes using a pipeline based exonerate (Slater et al. (2005)).  

Contents

orthologs:              clusters of orthologous groups

multiple_alignments:    multiple alignments of orthologous groups

codonbias:              codon bias analysis of predicted genes.

predictions_all:        predicted transcripts and genes.

predictions_aaa:        gene predictions as submitted to the consensus
                        annotation project.


Consensus annotations

The predictions use here constitute a super-set of the predictions submitted to the
consensus annotation project at http://rana.lbl.gov/drosophila/wiki/index.php/Main_Page.
For example, the set here includes pseudogenes. The mapping files in the directory
aaa_predictions map sequence identifiers used here to those submitted to the consensus
annotation project.

Transcript identifiers

The format for transcript identifiers is species|transcript_id|gene_id|quality:

species:        species name (TEXT). Species format is xxxx_vs_yyyyid, where
                xxxx is the species, yyyy is the template species (dmel in all cases)
                and id is the prediction run.

transcript_id:  transcript identifier (dmel: TEXT, others: INT)

gene_id:        gene identifier (dmel: TEXT, others: INT)

quality:        quality code. See http://wwwfgu.anat.ox.ac.uk:8080/flies/documentation.html
                for more information.

Note: I am aware that | is an unfortunate choice for a separator (I was looking for something readable,
that could be used in trees, belvu alignments, and others, which ruled out ':' and '/' and others).
I might clean up identifiers in the future.

Input data

Genomes are those in the comparative analysis freeze 1 (caf1). D. melanogaster annotations
are those from flybase obtained via ensembl.
"""

    outfile.write(readme + "\n")
    outfile.close()

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
