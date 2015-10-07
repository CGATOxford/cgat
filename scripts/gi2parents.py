##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
'''gi2parents.py
=============================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script takes a list of gi accession numbers - have to be
sub-species specific - from the NCBI and outputs a list of parent taxa
- phylum, class, order, family, genus, species.

Requires ncbi.map file downloaded from the MEGAN website:
http://www-ab2.informatik.uni-tuebingen.de/megan/taxonomy/ncbi.zip

Usage
-----

Example::

   python gi2parents.py -g gi_accessions.tsv -m ncbi.map -n gi_taxid_nucl.dmp.gz -c ncbi.lvl -t nodes.dmp

Type::

   python gi2parents.py --help

for command line help.

Code
----

'''

import sys
import optparse
import CGAT.IOTools as IOTools
import CGAT.Experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser(version="%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $",
                                   usage=globals()["__doc__"])

    parser.add_option("-g", "--gi-accessions", dest="gi_accessions", type="string",
                      help="list of gi accession numbers")
    parser.add_option("-m", "--ncbi-map", dest="ncbi_map", type="string",
                      help="ncbi.map file downloaded from the MEGAN website")
    parser.add_option("-n", "--nucl-map", dest="nucl_map", type="string",
                      help="gi mapping to tax id downloaded from ncbi website")
    parser.add_option("-c", "--taxa-code", dest="taxa_code", type="string",
                      help="code for different levels of the taxonomy downloaded from the MEGAN website")
    parser.add_option("-t", "--tree", dest="tree", type="string",
                      help="description of parents in the taxonomy")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    E.info("reading gi accession numbers")
    gi_accessions = set()
    for line in open(options.gi_accessions).readlines():
        gi_accessions.add(line[:-1])
    E.info("read gi accession numbers")

    E.info("building gi2taxid map")
    gi2taxid = {}
    c_gi = 0
    for line in IOTools.openFile(options.nucl_map).readlines():
        data = line[:-1].split("\t")
        if data[0] not in gi_accessions:
            continue
        else:
            c_gi += 1
            gi2taxid[data[0]] = data[1]
    E.info("built gi2taxid for %i gi accession numbers" % c_gi)

    E.info("building code map")
    code2taxa = {}

    for line in open(options.taxa_code).readlines():
        data = line[:-1].split("\t")
        code2taxa[data[0]] = data[1]
    E.info("built taxa code map")

    E.info("building taxa2name map")
    taxid2name = {}
    for line in open(options.ncbi_map).readlines():
        data = line[:-1].split("\t")
        # keep the taxa code
        taxid2name[data[0]] = (data[1], data[3])
    E.info("built taxa2name map")

    E.info("build taxid2parentmap")
    taxid2parents = {}
    for line in open(options.tree).readlines():
        data = line[:-1].split("\t")
        data = [x for x in data if x != "|"]
        taxid2parents[data[0]] = data[1]
    E.info("built taxid2parentmap")

    E.info("retrieving parents for each gi accession number")
    options.stdout.write(
        "gi\tsub_species\tspecies\tgenus\tfamily\torder\tclass\tphylum\n")
    for gi, taxid in gi2taxid.iteritems():
        # this will be the sub species id
        # walk through the parents
        parents = {}
        sub_species = taxid2name[taxid][0]
        for i in range(len(code2taxa.keys())):
            parent_taxid = taxid2parents[taxid]
            parent_name = taxid2name[parent_taxid][0]
            parent_code = taxid2name[parent_taxid][1]
            # ignore codes that we are not  interested in
            if parent_code not in code2taxa.keys():
                continue
            parent_taxa = code2taxa[parent_code]
            parents[parent_taxa] = parent_name
            taxid = parent_taxid

        if "genus" not in parents:
            genus = "NA"
        else:
            genus = parents["genus"]
        if "family" not in parents:
            family = "NA"
        else:
            family = parents["family"]
        if "order" not in parents:
            order = "NA"
        else:
            order = parents["order"]
        if "class" not in parents:
            _class = "NA"
        else:
            _class = parents["class"]
        if "phylum" not in parents:
            phylum = "NA"
        else:
            phylum = parents["phylum"]
            if phylum.find("<phylum>") != -1:
                phylum = phylum.replace(" <phylum>", "")
        if "species" not in parents:
            species = "NA"
        else:
            species = parents["species"]
        options.stdout.write("\t".join([gi, sub_species.replace(" ", "_"), species.replace(
            " ", "_"), genus, family, order, _class, phylum]) + "\n")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
