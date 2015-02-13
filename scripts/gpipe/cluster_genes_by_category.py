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
gpipe/cluster_genes_by_category.py - 
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

   python gpipe/cluster_genes_by_category.py --help

Type::

   python gpipe/cluster_genes_by_category.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import string
import getopt
import CGAT.Experiment as E
import CGAT.Genomics as Genomics

USAGE = """python %s [OPTIONS] < exonerate_output > filtered

Version: $Id: gpipe/cluster_genes_by_category.py 14 2005-08-09 15:24:07Z andreas $

Cluster genes belonging to the same category

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-g, --go=                       filename with go categories for genes
-d, --distance-method=                 distance for clustering
-r, --regions=                  filename with regions for genes
""" % sys.argv[0]


param_loglevel = 2

param_long_options = ["verbose=", "help",
                      "go=", "regions=",
                      "distance=",
                      "version"]
param_short_options = "v:hg:d:r:"

param_filename_go = None
param_filename_regions = None

param_distance = 1000000


# ------------------------------------------------------------

def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    try:
        optlist, args = getopt.getopt(
            sys.argv[1:], param_short_options, param_long_options)
    except getopt.error, msg:
        print USAGE, msg
        sys.exit(2)

    for o, a in optlist:
        if o in ("-v", "--verbose"):
            param_loglevel = int(a)
        elif o in ("--version", ):
            print "version="
            sys.exit(0)
        elif o in ("-h", "--help"):
            print USAGE
            sys.exit(0)
        elif o in ("-g", "--go"):
            param_filename_go = a
        elif o in ("-d", "--distance-method"):
            param_distance = int(a)
        elif o in ("-r", "--regions"):
            param_filename_regions = a

    print E.GetHeader()
    print E.GetParams()

    queries = {}
    for line in sys.stdin:
        if line[0] == "#":
            continue
        queries[line[:-1].split("\t")[0]] = 1

    if param_loglevel >= 1:
        print "# read %i queries" % len(queries)

    if param_filename_go:
        go = Genomics.ReadGo(open(param_filename_go, "r"))

    if param_filename_regions:
        locations = Genomics.ReadLocationsGFF(
            open(param_filename_regions, "r"))

    # compile a list of genes in the same go category. Save location
    map_go2genes = {}
    for gene_id in go:
        for assignment in go[gene_id]:
            if assignment.mGoid not in map_go2genes:
                map_go2genes[assignment.mGoid] = []
            if gene_id not in locations:
                continue
            map_go2genes[assignment.mGoid].append(locations[gene_id])

    for goid in map_go2genes:
        geneids = map_go2genes[goid]
        print "goid=", goid, "geneids=", len(geneids)
        geneids.sort(lambda x, y: cmp(
            (x.mSbjctToken, x.mSbjctStrand, x.mSbjctGenomeFrom),
            (y.mSbjctToken, y.mSbjctStrand, y.mSbjctGenomeFrom)))

        last = None
        for this in geneids:
            if last and \
               last.mSbjctToken == this.mSbjctToken and \
               last.mSbjctStrand == this.mSbjctStrand and \
               abs(min(last.mSbjctGenomeTo, this.mSbjctGenomeTo) -
                    max(last.mSbjctGenomeFrom, this.mSbjctGenomeFrom)) < param_distance:
                print "# removed duplicate:", str(last)
                print "# ", last.mQueryToken
                print string.join(map(str, go[last.mQueryToken]), "\n")
                print "# ", this.mQueryToken
                print string.join(map(str, go[this.mQueryToken]), "\n")
            else:
                print str(last)

            last = this

        if last:
            print str(last)

    print E.GetFooter()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
