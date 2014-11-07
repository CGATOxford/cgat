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
optic_optic/analyze_cluster_expansion.py - analyse values based on a clustering.
==========================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

given:
        a map of predictions to components

Usage
-----

Example::

   python optic_optic/analyze_cluster_expansion.py --help

Type::

   python optic_optic/analyze_cluster_expansion.py --help

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

import CGAT.Experiment as E

import numpy
import scipy
import scipy.stats

parser = E.OptionParser(
    version="%prog version: $Id: optic/analyze_cluster_expansion.py 2781 2009-09-10 11:33:14Z andreas $")


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser.add_option("-m", "--map", dest="map", type="string",
                      help="filename of map with ids to components.")
    parser.add_option("-e", "--methods", dest="methods", type="string",
                      help="methods to apply to dataset.")
    parser.add_option("-s", "--pattern-species", dest="species_pattern", type="string",
                      help="regular expression to extract species.")
    parser.add_option("-p", "--pattern-identifier", dest="pattern", type="string",
                      help="output pattern. Should contain two %s wildcards, one for the secion, the second for the method.")

    parser.set_defaults(
        map=None,
        species_pattern="^([^@:|]+)[@:|]",
        pattern="analysis_%s.data",
        methods="hists",
        separator="|",
        status='CG,SG,PG,RG,UG,CP,SP,PP,RP,UP,SF,CF,PF,UF,BF',
        species='pdmel_vs_dmel,pdsim_vs_dmel2,pdyak_vs_dmel4,pdere_vs_dmel3,pdana_vs_dmel3,pdpse_vs_dmel3,pdvir_vs_dmel3,pdmoj_vs_dmel3,pdgrim_vs_dmel2',
    )

    (options, args) = E.Start(parser)

    options.species = options.species.split(",")
    options.status = options.status.split(",") + ["all"]
    options.methods = options.methods.split(",")

    if options.map:
        map_id2component = {}
        map_component2id = {}
        lines = filter(
            lambda x: x[0] != "#", open(options.map, "r").readlines())
        for i, c in map(lambda x: x[:-1].split("\t"), lines):
            map_id2component[i] = c
            if c not in map_component2id:
                map_component2id[c] = []
            map_component2id[c].append(i)

    rx_species = re.compile(options.species_pattern)

    outfile_species = open(options.pattern % "species", "w")
    outfile_cluster = open(options.pattern % "cluster", "w")

    outfile_species.write("clus\tspecies")
    total_status_over_all = {}
    for s in options.status:
        total_status_over_all[s] = 0
        outfile_species.write("\t%s" % s)
    outfile_species.write("\n")

    for c, members in map_component2id.items():

        counts = {}
        for a in options.species:
            counts[a] = {}
            for s in options.status:
                counts[a][s] = []

        for m in members:
            species, prediction_id, gene_id, status = m.split(
                options.separator)
            counts[species][status].append((gene_id, prediction_id))

        total_status_over_species = {}
        for s in options.status:
            total_status_over_species[s] = 0

        tt = {}
        nspecies = 0
        for a in options.species:
            outfile_species.write("%s\t%s" % (c, a))
            t = {}
            for s in options.status:
                g = {}
                for x in counts[a][s]:
                    g[x[0]] = 1
                    t["%s_%s" % (a, x[0])] = 1
                    tt["%s_%s" % (a, x[0])] = 1

                outfile_species.write("\t%i" % len(g))
                total_status_over_species[s] += len(g)

            outfile_species.write("\t%i" % len(t))
            outfile_species.write("\n")
            if len(t) > 0:
                nspecies += 1

        outfile_species.write("%s\tall" % (c))
        for s in options.status:
            outfile_species.write("\t%i" % total_status_over_species[s])
            total_status_over_all[s] += total_status_over_species[s]

        outfile_species.write("\t%i" % len(tt))
        outfile_species.write("\n")

        print "# cluster=", c, len(tt), len(members), nspecies

    print outfile_species.write("all\tall")
    for s in options.status:
        outfile_species.write("\t%i" % total_status_over_all[s])

    outfile_species.close()
    outfile_cluster.close()

    print "# ncomponents=%i" % (len(map_component2id))
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
