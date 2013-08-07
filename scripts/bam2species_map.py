################################################################################
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
#################################################################################
'''
bam2species_map.py
=============================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Produce a mapping txt file between contigs a species based
on aligned reads.

Usage
-----

Example::

   python bam2species_map.py --help

Type::

   python bam2species_map.py --help

for command line help.

Documentation
-------------

This script would be used as a pre-step to using contigs2random_sample.py. It provides
a mapping between contigs and species that are represented in those contigs i.e. in 
a metagenomic simulation study the majority species for a contig will be returned with
the contig.

Code
----

'''

import os
import sys
import re
import optparse

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import pysam
import collections
import CGAT.FastaIterator as FastaIterator

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-b", "--bamfile", dest="bamfile", type="string",
                      help="supply bam file"  )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    # read in contigs
    E.info("reading in contig file")
    contigs = {}
    for fasta in FastaIterator.iterate(options.stdin):
        contigs[fasta.title] = (1, len(fasta.sequence) -1)
    E.info("read %i contigs" % len(contigs.keys()))

    # read in bamfile
    E.info("reading bam file")
    samfile = pysam.Samfile(options.bamfile)

    E.info("iterating over contigs")
    c = 0
    for contig, coords in contigs.iteritems():
        coords = list(coords)

        #################################
        # NB this is specific for my data!
        contig = contig.split(" ")[0]
        #################################

        species_counts = collections.defaultdict(int)
        for alignment in samfile.fetch(contig, coords[0], coords[1]):
            species_id = alignment.qname.split("|")[1]
            species_counts[species_id] += 1

        # at the moment ignore if there are no counts
        if len(species_counts.values()) == 0: 
            E.warn("no reads map to %s" % contig)
            continue

        for species, count in species_counts.iteritems():
            if species_counts[species] == max(species_counts.values()): 
                top_dog = species
                c += 1
                break
        E.info("species %s assigned to contig number %i" % (top_dog, c))
        options.stdout.write("%s\t%s\n" % (contig, top_dog))

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
