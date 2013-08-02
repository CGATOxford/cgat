################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
gtf2Overlap.py
=============================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Take two gtf files, merge transcripts, intersect them and draw a venn diagram of the overlap.
This is a crude way if assessing the overlap between two gtf files. Results in a locus
intersection. The intersection keeps all records in gtf-a that intersect gtf-b

Usage
-----

Example::

   python gtf2Overlap.py --help

Type::

   python gtf2Overlap.py --help

for command line help.

Documentation
-------------

Code
----

'''

import os
import sys
import re
import optparse
import itertools
import random
import CGAT.Pipeline as P
import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
from rpy2.robjects import r as R
from rpy2.robjects.vectors import FloatVector
from rpy2.robjects.vectors import StrVector

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-a", "--gtf-a", dest="gtf_a", type="string",
                      help="supply a gtf file - will compress uncompressed files"  )
    parser.add_option("-b", "--gtf-b", dest = "gtf_b", type = "string",
                      help="supply a second gtf file - will compress uncompressed files")
    parser.add_option("-s", "--scripts-dir", dest = "scripts_dir", type = "string",
                      help="supply a location for accessory scripts")
    parser.add_option( "--no-venn", dest = "no_venn", action="store_true", 
                      help="set if no venn is to be drawn")

    
    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    gtf_files = [options.gtf_a, options.gtf_b]

    merged_files = []
    prefices = []
    E.info("merging gtf files")
    for gtf in gtf_files:
        if gtf.endswith(".gtf.gz"):
            outfile = P.snip(gtf, ".gtf.gz") + ".merged.gtf.gz"
            prefices.append(P.snip(gtf, ".gtf.gz"))
            merged_files.append(outfile)
            statement = '''zcat %s | python %s/gtf2gtf.py --merge-transcripts --log=%s.log | gzip > %s''' % (gtf, options.scripts_dir, outfile, outfile)
            P.run()
        elif gtf.endswith(".gtf"):
            outfile = P.snip(gtf, ".gtf") + ".merged.gtf.gz"
            prefices.append(P.snip(gtf,".gtf"))
            merged_files.append(outfile)
            statement = '''cat %s | python %s/gtf2gtf.py --merge-transcripts --log=%s.log | gzip  > %s''' % (gtf, options.scripts_dir, outfile, outfile)
            P.run()
        else:
            raise ValueError("cannot perform merge on %s: is not a gtf file" % gtf)

    for prefix in prefices:
        if options.gtf_a.find(prefix) != -1:
            gtf_a = prefix + ".merged.gtf.gz"
            prefix_a = prefix
        elif options.gtf_b.find(prefix) != -1:
            gtf_b = prefix + ".merged.gtf.gz"
            prefix_b = prefix

    E.info("intersecting gtf files")
    # intersect the resulting merged files

    scriptsdir = options.scripts_dir
    intersection_out = "_vs_".join([prefix_a, prefix_b]) + ".intersection.gtf.gz" 
    statement = '''intersectBed -a %(gtf_a)s -b %(gtf_b)s -s -wa
                 | python %(scriptsdir)s/gtf2gtf.py --merge-transcripts --log=log | gzip > %(intersection_out)s'''
    P.run()

    if not options.no_venn:
        E.info("producing venn diagram for %s vs %s..." % (options.gtf_a, options.gtf_b))
        # produce the venn diagram
        intersection_file = intersection_out
        gtf_a_merged = gtf_a
        gtf_b_merged = gtf_b

        # create dictionary key
        gtf_pair = (gtf_a_merged, gtf_b_merged)

        # containers for counts
        count_gtf_merged_a = 0
        count_gtf_merged_b = 0
        count_intersection = 0

        # create GTF iterator objects
        gtf_iterator_a = GTF.iterator(IOTools.openFile(gtf_pair[0]))
        gtf_iterator_b = GTF.iterator(IOTools.openFile(gtf_pair[1]))
        gtf_iterator_intersection = GTF.iterator(IOTools.openFile(intersection_file))

        # do the counts for each file
        E.info("counting entries in %s" % gtf_a)
        for entry in gtf_iterator_a:
            count_gtf_merged_a += 1
        print "counts for gtf-a: ",count_gtf_merged_a

        E.info("counting entries in %s" % gtf_b)
        for entry in gtf_iterator_b:
            count_gtf_merged_b += 1
        print "counts for gtf-b: ",count_gtf_merged_b

        E.info("counting entries in %s" % intersection_file)
        for entry in gtf_iterator_intersection:
            count_intersection += 1
        print "counts for intersection: ", count_intersection

        # this is the important bit - basically take an arbitrary list of numbers to represent the list of lincrna in the refnoncoding set
        # then use the intersection count to represent the overlapping section in the lincrna set and add a set of random numbers to this 
        # set to make up the remaining - non-overlapping set

        result = {}
        E.info("assembling count lists")
        result[gtf_pair] = {"gtf-b" : map(str,xrange(count_gtf_merged_b))  , "gtf-a" : map(str,xrange(count_intersection)) + map(str, [random.random() for i in range(count_intersection,count_gtf_merged_a)]  )}

        R_source = os.path.join(os.path.abspath(options.scripts_dir), "venn_diagram.R")
        R.source(R_source)

        prefix_a = prefix_a.replace(".", "_").replace("-", "_")
        prefix_b = prefix_b.replace(".", "_").replace("-", "_")
        
        R('''prefix.a <- "%s"''' % prefix_a)
        R('''prefix.b <- "%s"''' % prefix_b) 
        E.info("drawing venn diagram to %s" % (prefix_a + "_vs_" + prefix_b + ".overlap.png"))
        
        R["venn.diagram2"](R.list( A = result[gtf_pair]["gtf-a"], B = result[gtf_pair]["gtf-b"])
        , prefix_a + "_vs_" + prefix_b + ".overlap.png"
        , **{'cat.cex': 1.5
             , 'main.fontfamily': "Arial"
             , 'cat.pos':FloatVector((0,0))
             , 'cat.fontfamily':"Arial"
             , 'main.cex':1.8                                                                                                                                                                                                              
             , 'height':1000
             , 'width':1000
             , 'cex':2                                                                                                                                                                                                                      
             , 'fontfamily':"Arial"                                                                                                                                                                                                         
             , 'lwd':R.c(1,1)                                                                                                                                                                                                               
             , 'fill':R.c(R.rgb(0,0,0.5,0.5), R.rgb(0.5,0,0,0.5))                                                                                                                                                         
             , 'category.names':R.c(prefix_a, prefix_b) 
             , 'margin' : R.c(0.1,0.1,0.1,0.1)
             })

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

