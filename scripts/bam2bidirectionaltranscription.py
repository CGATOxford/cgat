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
bam2bidirectionaltranscription.py
=============================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

The purpose of this script is two-fold.

1. Stranded RNA-seq libraries are becoming more and more commonplace. As such, a prudent 
   quality control step after alignment to assess the observed strandedness across a set of transcripts of
   interest (GTF file as input).

2. Bidirectional transcription at intergenic loci is an area of interest across biological science.
   This script allows the user to assess the level of bidirectional transcription across a set
   of intervals of interest (GTF file as input).


Usage
-----

Example::

   python bam2bidirectionaltranscription.py --help

Type::

   python bam2bidirectionaltranscription.py --help

for command line help.

Documentation
-------------

Requires both the NH and XS flags to be set in the bam file.
This is specific to the aligner used. TopHat will add an XS and NH
flag for example.

Doesn't use paired information.
 

Code
----

'''

import os
import sys
import re
import optparse
import pysam
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.Experiment as E
import CGAT.IndexedGenome as IndexedGenome
import numpy as np
import collections


def filterGtf(gtf_iterable, bamfile, tag_count = 10, length = 200, unique = True):
    '''
    filter the gtf file based on filtering options
    - This does not account for strand at this point
    '''
    filtered_set = set()
    total = 0
    filtered = 0
    for gtf in gtf_iterable:
        total += 1
        tag_number = 0
        if (gtf.end + 1)  - gtf.start < length:
            filtered += 1
            continue
        start, end = gtf.start-1, gtf.end
        for alignment in bamfile.fetch(gtf.contig, start, end):
            # get the number of hits
            nh = [x for x in alignment.tags if x[0] == "NH"]
            nh = nh[0][1]
            if unique:
                if nh != 1: continue
            tag_number += 1
        if tag_number < tag_count: 
            filtered += 1
            continue
        filtered_set.add(gtf.gene_id)
    E.info("total gtf entries in = %i...total filtered = %i" % (total, filtered))
    return filtered_set


def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )
    parser.add_option("-b", "--bam", dest="bam", type="string",
                      help="input bam file"  )
    parser.add_option("-a", "--annotation", dest="annotation", type="string",
                      help="input annotation gtf file"  )
    parser.add_option("-m", "--merge-transcripts", dest="merge_transcripts", action="store_true"
                      , help="merge transcripts with the same gene id"  )
    parser.add_option("-u", "--unique", dest="unique", action="store_true"
                      , help="only use unique alignments"  )
    parser.add_option("-t", "--tag-count", dest="tag_count", type="int"
                      , help="minimum tag count to allow inclusion in analysis"  )
    parser.add_option("-l", "--length", dest="tag_count", type="int"
                      , help="minimum length"  )
    parser.add_option("-n", "--analysis-type", dest = "analysis_type", type = "choice", choices = ("ratio", "profile", "shape"))
    parser.add_option("-o", "--outbase", dest="outbase", type = "string"
                      , help="basename for outfiles"  )
    parser.add_option("-s", "--bin-number", dest="bin_number", type = "int"
                      , help="number of bins to use - only active with --profile"  )

    parser.set_defaults(unique = True
                        , merge_transcripts = True
                        , tag_count = 10
                        , length = 200
                        , analysis_type = "ratio"
                        , bin_number = 100
                        , outbase = "bidirection"
                        , output_gtf = False)
                        

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    bam = pysam.Samfile(options.bam)
    gtffile = IOTools.openFile(options.annotation)

    if options.merge_transcripts:
        iterable = GTF.merged_gene_iterator(GTF.iterator(gtffile))
    else:
        iterable = GTF.iterator(gtffile)

    E.info("filtering gtf set")
    gtf_set = filterGtf(iterable, bam, tag_count = options.tag_count, length = options.length, unique = options.unique)

    # reinitiate iterable
    gtffile = IOTools.openFile(options.annotation)
    if options.merge_transcripts:
        iterable = GTF.merged_gene_iterator(GTF.iterator(gtffile))
    else:
        iterable = GTF.iterator(gtffile)
    
    ################################################
    # output the overall ratio of sense to antisense
    # reads for the gtf entries filtered in gtf_set
    ################################################
    if options.analysis_type == "ratio":
        E.info("calculating ratio of sense / antisense")
        outf = open(options.outbase + "_ratio.tsv", "w")
        outf.write("gene_id\ttranscript_id\tasratio\n")
        total_unmapped = 0
        total_no_xs = 0
        total_reads_in = 0
        for gtf in iterable:
            if gtf.gene_id not in gtf_set: continue
            transcription = {"sense": 0, "antisense": 0}
            start, end = gtf.start-1, gtf.end
            for alignment in bam.fetch(gtf.contig, start, end):
                total_reads_in += 1
                if alignment.is_unmapped:
                    total_unmapped += 1
                    continue
                # get the strand origin
                xs = [x for x in alignment.tags if x[0] == "XS"]
                if len(xs) == 0:
                    E.warn("No XS tag present, read %s not included in the analysis" % alignment.qname)
                    total_no_xs += 1
                    continue
                xs = xs[0][1]
                if xs == "+":
                    transcription["sense"] += 1
                elif xs == "-":
                    transcription["antisense"] += 1
                else:
                    raise ValueError("no such strand %s" % xs)
            if transcription["sense"] == 0 and transcription["antisense"] == 0:
                continue
            elif transcription["sense"] == 0 and transcription["antisense"] != 0:
                ratio = "-Inf"
            elif transcription["sense"] != 0 and transcription["antisense"] == 0:
                ratio = "Inf"
            else:
                gtf_set.add(gtf.gene_id)
                ratio = float(transcription["sense"]) / transcription["antisense"]
            outf.write("%s\t%s\t%s\n" % (gtf.gene_id, gtf.transcript_id, str(ratio)))
        outf.close()

    ####################################################
    # profiling average number of reads across intervals
    ####################################################
    elif options.analysis_type == "profile":
        profile = {}
        E.info("assessing profile over intervals in %i bins" % options.bin_number)
        total_unmapped = 0
        total_no_xs = 0
        total_reads_in = 0
        for gtf in iterable:
            if gtf.gene_id not in gtf_set: continue
            start, end = gtf.start-1, gtf.end
            bins = np.histogram(range(start, end), bins = options.bin_number)[1]
            bin_as = []
            for x in range(len(bins)):
                if x < len(bins)-1:
                    bin_size = int(bins[x+1]) - int(bins[x])
                    c_sense = 0
                    c_antisense = 0
                    for alignment in bam.fetch(gtf.contig, int(bins[x]), int(bins[x+1])):
                        total_reads_in += 1
                        if alignment.is_unmapped:
                            total_unmapped += 1
                            continue
                        # number of hits
                        nh = [x for x in alignment.tags if x[0] == "NH"]
                        nh = nh[0][1]
                        if options.unique:
                            if nh != 1: continue
                        # get the strand origin
                        xs = [x for x in alignment.tags if x[0] == "XS"]
                        if len(xs) == 0: 
                            E.warn("No XS tag present, read %s not included in the analysis" % alignment.qname)
                            total_no_xs += 1
                            continue
                        xs = xs[0][1]
                        if xs == "+":
                            c_sense += 1
                        elif xs == "-":
                            c_antisense += 1
                        else:
                            raise ValueError("no such strand %s" % xs)
                    if c_sense == 0 and c_antisense == 0:
                        # if bin is empty then add "NA"
                        ratio = "NA"
                    elif c_sense == 0 and c_antisense != 0:
                        ratio = "-Inf"
                    elif c_sense != 0 and c_antisense == 0:
                        ratio = "Inf"
                    else:
                        ratio = float(c_sense)/c_antisense
                    bin_as.append(ratio)
            profile[gtf.gene_id] = bin_as

        outf = open(options.outbase + "_profile.tsv", "w")
        outf.write("\t".join(profile.keys()) + "\n")
        for x in zip(*profile.values()):
            outf.write("\t".join(map(str,list(x))) + "\n")
        outf.close()

    ####################################################
    # provide table with sense and antisense
    # counts
    ####################################################
    elif options.analysis_type == "shape":
        sense = {}
        antisense = {}
        E.info("assessing shape over intervals in %i bins" % options.bin_number)
        total_unmapped = 0
        total_no_xs = 0 
        total_reads_in = 0
        for gtf in iterable:
            if gtf.gene_id not in gtf_set: continue
            start, end = gtf.start-1, gtf.end
            bins = np.histogram(range(start, end), bins = options.bin_number)[1]
            bin_as = []
            bin_s = []
            for x in range(len(bins)):
                if x < len(bins)-1:
                    bin_size = int(bins[x+1]) - int(bins[x])
                    c_sense = 0
                    c_antisense = 0
                    for alignment in bam.fetch(gtf.contig, int(bins[x]), int(bins[x+1])):
                        total_reads_in += 1
                        if alignment.is_unmapped:
                            total_unmapped += 1
                            continue
                        # number of hist
                        nh = [x for x in alignment.tags if x[0] == "NH"]
                        nh = nh[0][1]
                        if options.unique:
                            if nh != 1: continue
                        # get the strand origin
                        xs = [x for x in alignment.tags if x[0] == "XS"]
                        if len(xs) == 0: 
                            E.warn("No XS tag present, read %s not included in the analysis" % alignment.qname)
                            total_no_xs += 1
                            continue
                        xs = xs[0][1]
                        if xs == "+":
                            c_sense += 1
                        elif xs == "-":
                            c_antisense += 1
                        else:
                            raise ValueError("no such strand %s" % xs)
                    bin_s.append(c_sense)
                    bin_as.append(c_antisense)
            sense[gtf.gene_id] = bin_s
            antisense[gtf.gene_id] = bin_as

        outf = open(options.outbase + "_shape.tsv", "w")
        outf.write("\t".join(sense.keys()) + "\tstatus" + "\n")
        for x in zip(*sense.values()):
            outf.write("\t".join(map(str,list(x))) + "\tsense" + "\n")
        for x in zip(*antisense.values()):
            outf.write("\t".join(map(str,list(x))) + "\tantisense" + "\n")
        outf.close()

    E.info("number of reads in: %i" % total_reads_in)
    E.info("number of unmapped reads in bam file: %i" % total_unmapped)
    E.info("number of reads with no XS tag %i" % total_no_xs)

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
