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
bam2transcriptContribution.py
=============================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Computes the total number of reads from a BAM alignment file that contribute to
transcripts provided in a gtf file.



Usage
-----

Example::

   python bam2transcriptContribution.py --help

Type::

   python bam2transcriptContribution.py --help

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
import collections
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import pysam
import CGAT.GTF as GTF
import numpy
import pybedtools
import tempfile
import pyximport
pyximport.install(build_in_temp=False)
import CGAT.Pipeline as P

#########################################################
#########################################################
#########################################################

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-b", "--bam-file", dest="bam_file", type="string",
                      help="supply input bam file name"  )
   
    parser.add_option("-g", "--gtf-file", dest="gtf_file", type="string",
                      help="supply input gtf file name"  )
   
    parser.add_option("-o", "--outfile", dest = "outfile", type = "string",
                      help="supply output file name")

    parser.add_option("-G", "--reference-GTF", dest = "reference_gtf", type = "string",
                      help="supply reference gtf for context of reads not contributing to transcripts")


    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    ######################################################
    ######################################################
    # for all alignments
    ######################################################
    ######################################################

    # open outfile and prepare headers
    outf = open(options.outfile, "w")
    outf.write("\t".join(["total alignments"
                          , "aligments in transcripts"
                          , "percent alignments in transcripts"
                          , "total spliced alignments"
                          , "spliced alignments in transcripts"
                          , "percent spliced alignments in transcripts"]) + "\n")

    # calculate coverage over transcript file - NB split reads contribute twice to the transcript
    # use BedTool object
    pybedbamfile = pybedtools.BedTool(options.bam_file)
    
    # count alignments 
    E.info("counting total number of alignments and spliced alignments")
    total_alignments = 0
    spliced_alignments = 0
    
    for alignment in pybedbamfile:
        cigar = alignment[5]
        if cigar.find("N") != -1: # N signifies split read 
            total_alignments += 1
            spliced_alignments += 1
        else:
            total_alignments += 1

    # merge the gtf file to avoid double counting of exons in different transcripts - converts to a bed file
    gtffile = pybedtools.BedTool(options.gtf_file).merge()

    E.info("computing coverage of aligments in %s over intervals in %s" % (options.bam_file, options.gtf_file))
    cover = pybedbamfile.coverage(gtffile)

    # make sure that the exons aren't being counted twice - shouldn't be because of merge
    E.info("counting reads contributing to transcripts")
    c = 0
    for entry in cover:
        coverage = int(entry[3])
        if coverage > 0:
            c+=coverage
        
    # sum the coverage across exons from all transcripts
    coverage_in_transcripts = c

    ######################################################
    ######################################################
    # for spliced alignments
    ######################################################
    ######################################################

    # count total number of spliced alignments
    # requires that the CIGAR string 'N' is present
    
    # uses pysam to write out a bam file of the spliced reads only
    allreads = pysam.Samfile(options.bam_file)
    spliced_bamname = P.snip(options.bam_file, ".bam") + "_spliced_reads.bam"
    
    # open file for outputting spliced alignments
    splicedreads = pysam.Samfile(spliced_bamname, "wb", template=allreads)
    
    # cigar string in pysam for spliced alignment is (3, int)
    spliced = collections.defaultdict(list)
    for read in allreads:
        for cigar_tag in read.cigar:
            if cigar_tag[0] == 3:
                spliced[read].append(cigar_tag)

    # write out spliced alignments
    for read in spliced.keys():
        splicedreads.write(read)
    splicedreads.close()
    allreads.close()
    
    # index splice reads bam file
    pysam.sort(spliced_bamname, P.snip(spliced_bamname, ".bam"))
    pysam.index(spliced_bamname)

    # read in the spliced reads as a BedTool object
    splicedbam = pybedtools.BedTool(spliced_bamname)

    # perform coverage of spliced reads over intervals - will be twice as many as there should be
    # due to counting both exons overlapping
    spliced_coverage = splicedbam.coverage(gtffile)

    # avoid double counting exons
    E.info("counting spliced reads contributing to transcripts")
    spliced_exons = {}
    c = 0
    for entry in spliced_coverage:
        coverage = int(entry[3])
        if coverage > 0:
            c+=coverage
                
    spliced_coverage_in_transcripts = c

    # NOTE: the counting of spliced alignments is not accurate

    spliced_coverage_in_transcripts =  float(spliced_coverage_in_transcripts)/2
    
    ###########################
    ## write out the results ##
    ###########################

    outf.write(str(int(total_alignments)) + "\t")
    # remove half of the coverage assigned to spliced reads
    coverage_in_transcripts = (coverage_in_transcripts)-(spliced_coverage_in_transcripts)
    outf.write(str(int(coverage_in_transcripts)-int(spliced_coverage_in_transcripts)) + "\t")
    outf.write( str( int((coverage_in_transcripts/total_alignments) *100)) + "\t" )

    # write out spliced counts
    outf.write(str(int(spliced_alignments)) + "\t")
    outf.write(str(int(spliced_coverage_in_transcripts)) + "\t")
    outf.write( str( int( (spliced_coverage_in_transcripts/spliced_alignments) *100)) )

    outf.close()

    ############################
    # contextualise those that 
    # don't fall in transcripts
    ############################

    if options.reference_gtf:
        context_summary = open(P.snip(options.bam_file, ".bam") + ".excluded.context", "w")
        context_summary.write("\t".join(["Feature", "number"]) + "\n")

        # write out the read info as well
        context_file = open(P.snip(options.bam_file, ".bam") + ".excluded", "w")
    
        context_dict = collections.defaultdict(int)
        # intersect bam - write non-overlapping with transcripts - intersect with reference - write out
        context = pybedbamfile.intersect(gtffile, v=True, bed=True).intersect(pybedtools.BedTool(options.reference_gtf), wb=True)
        for entry in context:
            feature = entry[8]
            context_dict[feature] += 1
            context_file.write("\t".join([e for e in entry]) + "\n")

        for feature, value in context_dict.iteritems():
            context_summary.write("\t".join([feature, str(value)]) + "\n")

        context_file.close()
        context_summary.close()


    ## write footer and output benchmark information.
    E.Stop()


if __name__ == "__main__":
    sys.exit( main( sys.argv) )





    




    

