'''
bam_vs_bed.py - count context that reads map to
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics NGS Intervals BAM BED Counting

Purpose
-------

This script takes as input a :term:`BAM` file from an RNASeq experiment 
and a :term:`bed` formatted file. The :term:`bed` formatted file needs
at least four columns. The fourth (name) column is used to group counts.

It counts the number of alignments overlapping between the :term:`bam`
file and the :term:`bed` file. Annotations in the :term:`bed` file can
be overlapping - they are counted independently.

This scripts requires bedtools to be installed.

Usage
-----

Example::

   python bam_vs_bed.py in.bam in.bed.gz

Type::

   python bam_vs_bed.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import time
import subprocess
import tempfile
import collections
import itertools

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import pysam
import CGAT.Bed as Bed

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option( "-m", "--min-overlap", dest="min_overlap", type="float",
                       help = "minimum overlap [%default]" )

    parser.add_option( "-k", "--keep-temp", dest="keep_temp", action="store_true",
                       help = "do not delete temporary files [%default]" )

    parser.add_option( "-a", "--filename-bam", dest="filename_bam", metavar="bam", type="string",
                       help = "bam-file to use [%default]" )

    parser.add_option( "-b", "--filename-bed", dest="filename_bed", metavar="bam", type="string",
                       help = "bed-file to use [%default]" )

    parser.set_defaults(
        min_overlap = 0.5,
        keep_temp = False,
        filename_bam = None,
        filename_bed = None,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    filename_bam = options.filename_bam
    filename_bed = options.filename_bed

    if filename_bam == None and filename_bed == None:
        if len(args) != 2:
            raise ValueError( "please supply a bam and a bed file or two bed-files." )
    
        filename_bam, filename_bed = args

    if filename_bed == None:
        raise ValueError( "please supply a bed file to compare to." )

    if filename_bam == None:
        raise ValueError( "please supply a bam file to compare with." )

    E.info( "intersecting the two files" )

    tmpfile = tempfile.NamedTemporaryFile( delete=False )
    tmpfile.close()
    tmpfilename = tmpfile.name

    min_overlap = options.min_overlap

    options.stdout.write( "category\talignments\n" )

    # get number of columns of reference bed file
    for bed in Bed.iterator(IOTools.openFile( filename_bed )):
        ncolumns_bed = bed.columns
        break
    E.info( "assuming %s is bed%i format" % (filename_bed, ncolumns_bed))

    if ncolumns_bed < 4:
        raise ValueError("please supply a name attribute in the bed file")

    # get information about 
    if filename_bam.endswith(".bam"):
        format = "-abam"
        samfile = pysam.Samfile( filename_bam, "rb" )
        total = samfile.mapped
        # latest bedtools uses bed12 format when bam is input
        ncolumns_bam = 12
        # count per read
        sort_key = lambda x: x.name
    else: 
        format = "-a"
        total = IOTools.getNumLines( filename_bam )
        # get bed format
        ncolumns_bam = 0
        for bed in Bed.iterator(IOTools.openFile( filename_bam )):
            ncolumns_bam = bed.columns
            break

        if ncolumns_bam > 0:
            E.info( "assuming %s is bed%i fomat" % (filename_bam, ncolumns_bam ))
            if ncolumns_bam == 3:
                # count per interval
                sort_key = lambda x: (x.contig, x.start, x.end)
            else:
                # count per interval category
                sort_key = lambda x: x.name

    # use fields for bam/bed file (regions to count with)
    data_fields = [ "contig", "start", "end", "name",
                    "score", "strand", "thickstart", "thickend", "rgb",
                    "blockcount", "blockstarts", "blockends" ][:ncolumns_bam]

    # add fields for second bed (regions to count in)
    data_fields.extend( [ "contig2", "start2", "end2", "name2",
                          "score2", "strand2", "thickstart2", "thickend2", "rgb2",
                          "blockcount2", "blockstarts2", "blockends2" ][:ncolumns_bed] )

    # add bases overlap
    data_fields.append( "bases_overlap"  )

    data = collections.namedtuple( "data", data_fields )

    options.stdout.write( "total\t%i\n" % total )

    if total == 0:
        E.warn( "no data in %s" % filename_bam )
        return

    statement = """intersectBed %(format)s %(filename_bam)s -b %(filename_bed)s -bed -wo -f %(min_overlap)f > %(tmpfilename)s""" % locals()

    E.info( "running %s" % statement )
    retcode = E.run( statement )

    if retcode != 0:
        raise ValueError( "error while executing statement %s" % statement )

    infile = open( tmpfilename, "r")
    counts_per_alignment = collections.defaultdict(int)

    E.info( "counting" )

    take_columns = len(data._fields)

    def iter( infile ):
        for line in infile:
            if not line.strip(): continue
            yield data._make( line[:-1].split()[:take_columns] )

    for read, overlaps in itertools.groupby( iter(infile), key = sort_key ):
        annotations = [x.name2 for x in overlaps ]
        for anno in annotations:
            counts_per_alignment[anno] += 1
    infile.close()

    for key, counts  in counts_per_alignment.iteritems():
        options.stdout.write( "%s\t%i\n" % (key, counts) )

    if not options.keep_temp:
        os.unlink( tmpfilename )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

