################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: bam2wiggle.py 2832 2009-11-24 16:11:06Z andreas $
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
"""
bam2wiggle.py - convert bam to wig/bigwig file
==============================================

:Author: Andreas Heger
:Release: $Id: bam2wiggle.py 2832 2009-11-24 16:11:06Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

convert a bam file to a bigwig or bedgraph file.   

The script requires the executables :file:`wigToBigWig` and 
:file:`bedToBigBed` to be in the user's PATH.

Usage
-----

Type::

   python bam2wiggle.py --output-format=bigwig in.bam out.bigwig

to convert the :term:`bam` file file:`in.bam` to :term:`bigwig` format 
and save the result in :term:`out.bigwig`.

Type::

   python bam2wiggle.py --help

for command line help.

Code
----

""" 

import os, sys, re, optparse, itertools, tempfile, shutil, subprocess

import Experiment as E
import pysam
import IOTools

def main( argv = None ):
    """script main.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: bam2wiggle.py 2832 2009-11-24 16:11:06Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-o", "--output-format", dest="output_format", type="choice",
                      choices=("bedgraph", "wiggle", "bigbed", "bigwig", "bed"),
                      help="output format [default=%default]" )

    parser.add_option("-b", "--output-filename", dest="output_filename", type="string",
                      help="filename for output [default=%default]" )

    parser.add_option( "-s", "--shift", dest="shift", type = "int",
                       help = "shift reads by a certain amount (ChIP-Seq) [%default]" )

    parser.add_option( "-e", "--extend", dest="extend", type = "int",
                       help = "extend reads by a certain amount (ChIP-Seq) [%default]" )

    parser.set_defaults(
        samfile = None,
        output_format = "wiggle",
        output_filename = None,
        shift = 0,
        extend = 0,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if len(args) >= 1:
        options.samfile = args[0]

    if not options.samfile:
        raise ValueError("please provide a bam file")

    if len(args) == 2:
        options.output_filename = args[1]

    samfile = pysam.Samfile( options.samfile, "rb" )

    contig_sizes = dict( zip( samfile.references, samfile.lengths) )

    if options.shift or options.extend:
        if options.output_format != "bigwig":
            raise ValueError( "shift and extend only available for bigwig output" )

    if options.output_format in ("bigwig", "bigbed"):
        
        if not options.output_filename:
            raise ValueError("please output file for bigwig/bigbed computation.")

        if options.output_format == "bigwig":
            executable_name = "wigToBigWig"
        elif options.output_format == "bigbed":
            executable_name = "bedToBigBed"
        else:
            raise ValueError("unknown output format `%s`" % options.output_format)

        executable = IOTools.which( executable_name )

        if not executable:
            raise OSError( "could not find %s in path." % executable_name )

        tmpdir = tempfile.mkdtemp()
        E.debug( "temporary files are in %s" % tmpdir)

        tmpfile_wig = os.path.join( tmpdir, "wig" )
        tmpfile_sizes = os.path.join( tmpdir, "sizes" )

        # write contig sizes
        outfile_size = open( tmpfile_sizes, "w")
        for contig, size in contig_sizes.items():
            outfile_size.write("%s\t%s\n" % (contig, size) )
        outfile_size.close()    
        
        outfile = open( tmpfile_wig, "w" )
        E.info( "starting output to %s" % tmpfile_wig )

    else:
        outfile = options.stdout
        E.info( "starting output to stdout" )        

    if options.output_format in ("wiggle", "bigwig"):
        # wiggle is one-based, so add 1, also step-size is 1, so need to output all bases
        outf = lambda outfile, contig, start, end, val: \
            outfile.write( "".join( [ "%i\t%i\n" % (x, val) for x in xrange(start+1,end+1) ] ) )
    elif options.output_format in ("bed", "bigbed"):
        # bed is 0-based, open-closed
        outf = lambda outfile, contig, start, end, val: \
            outfile.write("%s\t%i\t%i\t%i\n" % (contig, start, end,val))

    output_filename = os.path.abspath( options.output_filename )

    ninput, nskipped, ncontigs = 0, 0, 0

    if options.shift > 0 or options.extend > 0:

        # create bed file with shifted tags
        shift, extend = options.shift, options.extend
        shift_extend = shift + extend

        for contig in samfile.references:
            E.debug("output for %s" % contig )
            lcontig = contig_sizes[contig]
            
            for read in samfile.fetch( contig ):
                pos = read.pos
                if read.is_reverse:
                    start = max(0, read.pos + read.alen - shift_extend )
                else: 
                    start = max(0, read.pos + shift)

                # intervals extending beyond contig are removed
                if start >= lcontig: continue

                end = min( lcontig, start + extend )
                outfile.write( "%s\t%i\t%i\n" % (contig, start, end))
                
        outfile.close()

        tmpfile_bed = os.path.join( tmpdir, "bed" )
        E.info("computing coverage")
        # calculate coverage - format is bedgraph
        statement = "genomeCoverageBed -bg -i %(tmpfile_wig)s -g %(tmpfile_sizes)s > %(tmpfile_bed)s" % locals()
        E.run( statement )
        
        E.info("converting to bigwig" )

        tmpfile_sorted = os.path.join( tmpdir, "sorted" )
        statement = "sort -k 1,1 -k2,2n %(tmpfile_bed)s > %(tmpfile_sorted)s; bedGraphToBigWig %(tmpfile_sorted)s %(tmpfile_sizes)s %(output_filename)s" % locals()
        E.run( statement)

        shutil.rmtree( tmpdir )

    else:
        for contig in samfile.references:
            # if contig != "chrX": continue
            E.debug("output for %s" % contig )
            lcontig = contig_sizes[contig]

            if options.output_format in ("wiggle", "bigwig"):
                outfile.write( "variableStep chrom=%s span=1\n" % contig )

            for v, iter in itertools.groupby( samfile.pileup(contig), lambda x: x.n ):
                l = list(iter)
                # print l[0].pos, l[-1].pos, l[0].n
                start, end, val = l[0].pos,l[-1].pos+1,l[0].n
                # patch: there was a problem with bam files and reads overextending at the end.
                # These are usually Ns, but need to check as otherwise wigToBigWig fails.
                if lcontig <= end: 
                    E.warn( "read extending beyond contig: %s: %i > %i" % (contig, end, lcontig))
                    end = lcontig 
                    if start >= end: continue

                if val > 0: 
                    outf( outfile, contig, start, end, val)
            ncontigs += 1

        E.info( "finished output" )

        E.info( "ninput=%i, ncontigs=%i, nskipped=%i" % (ninput, ncontigs, nskipped) )

        if options.output_format in ("bigwig", "bigbed"):
            outfile.close()

            E.info( "starting %s conversion" % executable )
            try:
                retcode = subprocess.call( " ".join( (executable,
                                                      tmpfile_wig,
                                                      tmpfile_sizes,
                                                      output_filename )),
                                           shell=True)
                if retcode != 0:
                    E.warn( "%s terminated with signal: %i" % (executable, -retcode))
                    return -retcode
            except OSError, msg:
                E.warn( "Error while executing bigwig: %s" % e)
                return 1

            shutil.rmtree( tmpdir )

            E.info( "finished bigwig conversion" )

    E.Stop()


if __name__ == "__main__":
    sys.exit( main( sys.argv) )

