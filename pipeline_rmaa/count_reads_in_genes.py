## part of pipeline_rmaa.py

import sys, optparse, itertools
import pysam
import Bed, IOTools

import Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    coords_file=args[0]

    bamfile=pysam.Samfile( args[1], 'rb' )  # bamfile

    options.stdout.write( "gene_id\tcounts\tlength\n" )

    iter = Bed.iterator( IOTools.openFile( coords_file ) )
    for gene_id, exons in itertools.groupby( iter, lambda x: x.name ):

        num_reads=0
        
        anames=set([])
        lgene = 0

        for bed in exons:
            lgene += bed.end - bed.start
            for alignedread in bamfile.fetch(bed.contig, bed.start, bed.end):
                anames.add((alignedread.qname, alignedread.is_read1))

        num_reads = len(anames)
        options.stdout.write( "\t".join( (gene_id,
                                          str(num_reads),
                                          str(lgene ) )) + "\n" )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
