'''
bed2fasta.py - get sequences from bed file
==========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics Intervals Sequences Conversion BED FASTA 

Purpose
-------

This script outputs nucleotide sequences for intervals within
a :term:`bed` formatted file using.

Usage
-----

A required input to bed2fasta.py is a CGAT indexed genome. To obtain an
idexed human reference genome we would type

Example::
   cat hg19.fasta | index_fasta.py hg19 > hg19.log

This file would then serve as the --genome-file when we wish to extract
sequences from a :term:`bed` formatted file.


For example we could now type::

   cat in.bed | python bed2fasta.py --genome-file hg19 > out.fasta

Where we take a set of genomic intervals (e.g. from a human ChIP-seq experiment)
and output their respective nucleotide sequences.


Type::

   python bed2fasta.py --help

for command line help.

Command line options
--------------------

'''
import sys
import string
import re
import optparse
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Bed as Bed
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Masker as Masker

############################################################
############################################################
############################################################
def maskSequences( sequences, masker = None):
    '''return a list of masked sequence.

    *masker* can be one of
        dust/dustmasker * run dustmasker on sequences
        softmask        * use softmask to hardmask sequences
    '''

    if masker in ("dust", "dustmasker"):
        masker_object = Masker.MaskerDustMasker()
    else:
        masker_object = None

    if masker == "softmask":
        # the genome sequence is repeat soft-masked
        masked_seq = sequences
    elif masker in ("dust", "dustmasker"):
        # run dust
        masked_seq = masker_object.maskSequences( [ x.upper() for x in sequences ] )
    elif masker == None:
        masked_seq = [x.upper() for x in sequences ]
    else:
        raise ValueError("unknown masker %s" % masker )

    # hard mask softmasked characters
    masked_seq = [re.sub( "[a-z]","N", x) for x in masked_seq ]

    return masked_seq

##------------------------------------------------------------------------
def main( argv = None ):
    if argv == None: argv = sys.argv

    parser = E.OptionParser( version = "%prog version: $Id: gff2fasta.py 2861 2010-02-23 17:36:32Z andreas $")

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome."  )

    parser.add_option("-m", "--masker", dest="masker", type="choice",
                      choices= ("dust", "dustmasker", "softmask", "none" ),
                      help="apply masker [%default]."  )

    parser.add_option( "-o", "--mode", dest="mode", type="choice",
                       choices = ("intervals", "leftright" ),
                       help="what to output [%default]" )

    parser.add_option( "--min-length", dest="min_length", type="int",
                       help="require a minimum sequence length [%default]" )

    parser.add_option( "--max-length", dest="max_length", type="int",
                       help="require a maximum sequence length [%default]" )

    parser.add_option( "--extend-at", dest="extend_at", type="choice",
                       choices=("none", "3", "5", "both", "3only", "5only" ),
                       help="extend at no, 3', 5' or both ends. If 3only or 5only are set, only the added sequence is returned [default=%default]" )

    parser.add_option( "--extend-by", dest="extend_by", type="int",
                       help="extend by # bases [default=%default]" )

    parser.add_option( "--use-strand", dest="ignore_strand", action="store_false",
                       help="use strand information and return reverse complement [default=%default]" )


    parser.set_defaults(
        genome_file = None,
        masker = None,
        mode = "intervals",
        min_length = 0,
        max_length = 0,
        extend_at = None,
        extend_by = 100,
        ignore_strand = True,
        )

    (options, args) = E.Start( parser )

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )
        contigs = fasta.getContigSizes()
        fasta.setConverter( IndexedFasta.getConverter( "zero-both-open") )

    counter = E.Counter()
    ids, seqs = [], []

    E.info( "collecting sequences" )
    for bed in Bed.setName( Bed.iterator( options.stdin )):
        counter.input += 1

        lcontig = fasta.getLength( bed.contig )

        if options.ignore_strand:
            strand = "+"
        else:
            strand = bed.strand

        if options.mode == "intervals":
            ids.append( "%s %s:%i..%i (%s)" % (bed.name, bed.contig, bed.start, bed.end, strand) )
            seqs.append( fasta.getSequence( bed.contig, strand, bed.start, bed.end) )

        elif options.mode == "leftright":
            l = bed.end - bed.start

            start, end = max(0,bed.start-l), bed.end-l
            ids.append( "%s_l %s:%i..%i (%s)" % (bed.name, bed.contig, start, end, strand) )
            seqs.append( fasta.getSequence( bed.contig, strand, start, end) )
            
            start, end = bed.start+l, min(lcontig,bed.end+l)
            ids.append( "%s_r %s:%i..%i (%s)" % (bed.name, bed.contig, start, end, strand) )
            seqs.append( fasta.getSequence( bed.contig, strand, start, end) )
            
    E.info( "collected %i sequences" % len(seqs) )

    masked = maskSequences( seqs, options.masker )
    options.stdout.write("\n".join( [ ">%s\n%s" % (x,y) for x,y in zip(ids, masked) ] ) + "\n" )

    E.info( "masked %i sequences" % len(seqs) )

    counter.output = len(seqs )

    E.info( "%s" % counter )

    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )


