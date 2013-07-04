################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: run_nubiscan.py 2861 2010-02-23 17:36:32Z andreas $
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

:Author: Andreas Heger
:Release: $Id: run_nubiscan.py 2861 2010-02-23 17:36:32Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

""" 

import os
import sys
import re
import optparse
import numpy
import CGAT.Experiment as E
import CGAT.Nubiscan as Nubiscan
import CGAT.FastaIterator as FastaIterator
import CGAT.Stats as Stats
import CGAT.Genomics as Genomics
import CGAT.Masker as Masker

# from MEME
RXRVDR = numpy.matrix( (( 0.02,  0.16,  0.82,  0.66,  0.00,  0.00 ),
                        ( 0.10,  0.00,  0.12,  0.34,  1.00,  0.32 ),
                        ( 0.02,  0.84,  0.04,  0.00,  0.00,  0.00 ),
                        ( 0.86,  0.00,  0.02,  0.00,  0.00,  0.68 ),
                        (    0,     0,     0,     0,     0,     0 ),
                        ) )

# from Bioprospector - first half site
RXRVDR1 = numpy.matrix( ((  0.09,  0.22,  0.76,  0.70,  0.00,  0.00 ),
                        ( 0.15,  0.00,  0.17,  0.30,  1.00,  0.36 ),
                        ( 0.00,  0.77,  0.00,  0.00,  0.00,  0.00 ),
                        ( 0.76,  0.00,  0.07,  0.00,  0.00,  0.64 ),
                        (    0,     0,     0,     0,     0,     0 ),
                        ) )


# from Bioprospector - second half site
RXRVDR2 = numpy.matrix( ((  0.02,  0.02,  0.88,  0.27,  0.05,  0.03 ),
                        ( 0.05,  0.15,  0.07,  0.73,  0.92,  0.42 ),
                        ( 0.00,  0.83,  0.04,  0.00,  0.00,  0.00 ),
                        ( 0.93,  0.00,  0.01,  0.00,  0.04,  0.55 ),
                        (    0,     0,     0,     0,     0,     0 ),
                        ) )

# from Bioprospector - combined half-sites
RXVDR = (RXRVDR1 + RXRVDR2) / 2.0

NR =  numpy.matrix( (( 0.72,     0,  0.09,  0.27,     0,  0.90 ),
                     (    0,     0,     0,     0,  0.90,  0.09 ),
                     ( 0.27,  0.90,  0.54,  0.09,  0.09,     0 ),
                     (    0,  0.09,  0.36,  0.63,     0,     0 ),
                     (    0,     0,     0,     0,     0,     0 ),
                     ) )


def maskSequences( sequences, masker ):

    if masker == "repeatmasker":
        # the genome sequence is repeat masked
        masked_seq = sequences
    elif masker == "dust":
        masker_object = Masker.MaskerDustMasker()
        masked_seq = [ masker_object( x.upper() ) for x in sequences ]
    else:
        masked_seq = [x.upper() for x in sequences ]

    # hard mask softmasked characters
    masked_seq = [re.sub( "[a-z]","N", x) for x in masked_seq ]
    return masked_seq

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: run_nubiscan.py 2861 2010-02-23 17:36:32Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-i", "--iterations", dest="iterations", type="int",
                      help="number of iterations for sampling [default=%default]."  )

    parser.add_option("-q", "--qvalue", dest="qvalue_threshold", type="float",
                      help="qvalue threshold [default=%default]."  )

    parser.add_option("--without-combine", dest="combine", action = "store_false",
                      help="combine overlapping motifs [default=%default]."  )

    parser.add_option("-f", "--fdr-control", dest="fdr_control", type="choice",
                      choices = ("per-sequence", "all", "xall"),
                      help="qvalue threshold [default=%default]."  )

    parser.add_option("-m", "--motif", dest="motif", type="choice",
                      choices=("rxrvdr", "rxrvdr1", "rxrvdr2", "nr"),
                      help="qvalue threshold [default=%default]."  )

    parser.add_option("-a", "--arrangements", dest="arrangements", type="string",
                      help ="',' separated list of repeat arrangements [default=%default]")

    parser.add_option("-x", "--mask", dest="mask", type="choice",
                      choices=("dust","repeatmasker"),
                      help ="mask sequences before scanning [default=%default]")

    parser.add_option("--output-stats", dest="output_stats", action = "store_true",
                      help="output stats [default=%default]."  )

    parser.add_option("--add-sequence", dest="add_sequence", action = "store_true",
                      help="add sequence information [default=%default]."  )

    parser.set_defaults(
        iterations = 100,
        qvalue_threshold = 0.05,
        motif = "rxrvdr",
        fdr_control = "all",
        combine = True,
        arrangements = None,
        mask = None,
        output_stats = False,
        add_sequence = False,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

    ## do sth
    ninput, nskipped, noutput = 0, 0, 0

    if options.arrangements == None:
        options.arrangements = [ "DR%s" % x for x in range(0,15) ] + [ "ER%s" % x for x in range(0,15) ]
    else:
        options.arrangements = options.arrangements.split(",")
        
    options.stdout.write( "%s" % "\t".join(Nubiscan.NubiscanMatch._fields) )
    if options.add_sequence: options.stdout.write( "\tsequence" )
    options.stdout.write("\n")

    if options.motif == 'nr': sense_matrix = NR
    elif options.motif == "rxrvdr": sense_matrix = RXRVDR
    elif options.motif == "rxrvdr1": sense_matrix = RXRVDR1
    elif options.motif == "rxrvdr2": sense_matrix = RXRVDR2
    else:
        raise ValueError("unknown matrix %s" % options.motif)

    if options.fdr_control == "all":

        seqs = list(FastaIterator.iterate(options.stdin))

        if options.mask:
            masked_seqs = maskSequences( [x.sequence for x in seqs], options.mask )
        else:
            masked_seqs = [x.sequence for x in seqs]
            
        ninput = len(seqs)
        map_id2title = dict( enumerate( [re.sub("\s.*", "", x.title) for x in seqs] ) )
        matcher = Nubiscan.MatcherRandomisationSequences( sense_matrix,
                                                          samples = options.iterations )
        
        results = matcher.run( masked_seqs,
                               options.arrangements,
                               qvalue_threshold = options.qvalue_threshold )

        if options.combine:
            results =  Nubiscan.combineMotifs( results )
        
        for r in results:

            if r.alternatives:
                alternatives = ",".join( [x.arrangement for x in r.alternatives ] )
            else:
                alternatives = ""

            options.stdout.write( "\t".join( ( 
                map_id2title[r.id],
                "%i" % r.start,
                "%i" % r.end,
                r.strand,
                r.arrangement,
                "%6.4f" % r.score,
                "%6.4f" % r.zscore,
                "%6.4e" % r.pvalue,
                "%6.4e" % r.qvalue,
                alternatives) ) )

            if options.add_sequence:
                s = masked_seqs[int(r.id)][r.start:r.end]
                if r.strand == "-": s = Genomics.complement( s )
                s = s[:6].upper() + s[6:-6].lower() + s[-6:].upper()
                options.stdout.write( "\t%s" % s )
                
            options.stdout.write("\n")
            noutput += 1

        # output stats
        if options.output_stats:
            outfile = E.openOutputFile( "fdr" )
            outfile.write("bin\thist\tnobserved\n" )
            for bin, hist, nobs in zip(matcher.bin_edges, matcher.hist, matcher.nobservations):
                outfile.write( "%f\t%f\t%f\n" % (bin, hist, nobs))
            outfile.close()


    elif options.fdr_control == "xall":

        matcher = Nubiscan.MatcherRandomisationSequence( sense_matrix,
                                                         samples = options.iterations )
    

        # collect all results
        matches = []
        for seq in FastaIterator.iterate(options.stdin):
            ninput += 1
            mm = matcher.run( seq.sequence,
                              options.arrangements,
                              qvalue_threshold = None )
            for m in mm:
                matches.append( m._replace( sequence = seq.title ) )

        # estimate qvalues for all matches across all sequences
        pvalues = [ x.pvalue for x in matches ]
        fdr = Stats.doFDR( pvalues )
        qvalues = fdr.mQValues
        results = []
        for m, qvalue in zip(matches, qvalues):
            if qvalue > options.qvalue_threshold: continue
            results.append( m._replace( qvalue = qvalue ) )

        if options.combine:            
            results =  Nubiscan.combineMotifs( results )

        # output
        for r in results:
            options.stdout.write( "\t".join( ( 
                r.id,
                "%i" % r.start,
                "%i" % r.end,
                r.strand,
                r.arrangement,
                "%6.4f" % r.score,
                "%6.4f" % r.zscore,
                "%6.4e" % r.pvalue,
                "%6.4e" % r.qvalue ) ) + "\n" )

            noutput += 1

    elif options.fdr_control == "per-sequence":
        matcher = Nubiscan.MatcherRandomisationSequence( sense_matrix,
                                                         samples = options.iterations )
    

        for seq in FastaIterator.iterate(options.stdin):
            ninput += 1
            result = matcher.run( seq.sequence,
                                  options.arrangements,
                                  qvalue_threshold = options.qvalue_threshold )
            
            if options.combine:
                result =  Nubiscan.combineMotifs( result )

            t = re.sub(" .*","",  seq.title)
            for r in result:
                options.stdout.write( "\t".join( ( 
                    t,
                    "%i" % r.start,
                    "%i" % r.end,
                    r.strand,
                    r.arrangement,
                    "%6.4f" % r.score,
                    "%6.4f" % r.zscore,
                    "%f" % r.pvalue,
                    "%f" % r.qvalue ) ) + "\n" )

            noutput += 1
    
    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput,nskipped) )

    ## write footer and output benchmark information.
    E.Stop()


        

if __name__ == "__main__":
    sys.exit( main( sys.argv) )



