################################################################################
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
#################################################################################
'''
psl2table.py - output stats for psl formatted alignments
================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python psl2table.py --help

Type::

   python psl2table.py --help

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
import tempfile
import subprocess
import optparse
import math

USAGE = \
"""analyze sequence pairs from a psl formatted table.

The sequences are assumed to be nucleotide sequences.

Methods available are:

baseml: compute baseml rates 
counts: compute residue counts (percent G+C, ...)
match:  compute match statistics (pid, coverage)
"""

import CGAT.Experiment as E
import CGAT.Blat as Blat
import CGAT.SequenceProperties as SequenceProperties
import CGAT.SequencePairProperties as SequencePairProperties
import CGAT.IOTools as IOTools
import CGAT.WrapperCodeML as WrapperCodeML

##---------------------------------------------------------------------------------------------
class Counter:
    def __init__( self, options):
        pass

class CounterMatch( Counter ):
    
    def __init__(self,*args,**kwargs):
        Counter.__init__(self,*args,**kwargs)

    def __call__(self, match):
        
        self.mPid = 100.0 * match.mNMatches / (match.mNMatches + match.mNMismatches )
        self.mQueryCoverage = 100.0 * (match.mNMatches + match.mNMismatches) / match.mQueryLength

    def __str__( self ):
        return "%6.4f\t%6.4f" % (self.mPid, self.mQueryCoverage)

    def getHeaders(self):
        return [ "pid", "qCov" ]

class QueriesCounter( SequenceProperties.SequencePropertiesNA ):
    def __init__(self, *args, **kwargs):
        SequenceProperties.SequencePropertiesNA.__init__(self, *args, **kwargs)

    def __call__( self, seq1, seq2 ):
        SequenceProperties.SequencePropertiesNA.loadSequence( self, seq1 )

class SbjctsCounter( SequenceProperties.SequencePropertiesNA ):
    def __init__(self, *args, **kwargs):
        SequenceProperties.SequencePropertiesNA.__init__(self, *args, **kwargs)

    def __call__( self, seq1, seq2 ):
        SequenceProperties.SequencePropertiesNA.loadSequence( self, seq2 )



##---------------------------------------------------------------------------------------------
    
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: psl2table.py 2891 2010-04-07 08:59:18Z andreas $",
                                    usage = globals()["__doc__"])

    parser.add_option( "--mask-lowercase", dest="mask_lowercase", action="store_true",
                      help="mask lowercase characters before computing properties [default=%default]" )

    parser.add_option( "--with-match", dest="with_match", action="store_true",
                      help="echo the match in output [default=%default]" )

    parser.add_option( "--without-match", dest="with_match", action="store_false",
                      help="do not echo the match in output [default=%default]" )

    parser.add_option("-m", "--method", dest="methods", type="choice", action="append",
                      choices=("counts", "baseml", "match", "query-counts", "sbjct-counts" ),
                      help="methods to compute properties between sequence pairs." )

    WrapperCodeML.BaseML().AddOptions( parser )

    parser.set_defaults(
        methods = [],
        mask_lowercase = False,
        is_pslx = True,
        with_match = True,
        )
    
    (options, args) = E.Start( parser )

    counters_plain = []
    counters = []

    for method in options.methods:
        if method == "counts":
            counters.append( SequencePairProperties.SequencePairPropertiesCountsNa() )
        elif method == "query-counts":
            counters.append( QueriesCounter() )
        elif method == "sbjct-counts":
            counters.append( SbjctsCounter() )
        elif method == "baseml":
            counters.append( SequencePairProperties.SequencePairPropertiesBaseML( options ) )
        elif method == "match":
            counters_plain.append( CounterMatch( options ) )
            
    if counters:
        iterator = Blat.iterator_pslx( options.stdin )
        header = "\t".join(Blat.MatchPSLX().getHeaders())
    else:
        iterator = Blat.iterator( options.stdin )
        header = "\t".join(Blat.Match().getHeaders())

    if not options.with_match:
        header = "qName"

    options.stdout.write( "\t".join( 
            [header,] + 
            [ "\t".join(x.getHeaders()) for x in counters] +
            [ "\t".join(x.getHeaders()) for x in counters_plain] ) + "\n" )

    ninput, noutput, nskipped = 0, 0, 0

#     ## setup totals
#     totals = {}
#     for section in options.sections:
#         if section == "length":
#             s = SequencePropertiesLength()
#         elif section == "na":
#             s = SequencePropertiesNA()
#         elif section == "aa":
#             s = SequencePropertiesAA()
#         elif section == "degeneracy":
#             s = SequencePropertiesDegeneracy()
#         elif section == "bias":
#             s = SequencePropertiesBias( reference_codons )
#         elif section == "codons":
#             s = SequencePropertiesCodons()
#         elif section == "codon-usage":
#             s = SequencePropertiesCodonUsage()
#         elif section == "codon-translator":
#             s = SequencePropertiesCodonTranslator()
#         else:
#             raise "unknown section %s" % section
        
#         totals[section] = s


    for match in iterator:
        ninput += 1

        if options.with_match:
            options.stdout.write( str(match) )
        else:
            options.stdout.write( match.mQueryId )
            
        if counters:
        
            qseq = match.mQuerySequence 
            sseq = match.mSbjctSequence 

            # mask non printable characters - sometimes
            # appear after using pslToPslX
            qseq = [ re.sub( "[^a-zA-Z]", "N", x ) for x in qseq ]
            sseq = [ re.sub( "[^a-zA-Z]", "N", x ) for x in sseq ]

            if options.mask_lowercase:
                qseq = [ re.sub( "[a-z]", "N", x) for x in qseq ]
                sseq = [ re.sub( "[a-z]", "N", x) for x in sseq ]
                
            match.mQuerySequence = qseq
            match.mSbjctSequence = sseq

            qseq = "".join( match.mQuerySequence ).upper()
            sseq = "".join( match.mSbjctSequence ).upper()

            if len(qseq) != len(sseq):
                if options.loglevel >= 1:
                    options.stdlog.write( "# WARNING: two sequences of unequal length in match\n# %s\n" % str(match) )
                nskipped += 1
                continue

            for counter in counters: counter( qseq, sseq )

            options.stdout.write( "\t" + 
                                  "\t".join( \
                    [str( counter ) for counter in counters] ) )
        
        if counters_plain:

            for counter in counters_plain: counter( match )

            options.stdout.write( "\t" + 
                                  "\t".join( \
                    [str( counter ) for counter in counters_plain] ) )

        options.stdout.write( "\n" )

        noutput += 1

    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nskipped=%i\n" % (ninput, noutput, nskipped) )

    E.Stop()
