################################################################################
#   Gene prediction pipeline 
#
#   $Id: sequences2graph.py 1799 2008-03-28 11:44:19Z andreas $
#
#   Copyright (C) 2004 Andreas Heger
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
import os, sys, string, re, getopt

USAGE="""python %s < sequences > graph

Version: $Id: sequences2graph.py 1799 2008-03-28 11:44:19Z andreas $

align sequences and put them into a graph.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
""" % sys.argv[0]

import Experiment
import Genomics
import Intervalls
import PredictionParser
import alignlib

param_long_options=["verbose=", "help"]
param_short_options="v:ho:t:"

param_loglevel = 1

param_overlap_residues = 90
param_filename_filter_tokens = None

param_gop = -10.0
param_gep = -2.0

##------------------------------------------------------------
if __name__ == '__main__':

    try:
        optlist, args = getopt.getopt(sys.argv[1:], param_short_options, param_long_options)
    except getopt.error, msg:
        print USAGE, msg
        sys.exit(2)
                                                                                                                                                             
    for o,a in optlist:
        if o in ( "-v", "--verbose" ):
            param_loglevel = int(a)
        elif o in ( "-h", "--help" ):
            print USAGE
            sys.exit(0)

    alignator = alignlib.makeAlignatorDPFull( alignlib.ALIGNMENT_LOCAL, param_gop, param_gep )
    map_query2token = alignlib.makeAlignmentVector()
    
    for line in sys.stdin:
        if line[0] == "#": continue

        query_token, sbjct_token, query_sequence, sbjct_sequence = string.split(line[:-1], "\t")

        map_query2token.clear()
        row = alignlib.makeSequence(query_sequence)
        col = alignlib.makeSequence(sbjct_sequence)
        alignator.align( map_query2token, row, col )

        pidentity = 100.0 * alignlib.calculatePercentIdentity( map_query2token, row, col )
        psimilarity = 100.0 * alignlib.calculatePercentSimilarity( map_query2token )        
        print string.join( map(str, (
            query_token, sbjct_token,
            map_query2token.getScore(),
            alignlib.AlignmentFormatEmissions( map_query2token ),
            pidentity,
            psimilarity,
            map_query2token.getNumGaps()) ), "\t" )
            
            
            
        
