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
MACS.py - Parser for MACS output
================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

The :mod:`Pipeline` module contains various utility functions
for parsing MACS output.

API
----

'''

import collections

MacsPeak = collections.namedtuple( "MacsPeak", "contig start end length summit tags pvalue fold fdr" )

def iterateMacsPeaks( infile ):
    '''iterate over peaks.xls file and return parsed data.
    The fdr is converted from percent to values between 0 and 1.
    '''
 
    for line in infile:
        if line.startswith("#"): continue
        if line.startswith("chr\tstart"): continue
        # skip empty lines
        if line.startswith("\n"): continue

        data = line[:-1].split("\t")

        if len(data) == 9:
            # convert % fdr 
            data[8] = float(data[8]) / 100.0
        elif len(data) == 8:
            ## if no fdr given, set to 0
            # data.append( 0.0 )
            ## Steve - I don't understand this so I'm commenting it out and raising an error
            raise ValueError( "FDR value not set line %s" % line )
        else:
            raise ValueError( "could not parse line %s" % line )

        # these are 1-based coordinates
        # macs can have negative start coordinates
        # start
        data[1] = max(int(data[1]) - 1, 0)
        # end
        data[2] = int(data[2])
        # length
        data[3] = int(data[3])
        # summit
        data[4] = int(data[4])
        # ntags
        data[5] = int(data[5])
        # -10log10(pvalue)
        data[6] = float(data[6])
        # fold
        data[7] = float(data[7])
        
        yield MacsPeak._make( data )

Macs2Peak = collections.namedtuple( "Macs2Peak", "contig start end length summit pileup pvalue fold fdr name" )

def iterateMacs2Peaks( infile ):
    '''iterate over peaks.xls file and return parsed data.
    The fdr is converted from percent to values between 0 and 1.
    '''
 
    for line in infile:
        if line.startswith("#"): continue
        if line.startswith("chr\tstart"): continue
        # skip empty lines
        if line.startswith("\n"): continue

        data = line[:-1].split("\t")

        if len(data) == 10:
            # convert % fdr 
            data[8] = 10**-float(data[8])
        elif len(data) == 8:
            # if no fdr given, set to 0
            #data.append( 0.0 )
            raise ValueError( "No FDR found line %s" % line )
        else:
            raise ValueError( "could not parse line %s" % line )

        # these are 1-based coordinates
        # macs can have negative start coordinates
        # start
        data[1] = max(int(data[1]) - 1, 0)
        # end
        data[2] = int(data[2])
        # length
        data[3] = int(data[3])
        # summit
        data[4] = int(data[4])
        # ntags
        data[5] = float(data[5])
        # -10log10(pvalue)
        data[6] = float(data[6])
        # fold
        data[7] = float(data[7])
        
        yield Macs2Peak._make( data )

