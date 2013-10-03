################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: gff2bed.py 2861 2010-02-23 17:36:32Z andreas $
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
gff2bed.py - convert from gff to bed
====================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics Intervals

Purpose
-------

convert gff or gtf to bed formatted files.

Usage
-----

Type::

   python gff2bed.py --help

for command line usage.

Command line options
--------------------
'''

import sys
import re
import string
import optparse
import time
import os
import itertools
import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.Bed as Bed

def main( argv = sys.argv ):

    parser = E.OptionParser( version = "%prog version: $Id: gff2bed.py 2861 2010-02-23 17:36:32Z andreas $", 
                                    usage=globals()["__doc__"])
    
    parser.add_option( "--is-gtf", dest="is_gtf", action="store_true",
                       help="input file is gtf [default=%default] ")

    parser.add_option( "--name", dest="name", type="choice",
                       help = "field to use as the name field [%default]",
                       choices = ("gene_id", "transcript_id", "class", "family", "feature", "source" ) )
    
    parser.add_option( "--track", dest="track", type="choice",
                       choices=("feature", "source", None),
                       help="use feature/source field to define tracks [default=%default] ")

    parser.set_defaults(
        track = None,
        name = "gene_id",
        is_gtf = False )
    
    (options, args) = E.Start( parser, add_pipe_options = True )

    ninput, noutput = 0, 0
    
    is_gtf, name = options.is_gtf, options.name
    iterator = GTF.iterator( options.stdin )

    if options.track:
        all_input = list( iterator )

        if options.track == "feature":
            grouper = lambda x: x.feature
        elif options.track == "source":
            grouper = lambda x: x.source

        all_input.sort( key = grouper )

        bed = Bed.Bed()
        for key, vals in itertools.groupby( all_input, grouper ):
            options.stdout.write( "track name=%s\n" % key )
            for gff in vals:
                ninput += 1
                bed.fromGTF( gff, is_gtf = is_gtf, name = name )
                options.stdout.write( str(bed) + "\n" )
                noutput += 1

    else:
        bed = Bed.Bed()
        for gff in iterator:
            ninput += 1
            bed.fromGTF( gff, is_gtf = is_gtf, name = name )
            options.stdout.write( str(bed) + "\n" )

            noutput += 1

    E.info( "ninput=%i, noutput=%i" % (ninput, noutput) )
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
