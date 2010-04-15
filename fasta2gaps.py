################################################################################
#   Gene prediction pipeline 
#
#   $Id: fasta2gaps.py 2781 2009-09-10 11:33:14Z andreas $
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
import os, sys, string, re, getopt, tempfile, time, optparse, math, glob

import Experiment
import IndexedFasta

USAGE="""python %s [OPTIONS] 

Version: $Id: fasta2gaps.py 2781 2009-09-10 11:33:14Z andreas $

extract sequences/sequence regions from a fasta file.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
""" % sys.argv[0]


if __name__ == "__main__":

    parser = optparse.OptionParser( version = "%prog version: $Id: fasta2gaps.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="pattern to look for sequence filename."  )

    parser.set_defaults(
        genome_file = None,
        gap_char = "NnXx"
        )

    (options, args) = Experiment.Start( parser )

    fasta = IndexedFasta.IndexedFasta( options.genome_file )
    contigs = fasta.getContigSizes()

    options.stdout.write( "contig\tstart\tend\n" )

    for contig in contigs.keys():

        s = fasta.getSequence( contig, "+", 0, 0 )
        
        first_res = -1
        x = 0
        xx = len(s)
        while x < xx:
            if s[x] not in options.gap_char:
                if first_res > 0:
                    options.stdout.write( "%s\t%i\t%i\n" % (contig, first_res, x))
                first_res = -1
            else:
                if first_res < 0:
                    first_res = x

            x += 1
            
        if first_res > 0:
            options.stdout.write( "%s\t%i\t%i\n" % (contig, first_res, x))

    Experiment.Stop()
