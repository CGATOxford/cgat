################################################################################
#   Gene prediction pipeline 
#
#   $Id: patch_translations.py 1841 2008-05-08 12:07:13Z andreas $
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

USAGE="""recalculated translations from a predictions file.
"""

import sys, os, string, re, optparse

import Experiment
import PredictionParser
import IndexedFasta
import Genomics

##------------------------------------------------------------
if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id: patch_translations.py 1841 2008-05-08 12:07:13Z andreas $",
                                    usage = USAGE )

    parser.add_option( "-g", "--genome-file", dest="genome_file", type="string",
                       help="filename with genome."  )

    parser.set_defaults(
        genome_file = None,
        )
    
    (options, args) = Experiment.Start( parser, add_pipe_options = True )

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(2)

    entry = PredictionParser.PredictionParserEntry()

    fasta = IndexedFasta.IndexedFasta( options.genome_file )

    ninput, noutput = 0, 0
    for line in sys.stdin:
        if line[0] == "#":
            print line[:-1]
            continue
        
        entry.Read( line )
        
        ninput += 1
        
        ## get genomic sequence
        genomic_sequence = fasta.getSequence( entry.mSbjctToken,
                                              entry.mSbjctStrand,
                                              entry.mSbjctGenomeFrom,
                                              entry.mSbjctGenomeTo )

        entry.mMapPeptide2Translation, entry.mTranslation = Genomics.Alignment2PeptideAlignment( \
               entry.mMapPeptide2Genome, entry.mQueryFrom, 0, genomic_sequence )

        options.stdout.write( str(entry) + "\n" )

        noutput += 1

    if options.loglevel >= 1:
        options.stdlog.write("# ninput=%i, noutput=%i\n" % (ninput, noutput))
        

    Experiment.Stop()
