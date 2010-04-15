################################################################################
#   Gene prediction pipeline 
#
#   $Id: links2exons.py 1799 2008-03-28 11:44:19Z andreas $
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
import os, sys, string, re, getopt, tempfile

USAGE="""python %s [OPTIONS] < orthologs > genes

Version: $Id: links2exons.py 1799 2008-03-28 11:44:19Z andreas $

Transform a list of blastp alignments between sequences to
a list of alignments between exons.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-m, --map=                      map alignments between peptides to nucleotides
-c, --cds=                      files with coding sequences (exon boundaries)
-e, --expand                    write in nucleic acid coordinates (default: peptide coordinates)
""" % sys.argv[0]

import Experiment
import Genomics
import Exons
import alignlib
import BlastAlignments

param_long_options=["verbose=", "help",
                    "cds=", "map="]

param_short_options="v:hc:m:"

param_loglevel = 2

param_report_step=100000

param_filename_cds = None
param_expand = False

def ScaleAlignment( alignment, factor ):
    """scale alignment string."""

    data = re.split("[+-]", alignment[1:])
    
    data = map( lambda x: int(x) * factor, data )
    signs = [ "+", "-" ] * (1 + len(data) / 2)
    
    if alignment[0] == "+":
        del signs[-1]
    else:
        del signs[0]

    s = map( lambda x,y: "%s%i" % (x,y), signs, data)
    return string.join(s, "")

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
        elif o == "--cds":
            param_filename_cds = a
        elif o in ("-e", "--expand"):
            param_expand = True

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(2)

    print Experiment.GetHeader()
    print Experiment.GetParams()
    sys.stdout.flush()

    if param_loglevel >= 1:
        print "# reading exon boundaries."
        sys.stdout.flush()

    cds = Exons.ReadExonBoundaries( open(param_filename_cds, "r") )

    if param_loglevel >= 1:
        print "# read %i cds" % (len(cds))
        sys.stdout.flush()

    ninput, npairs, nskipped = 0, 0, 0

    map_row2col = alignlib.makeAlignmentVector()
    tmp_map_row2col = alignlib.makeAlignmentVector()    
    
    for line in sys.stdin:
        if line[0] == "#": continue
        ninput += 1
        link = BlastAlignments.Link()

        link.Read(line)
        
        if link.mQueryToken == link.mSbjctToken: continue
        
        if link.mQueryToken in cds and \
               link.mSbjctToken in cds:
            
            ## expand to codons
            if param_expand:
                link.mQueryFrom = (link.mQueryFrom - 1) * 3 + 1
                link.mSbjctFrom = (link.mSbjctFrom - 1) * 3 + 1        
                link.mQueryAli = ScaleAlignment( link.mQueryAli, 3 )
                link.mSbjctAli = ScaleAlignment( link.mSbjctAli, 3 )            

            map_row2col.clear()
            alignlib.AlignmentFormatExplicit( 
                link.mQueryFrom, link.mQueryAli,
                link.mSbjctFrom, link.mSbjctAli ).copy( map_row2col )

            ## test all combinations, the alignment might be a suboptimal alignment in case
            ## of repeats.
            for e1 in cds[link.mQueryToken]:
                for e2 in cds[link.mSbjctToken]:
                    tmp_map_row2col.clear()
                    if param_expand:
                        alignlib.copyAlignment( tmp_map_row2col, map_row2col,
                                               e1.mPeptideFrom + 1, e1.mPeptideTo,
                                               e2.mPeptideFrom + 1, e2.mPeptideTo,                                               
                                               )
                    else:
                        alignlib.copyAlignment( tmp_map_row2col, map_row2col,
                                               e1.mPeptideFrom / 3 + 1, e1.mPeptideTo / 3 + 1,
                                               e2.mPeptideFrom / 3 + 1, e2.mPeptideTo / 3 + 1)

                    ## in case of split codons, there is an alignment of length 1. Skip that.
                    if tmp_map_row2col.getLength() > 1:
                        
                        print string.join( map(str, (link.mQueryToken, e1.mRank,
                                                     link.mSbjctToken, e2.mRank,
                                                     link.mEvalue,
                                                     alignlib.AlignmentFormatEmissions( tmp_map_row2col ))), "\t" )
                        
                        npairs += 1
        else:
            if param_loglevel >= 2: print "# SKIPPED: %s" % str(link)
            nskipped += 1

        if (ninput % param_report_step) == 0:
            if param_loglevel >= 1:
                print "# ninput=%i, noutput=%i, nskipped=%i" % (ninput, npairs, nskipped)
            sys.stdout.flush()

    if param_loglevel >= 1:            
        print "# ninput=%i, noutput=%i, nskipped=%i" % (ninput, npairs, nskipped)
    
    print Experiment.GetFooter()
