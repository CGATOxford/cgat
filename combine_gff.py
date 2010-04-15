################################################################################
#   Gene prediction pipeline 
#
#   $Id: combine_gff.py 2781 2009-09-10 11:33:14Z andreas $
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
import sys, string, re, optparse

USAGE="""python %s [OPTIONS] < in > out

combine gff features as overlapping regions.

Version: $Id: combine_gff.py 2781 2009-09-10 11:33:14Z andreas $

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
""" % sys.argv[0]


import Experiment
import PredictionParser
import GFF

if __name__ == "__main__":

    parser = optparse.OptionParser( version = "%prog version: $Id: combine_gff.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-f", "--format", dest="format",
                      help="output format.", type="choice", choices=("flat", "full", "first") )

    parser.set_defaults(
        format="flat"
        )

    (options, args) = Experiment.Start( parser )

    last_e = None
    for line in sys.stdin:
        if line[0] == "#" : continue
        if options.format in ("full", "first" ):
            last_e = GFF.Entry()
        else:
            last_e = GFF.Entry()
        last_e.Read(line)
        break
    
    for line in sys.stdin:
        
        if line[0] == "#" : continue

        if options.format in ("full", "first" ):
            e = GFF.Entry()
        else:
            e = GFF.Entry()
        e.Read(line)

        if not GFF.Overlap( last_e, e):
            print str(last_e)
            last_e = e
        else:
            last_e.start = min( last_e.start, e.start )
            last_e.end = max( last_e.end, e.end )
            if options.format == "full":
                last_e.mInfo += " ; " + e.mInfo
            continue
        
    print str(last_e)

    Experiment.Stop()

        
