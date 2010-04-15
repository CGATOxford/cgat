################################################################################
#   Gene prediction pipeline 
#
#   $Id: evaluate_clusters.py 2781 2009-09-10 11:33:14Z andreas $
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
import os, sys, string, re, getopt, time, sets, optparse, math, tempfile

import Experiment
import csv

parser = optparse.OptionParser( version = "%prog version: $Id: evaluate_clusters.py 2781 2009-09-10 11:33:14Z andreas $")

def ConvertDictionary( d ):
    """tries to convert values in a dictionary.
    """

    rx_int = re.compile("^[+-]*[0-9]+$")
    rx_float = re.compile("^[+-]*[0-9.+-eE]+$")
    for k,v in d.items():

        if rx_int.match( v ):
            d[k] = int(v)
        elif rx_float.match( v ):
            d[k] = float(v)

    return d
    
if __name__ == "__main__":

    parser.add_option("-c", "--file-components", dest="filename_components", type="string",
                      help="filename to write compoments into." )

    parser.set_defaults(
        file_headers = None,
        min_identical_exons = 0.50,
        min_num_taxa = 1,
        max_num_fragments = 0,
        max_num_nonoverlaps = 0,
        min_percent_conserved_perfect = 0.0,
        min_percent_conserved_partially = 0.0,
        min_percent_identical_exons = 0.5,
        )

    (options, args) = Experiment.Start( parser )

    reader = csv.DictReader( sys.stdin, dialect='excel-tab' )
    
    for row in reader:
        
        row = ConvertDictionary( row )

        
        

        
    Experiment.Stop()
