################################################################################
#   Gene prediction pipeline 
#
#   $Id: compare_subset.py 2782 2009-09-10 11:40:29Z andreas $
#
#   Copyright (C) 2006 Andreas Heger
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
import os, sys, string, re, optparse, time, random

"""check if a subset if equally sampled from a background

Input: a sample and a background. These can be either given
as values directly or as categories.

Question: is the distribution of values sampled from Input3
different from the background sample Input2?

Test:           Mann-Whitney U test
Implementation: Use R implemenation

"""

import Experiment
import pgdb

if __name__  == "__main__":

    parser = optparse.OptionParser( version = "%prog version: $Id: compare_subset.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option( "-s", "--schemas", dest="schemas", type="string" ,
                       help="schemas to lookup codon bias from.")

    parser.set_defaults(
        schemas="",
        filename_input_sample = None,
        filename_input_background = None,
        filename_input_map = None,
        )
    
    (options, args) = Experiment.Start( parser,
                                        add_pipe_options = True,
                                        add_psql_options = True,)


    

    
    Experiment.Stop()
