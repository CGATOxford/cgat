################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
runExpression.py - wrap various differential expression tools
=============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script provides a convenience wrapper for differential expression analysis 
for a variety of methods.

Methods implemented are:

   DESeq
   EdgeR

The aim of this script is to provide a common tabular output format that is consistent
between the different methods. The columns in the output are:

+--------------+------------------------------------------------------+
|*Column name* |*Content*                                             |
+--------------+------------------------------------------------------+
|test_id       |Name of the test (gene name, ...                      |
+--------------+------------------------------------------------------+
|treatment_name|Name of the treatment condition                       |
+--------------+------------------------------------------------------+
|treatment_mean|Estimated expression value for treatment              |
+--------------+------------------------------------------------------+
|treatment_std |Standard deviation                                    |
+--------------+------------------------------------------------------+
|control_name  |Name of the control condition                         |
+--------------+------------------------------------------------------+
|control_mean  |Estimated expression value for control                |
+--------------+------------------------------------------------------+
|control_std   |Standard deviation                                    |
+--------------+------------------------------------------------------+
|pvalue        |The p value for rejecting the null hypothesis         |
+--------------+------------------------------------------------------+
|qvalue        |Multiple testing correction                           |
+--------------+------------------------------------------------------+
|l2fold        |log2 foldchange of treatment/control                  |
+--------------+------------------------------------------------------+
|fold          |foldchange of treatment/control                       |
+--------------+------------------------------------------------------+
|significant   |Flag, 1 if test called significant according to FDR   |
+--------------+------------------------------------------------------+
|status        |test status (OK|FAIL)                                 |
+--------------+------------------------------------------------------+

The script will call each of the method and output a variety of diagnostic plots.

Usage
-----

Documentation
-------------

Command line options
--------------------

'''

import math
import numpy
import sys, os, optparse
import collections
import itertools

from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri

try:
    import CGAT.Experiment as E
    import CGAT.Pipeline as P
    import CGAT.Database as Database
    import CGAT.IOTools as IOTools
    import CGAT.Stats as Stats
    import CGAT.Expression as Expression
except ImportError:
    import Experiment as E
    import Pipeline as P
    import Database,IOTools,Stats
    import Expression

import sqlite3

try:
    PARAMS = P.getParameters()
except IOError:
    pass


def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-t", "--filename-tags", dest="input_filename_tags", type="string",
                      help="input file with tag counts [default=%default]."  )

    parser.add_option("-d", "--filename-design", dest="input_filename_design", type="string",
                      help="input file with experimental design [default=%default]."  )

    parser.add_option("-o", "--outfile", dest="output_filename", type="string",
                      help="output filename [default=%default]."  )

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices = ("deseq", "edger", "cuffdiff", "summary" ),
                      help="differential expression method to apply [default=%default]."  )

    parser.add_option( "--deseq-dispersion-method", dest="deseq_dispersion_method", type="choice",
                      choices = ("pooled", "per-condition", "blind"),
                      help="dispersion method for deseq [default=%default]."  )

    parser.add_option( "--deseq-fit-type", dest="deseq_fit_type", type="choice",
                      choices = ("parametric", "local"),
                      help="fit type for deseq [default=%default]."  )

    parser.add_option("-f", "--fdr", dest="fdr", type="float",
                      help="fdr to apply [default=%default]."  )

    parser.add_option("-r","--reference-group", dest="ref_group", type="string",
                      help="Group to use as reference to compute fold changes against [default=$default]")


    parser.set_defaults(
        input_filename_tags = "-",
        input_filename_design = None,
        output_filename = sys.stdout,
        method = "deseq",
        fdr = 0.1,
        deseq_dispersion_method = "pooled",
        deseq_fit_type = "local",
        ref_group = None
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

    if options.input_filename_tags == "-":
        fh = P.getTempFile()
        fh.write( "".join( [ x for x in options.stdin ] ) )
        fh.close()
        options.input_filename_tags = fh.name
    else:
        fh = None

    if options.method == "deseq":
        assert options.input_filename_tags and os.path.exists(options.input_filename_tags)
        assert options.input_filename_design and os.path.exists(options.input_filename_design)
        Expression.runDESeq( options.input_filename_tags,
                             options.input_filename_design,
                             options.output_filename,
                             options.output_filename_pattern,
                             fdr = options.fdr,
                             dispersion_method = options.deseq_dispersion_method,
                             fit_type = options.deseq_fit_type,
                             ref_group = options.ref_group)

    elif options.method == "edger":
        assert options.input_filename_tags and os.path.exists(options.input_filename_tags)
        assert options.input_filename_design and os.path.exists(options.input_filename_design)
        Expression.runEdgeR( options.input_filename_tags,
                             options.input_filename_design,
                             options.output_filename,
                             options.output_filename_pattern,
                             fdr = options.fdr,
                             ref_group = options.ref_group)

    elif options.method == "summary":
        assert options.input_filename_tags and os.path.exists(options.input_filename_tags)
        Expression.outputTagSummary( options.input_filename_tags,
                                     options.stdout,
                                     options.output_filename_pattern )

    if fh and os.path.exists( fh.name): os.unlink( fh.name )

    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
    
