################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
#
#   Copyright (C) 2013 Stephen Sansom
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
run_function.py - run a function within a python module remotely
=================================================================

:Author: Stephen Sansom, with credit to Tim for implementation suggestions.
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script allows a specificied function from a specified python module 
to be executed on the cluster with specified parameters.

Usage
-----

A statement is specified in the normal way i.e.::

   statement = """python %(scriptsdir)s/run_function.py 
                   -p infile,outfile,additional_param1
                   -m modulefile 
                   -f function"""

   P.run()

If the module is in your $PYTHONPATH you can just name it 
directly. i.e "Pipeline" would specifiy Pipeline.py

The script has currently only been tested with single input/output.

Documentation
-------------

Code
----

'''

import sys
import os
import importlib

import CGAT.Experiment as E

def main(argv = None):
    
    # Parse the options
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $",
                             usage = globals()["__doc__"] )

    parser.add_option("-p", "--params", dest="params", type="string",
                      help="comma separated list of addtional parameter strings")

    parser.add_option("-m", "--module", dest="module", type="string",
                      help="the full path to the module file", default = None)

    parser.add_option("-i", "--input", dest="input_filenames", type="string", action="append",
                      help="input filename")

    parser.add_option("-o", "--output", dest="output_filenames", type="string", action="append",
                      help="output filename")

    parser.add_option("-f", "--function", dest="function", type="string",
                      help="the module function", default = None)
    
    parser.set_defaults( 
        input_filenames = [],
        output_filenames = [],
        params = None
        )

    (options, args) = E.Start( parser )

    # Check a module and function have been specified
    if not options.module or not options.function:
        raise ValueError("Both a function and Module must be specified")

    # If a full path was given, add this path to the system path
    location = os.path.dirname(options.module)
    if location != "":
        sys.path.append(location)
   
    # Establish the module name, accomodating cases where the
    # .py extension has been included in the module name
    module_name = os.path.basename(options.module)
    if module_name.endswith(".py"):
        module_base_name = module_name[:-3]
    else:
        module_base_name = module_name

    # Import the specified module and map the specified fuction
    E.info( "importing module '%s' " % module_base_name)
    E.debug( "sys.path is: %s" % sys.path )

    module = importlib.import_module( module_base_name )
    try:
        function = getattr(module, options.function)
    except AttributeError as msg:
        raise AttributeError( msg.message + "unknown function, available functions are: %s" % \
                                  ",".join( [ x for x in dir(module) if not x.startswith("_") ] ) )

    if options.input_filenames:
        infiles = options.input_filenames
    else:
        infiles = False

    if options.output_filenames:
        outfiles = options.output_filenames
    else:
        outfiles = False

    # Parse the parameters into an array
    if options.params:
        params = [param.strip() for param in options.params.split(",")]
    else:
        params = False
              
    # deal with single file case
    if len(infiles) == 1: infiles = infiles[0]
    if len(outfiles) == 1: outfiles = outfiles[0]

    # Make the function call
    if infiles and outfiles and params:
        function(infiles, outfiles, params)
    elif infiles and outfiles and not params:
        function(infiles, outfiles)
    elif params:
        function(params)
    else:
        raise ValueError("Expecting infile+outfile+params or infile+outfile or params")

if __name__ == "__main__":
    sys.exit(main(sys.argv))
