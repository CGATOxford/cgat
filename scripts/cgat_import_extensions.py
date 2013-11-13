'''
cgat_import_extensions.py - test importing all scripts
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script attempts to import all the python scripts in the current
directory.

Importing a script is a pre-requisite for building documentation with
sphinx. A script/module that can not be imported will fail within sphinx.

Usage
-----

Example::

   python cgat_import_extensions.py

Type::

   python cgat_import_extensions.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import glob
import traceback
import imp

import CGAT.Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.set_defaults(
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    files = glob.glob( os.path.join( os.path.dirname(__file__), "*.py" ))

    files.sort()

    ## do sth
    c = E.Counter()

    for f in files:
        # if f != "pipeline_ancestral_repeats.py" : continue
        E.debug( "importing %s" % f )
        c.input += 1
        prefix, suffix = os.path.splitext( f )
            
        dirname, basename = os.path.split( prefix )

        if os.path.exists( prefix + ".pyc"):
            os.remove( prefix + ".pyc" )

        success = False
        try:
            __import__(basename, globals(), locals() )
            c.success += 1
            success = True
            options.stdout.write("PASS %s\n" % basename )
            options.stdout.flush()
        except ImportError, msg:
            c.import_fail += 1
            options.stdout.write( "FAIL %s\n%s\n" % (basename, msg) )
            options.stdout.flush()
            traceback.print_exc()
        except Exception, msg:
            c.other_fail += 1
            options.stdout.write( "FAIL %s\n%s\n" % (basename, msg))
            options.stdout.flush()
            traceback.print_exc()

        if success:
            # check for main
            path = os.path.abspath( os.path.dirname(__file__) )
            (file, pathname, description ) = imp.find_module( basename, [path,] )
            module = imp.load_module( basename, file, pathname, description)
            if "main" in dir(module):
                c.has_main +=1 
            else:
                options.stdout.write( "FAIL %s - %s\n" % (basename, "no main"))
                options.stdout.flush()
                c.no_main += 1

        c.output += 1

    E.info( c )
    
    ## write footer and output benchmark information.
    E.Stop()
    
if __name__ == "__main__":
    sys.exit( main( sys.argv) )
