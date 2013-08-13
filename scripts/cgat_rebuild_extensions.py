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
cgat_rebuild_extensions.py - rebuild all cython extensions
==========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script rebuilds all cython extensions in the source directory.

Some scripts in the repository make use of ``pyximport`` to compile
associated cython scripts with embedded C code. Theses scripts are
automatically re-compiled if the script has changed, but this process
can fail if:

   * the script is executed on a machine without a C-compiler
   * some underlying libraries have changed.

Thus, it is safer to rebuild all scripts on a machine with a C compiler
before running a script in production on a cluster, where not all nodes
might be fully configured for compilation.

Usage
-----

Example::

   python cgat_rebuild_extensions.py

Type::

   python cgat_rebuild_extensions.py --help

for command line help.

Documentation
-------------

Code
----

'''

import os
import sys
import re
import optparse
import glob

import CGAT.Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-i", "--test-option", dest="test_option", type="string",
                      help="test option [default=%default]."  )

    parser.set_defaults(
        test_option = "test"
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    files = glob.glob( os.path.join( os.path.dirname(__file__), "*.pyx" ))

    ## do sth
    ninput, nskipped, noutput = 0, 0, 0

    for f in files:
        E.info( "rebuilding %s" % f )
        ninput += 1
        prefix, suffix = os.path.splitext( f )
        for ext in (".c", ".pyxbldc" ):
            try: os.remove( prefix + ext )
            except OSError: pass
            
        dirname, basename = os.path.split( prefix )
        assert basename.startswith("_")
        
        scriptname = os.path.join( dirname, basename[1:] ) + ".py" 
        if not os.path.exists( scriptname ):
            E.warn("script %s does not exist - skipped" % scriptname)
            nskipped += 1
            continue
        
        E.info( "compiling %s" % scriptname )
        os.system( "%s %s --help > /dev/null" % (sys.executable, scriptname ) )

        noutput += 1

    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput,nskipped) )
    
    ## write footer and output benchmark information.
    E.Stop()
    
if __name__ == "__main__":
    sys.exit( main( sys.argv) )
