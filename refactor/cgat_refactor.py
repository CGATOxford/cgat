################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
cgat_refactor.py - refactor CGAT Code
=====================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Usage
-----

Example::

   python cgat_refactor.py --rename=rename.txt

Type::

   python cgat_refactor.py --help

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
import CGAT.IOTools as IOTools

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-r", "--rename", dest="rename", type="string",
                      help="rename scripts"  )

    parser.add_option( "--split-prefix", dest="split_prefix", type="string",
                      help="move scripts with prefix to subdirectory"  )

    parser.add_option("-n", "--dry-run", dest="dry_run", action = "store_true",
                      help="dry run, do not implement any changes"  )


    parser.set_defaults(
        scriptsdir = "scripts", 
        dirs = ["CGAT", "CGATPipelines", "scripts", "makefiles" ],
        dry_run = False )
    

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    scriptsdir = options.scriptsdir

    counter = E.Counter()

    map_old2new = {}

    if options.rename:
        with IOTools.openFile( options.rename, "r") as inf:
            for line in inf: 
                if line.startswith("#"): continue
                if line.startswith("old"): continue
                try:
                    old, new = line[:-1].split("\t")
                except ValueError:
                    continue
                if not os.path.exists( os.path.join( scriptsdir, old )):
                    E.warn( "%s does not exist - no renaming" % old )
                    continue
                map_old2new[old] = new


    elif options.split_prefix:

        if not os.path.exists( os.path.join( scriptsdir, options.split_prefix )):
            E.warn( "destination %s does not exist - no renaming" % options.split_prefix )
            return

        scripts = glob.glob( "%s/%s_*.py" % (scriptsdir, options.split_prefix ))
        if len(scripts) == 0:
            E.info("nothing to change")
            return
        
        for script in scripts:
            scriptname = os.path.basename( script )
            newname = scriptname[len(options.split_prefix)+1:]
            map_old2new[ scriptname ] = "%s/%s" % (options.split_prefix, newname )

    if len(map_old2new) == 0: 
        E.info("nothing to change")
        return

    for old, new in map_old2new.items():
        statement = "hg mv %(scriptsdir)s/%(old)s %(scriptsdir)s/%(new)s" % locals()
        counter.renamed += 1
        
        if options.dry_run:
            E.info( statement )
        else:
            E.run( statement )

    for d in options.dirs:
        for root, dirs, files in os.walk(d):
            for f in files:
                if f.endswith(".pyc"): continue
                fn = os.path.join( root, f )
                with IOTools.openFile( fn, "r") as inf:
                    old_data = inf.read()

                changed = False
                for old_name, new_name in map_old2new.items():
                    new_data = re.sub( old_name, new_name, old_data )
                    if old_data != new_data:
                        changed = True
                        E.info( "changed: %s : %s to %s" % (fn, old_name, new_name))
                    old_data = new_data

                if changed:
                    counter.changed += 1

                    if not options.dry_run:
                        with IOTools.openFile( fn, "w" ) as outf:
                            outf.write( new_data )

    E.info( str(counter) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
