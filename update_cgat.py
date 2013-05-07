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
update_cgat.py - change CGAT scripts to new package layout
==========================================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

scripts to update CGAT files such that they conform to
new package layout.

The script works by building a collection of all CGAT modules
and will replace import statements in all python scripts in
a directory. For example::

   import IOTools
   import Experiment as E
   import GTF, GFF

will become::

   import CGAT.IOTools as IOTools
   import CGAT.Experiment as E
   import CGAT.GTF as GTF
   import CGAT.GFF as GFF

The script will create ".bak" files for each python script
it processes (even if they are not changed).

The script is not recursive, if you have python scripts in subdirectories,
please execute the script on each directory individually. To automate this,
you can use the ``find`` command.

Usage
-----

To update all python scripts in the current directory, type::

   python update_cgat.py -d .

Type::

   python update_cgat.py --help

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
import shutil

import CGAT.Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-d", "--directory", dest="directory", type="string",
                      help="supply help"  )

    parser.set_defaults( directory = ".",
                         # using Andreas' repository in order to delay 
                         # changes to main repository in /ifs/devel/cgat
                         basename = "/ifs/devel/andreas/cgat/",                         
                         )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    # collect a list of python modules
    module_names = [ os.path.basename(x)[:-3] for x in glob.glob(options.basename + "CGAT/*.py" )]
    pipeline_names = [ os.path.basename(x)[:-3] for x in glob.glob(options.basename + "CGATPipelines/Pipeline*.py" )]

    if options.directory == "CGAT":
        files_to_update = glob.glob( os.path.join( basename, "scripts/*.py" )) +\
            glob.glob( os.path.join( options.basename, "scripts/*.pyx" )) +\
            glob.glob( os.path.join( options.basename, "CGATPipelines/pipeline_*.py" ) )+\
            glob.glob( os.path.join( options.basename, "CGATPipelines/Pipeline*.py" ) )
    else:
        files_to_update = glob.glob( os.path.join( options.directory, "*.py" ))

    E.info("updating %i python scripts/modules" % len(files_to_update))

    counter = E.Counter()

    for script in files_to_update:

        counter.input += 1
        print "working on", script

        inf = open( script )
        lines = inf.readlines()
        inf.close()

        # create a backup copy
        shutil.move( script, script + ".bak" )
        
        outf = open( script, "w" )
        updated = False
        for line in lines:
            if re.match( "import ", line ):
                if " as " in line:
                    try:
                        module, name = re.match( "import (\S+) as (\S+)\s*$", line).groups()
                    except AttributeError as msg:
                        raise AttributeError( "parsing error in line '%s': '%s'" % (line[:-1], msg))
                    if module in module_names:
                        line = "import CGAT.%s as %s\n" % (module, name )
                        updated = True
                else:
                    try:
                        modules = re.match( "import (.+)", line ).groups()[0]
                    except AttributeError as msg:
                        raise AttributeError( "parsing error in line '%s': '%s'" % (line[:-1], msg))
                    modules = [ x.strip() for x in modules.split(",")]
                    for module in modules:
                        if module in module_names:
                            outf.write( "import CGAT.%s as %s\n" % (module, module ))
                            updated = True
                        elif module in pipeline_names:
                            outf.write( "import CGATPipelines.%s as %s\n" % (module, module ))
                            updated = True
                        else:
                            outf.write( "import %s\n" % module )
                    continue

            outf.write( line )
        outf.close()
        
        if updated: counter.updated += 1

    E.info( "summary: %s" % str(counter))

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )





