#!/bin/env python
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
pipeline_quickstart.py - setup a new pipeline
=============================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Usage
-----

Example::

   python pipeline_quickstart.py --name=chipseq

Type::

   python pipeline_quickstart.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import re
import os
import optparse

def main( argv = sys.argv ):

    parser = optparse.OptionParser( version = "%prog version: $Id$", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-d", "--dest", dest="destination", type="string",
                      help="destination directory." )

    parser.add_option("-n", "--name", dest="name", type="string",
                      help="name of this pipeline. 'pipeline_' will be prefixed." )

    parser.add_option("-f", "--force", dest="force", action="store_true",
                      help="overwrite existing files." )
    
    parser.set_defaults(
        destination = ".",
        name = None,
        force = False,
        )
    
    (options, args) = parser.parse_args()

    if not options.name: raise ValueError( "please provide a pipeline name" )

    reportdir = os.path.abspath( "src/pipeline_docs/pipeline_%s" % options.name )
    confdir = os.path.abspath( "src/pipeline_%s" % (options.name))

    dest = options.destination
    name = options.name

    # create directories
    for d in ("", "src", "report",
              "src/pipeline_docs", 
              "src/pipeline_%s" % options.name,
              reportdir,
              "%s/_templates" % reportdir,
              "%s/pipeline" % reportdir,
              "%s/trackers" % reportdir ):

        dd = os.path.join( dest, d )
        if not os.path.exists( dd ): os.makedirs( dd )
        
    # copy files
    # replaces all instances of template with options.name within
    # filenames and inside files.
    rx_file = re.compile( "template" )
    rx_template = re.compile( "@template@" )
    rx_reportdir = re.compile( "@reportdir@" )

    srcdir = os.path.dirname( __file__ )

    def copy( src, dst ):
        fn_dest = os.path.join( dest, dst, rx_file.sub(name, src) )
                                    
        fn_src = os.path.join( srcdir, "pipeline_template_data", src)

        if os.path.exists( fn_dest ) and not options.force:
            raise OSError( "file %s already exists - not overwriting." % fn_dest )

        outfile = open(fn_dest, "w")
        infile = open(fn_src )
        for line in infile:
            outfile.write( rx_reportdir.sub( reportdir, 
                                             rx_template.sub( name, line ) ))
            
        outfile.close()
        infile.close()

    for f in ( "sphinxreport.ini",
               "conf.py",
               "pipeline.ini" ):
        copy( f, 'src/pipeline_%s' % options.name )
        
    for f in ( "pipeline_template.py", ):
        copy( f, 'src' )

    # create links
    for src,dest in ( ("sphinxreport.ini", "sphinxreport.ini"),
                      ("conf.py", "conf.py"),
                      ( "pipeline.ini", "pipeline.ini" )):
        os.symlink( os.path.join( confdir, src), os.path.join( "report", dest ) )

    for f in ( "cgat_logo.png",
               "index.html",
               "gallery.html" ):
        copy( f, "%s/_templates" % reportdir )

    for f in ( "contents.rst",
               "pipeline.rst",
               "analysis.rst",
               "__init__.py" ):
        copy( f, reportdir )

    for f in ( "Dummy.rst",
               "Methods.rst"):
        copy( f, "%s/pipeline" % reportdir )

    for f in ( "TemplateReport.py", ):
        copy( f, "%s/trackers" % reportdir )

    absdest = os.path.abspath( dest )

    print """
Welcome to your new %(name)s CGAT pipeline.

All files have been successfully copied to `%(dest)s`. In order to start
the pipeline, go to `%(dest)s/report`

   cd %(dest)s/report

You can start the pipeline by typing:

   python ../src/pipeline_%(name)s.py -v 5 -p 5 make full

To build the report, type:

   python ../src/pipeline_%(name)s.py -v 5 -p 5 make build_report

The report will be in file:/%(absdest)s/report/report/html/index.html.

The source code for the pipeline is in %(dest)s/src.

""" % locals()

   
if __name__ == "__main__":
    sys.exit(main())
