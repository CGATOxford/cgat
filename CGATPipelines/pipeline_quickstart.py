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
import shutil

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

    destination_dir = options.destination
    name = options.name

    # create directories
    for d in ("", "src", "report",
              "src/pipeline_docs", 
              "src/pipeline_%s" % options.name,
              reportdir,
              "%s/_templates" % reportdir,
              "%s/pipeline" % reportdir,
              "%s/trackers" % reportdir ):

        dd = os.path.join( destination_dir, d )
        if not os.path.exists( dd ): os.makedirs( dd )
        
    # copy files
    # replaces all instances of template with options.name within
    # filenames and inside files.
    rx_file = re.compile( "template" )
    rx_template = re.compile( "@template@" )
    rx_reportdir = re.compile( "@reportdir@" )

    srcdir = os.path.dirname( __file__ )

    def copy( src, dst ):

        fn_dest = os.path.join( destination_dir, dst, rx_file.sub(name, src) )
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

    def copytree( src, dst ):

        fn_dest = os.path.join( destination_dir, dst, rx_file.sub(name, src) )
        fn_src = os.path.join( srcdir, "pipeline_template_data", src)

        print fn_dest, fn_src
        if os.path.exists( fn_dest ) and not options.force:
            raise OSError( "file %s already exists - not overwriting." % fn_dest )

        shutil.copytree( fn_src, fn_dest )

    for f in ( "conf.py",
               "pipeline.ini" ):
        copy( f, 'src/pipeline_%s' % options.name )
        
    for f in ( "pipeline_template.py", ):
        copy( f, 'src' )

    # create links
    for src, dest in ( ("conf.py", "conf.py"), 
                       ( "pipeline.ini", "pipeline.ini" )):
        d = os.path.join( "report", dest )
        if os.path.exists( d ) and options.force: 
            os.unlink( d )
        os.symlink( os.path.join( confdir, src), d )

    for f in ( "cgat_logo.png",
               "index.html" ):
        copy( f, "%s/_templates" % reportdir )

    for f in ("themes",):
        copytree( f, "src/pipeline_docs" )
        
    for f in ( "contents.rst",
               "pipeline.rst",
               "__init__.py" ):
        copy( f, reportdir )

    for f in ( "Dummy.rst",
               "Methods.rst"):
        copy( f, "%s/pipeline" % reportdir )

    for f in ( "TemplateReport.py", ):
        copy( f, "%s/trackers" % reportdir )

    absdest = os.path.abspath( destination_dir )

    print """
Welcome to your new %(name)s CGAT pipeline.

All files have been successfully copied to `%(destination_dir)s`. In order to start
the pipeline, go to `%(destination_dir)s/report`

   cd %(destination_dir)s/report

You can start the pipeline by typing:

   python ../src/pipeline_%(name)s.py -v 5 -p 5 make full

To build the report, type:

   python ../src/pipeline_%(name)s.py -v 5 -p 5 make build_report

The report will be in file:/%(absdest)s/report/report/html/index.html.

The source code for the pipeline is in %(destination_dir)s/src.

""" % locals()

   
if __name__ == "__main__":
    sys.exit(main())
