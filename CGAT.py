################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
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
CGAT.py - CGAT project specific functions
=========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

The :mod:`CGAT` module contains various utility functions
for working on CGAT projects.

API
----
'''

from ruffus import *

PROJECT_ROOT = '/ifs/projects'

def getProjectId():
    '''cgat specific method: get the (obfuscated) project id
    based on the current working directory.
    '''
    curdir = os.path.abspath(os.getcwd())
    if not curdir.startswith( PROJECT_ROOT ):
        raise ValueError( "method getProjectId no called within %s" % PROJECT_ROOT )
    prefixes = len(PROJECT_ROOT.split("/"))
    # patch for projects in mini directory
    if "/mini/" in curdir: prefixes += 1
    rootdir = "/" + os.path.join( *(curdir.split( "/" )[:prefixes+1]) )
    f = os.path.join( rootdir, "sftp", "web" )
    if not os.path.exists(f):
        raise OSError( "web directory at '%s' does not exist" % f )
    target = os.readlink( f )
    return os.path.basename( target )


def getProjectName():
    '''cgat specific method: get the name of the project 
    based on the current working directory.'''

    curdir = os.path.abspath(os.getcwd())
    if not curdir.startswith( PROJECT_ROOT ):
        raise ValueError( "method getProjectName no called within %s" % PROJECT_ROOT )
    prefixes = len(PROJECT_ROOT.split("/"))
    return curdir.split( "/" )[prefixes]

def publish_report( prefix = "", 
                    patterns = [], 
                    project_id = None,
                    prefix_project = "/ifs/projects",
                    ):
    '''publish report into web directory.

    Links export directory into web directory.

    Copies html pages and fudges links to the pages in the
    export directory.

    If *prefix* is given, the directories will start with prefix.

    *patterns* is an optional list of two-element tuples (<pattern>, replacement_string).
    Each substitutions will be applied on each file ending in .html.

    If *project_id* is not given, it will be looked up. This requires
    that this method is called within a subdirectory of PROJECT_ROOT.

    .. note::
       This function is CGAT specific.

    '''

    if not prefix:
        prefix = PARAMS.get( "report_prefix", "" )

    web_dir = PARAMS["web_dir"]
    if project_id == None: 
        project_id = getProjectId()

    src_export = os.path.abspath( "export" )
    dest_report = prefix + "report"
    dest_export = prefix + "export"

    curdir = os.path.abspath( os.getcwd() )

    def _link( src, dest ):
        '''create links. 

        Only link to existing targets.
        '''
        dest = os.path.abspath( os.path.join( PARAMS["web_dir"], dest ) )
        if os.path.exists( dest ):
            os.remove(dest)
        
        if os.path.exists( src ):
            #IMS: check if base path of dest exists. This allows for prefix to be a 
            #nested path structure e.g. project_id/
            if not os.path.exists(os.path.dirname(os.path.abspath(dest))):
                os.mkdir(os.path.dirname(os.path.abspath(dest)))

            os.symlink( os.path.abspath(src), dest )

    def _copy( src, dest ):
        dest = os.path.abspath( os.path.join( PARAMS["web_dir"], dest ) )
        if os.path.exists( dest ): shutil.rmtree( dest )
        shutil.copytree( os.path.abspath(src), dest ) 

    # publish export dir via symlinking
    E.info( "linking export directory in %s" % dest_export )
    _link( src_export, dest_export )

    # publish web pages by copying
    E.info( "publishing web pages in %s" % os.path.abspath( os.path.join( PARAMS["web_dir"], dest_report )))
    _copy( os.path.abspath("report/html"), dest_report )

    # substitute links to export and report
    _patterns = [ (re.compile( src_export ), 
                   "http://www.cgat.org/downloads/%(project_id)s/%(dest_export)s" % locals() ), 
                  (re.compile( '(%s)/report' % os.path.join( prefix_project, getProjectName() ) ),
                   "http://www.cgat.org/downloads/%(project_id)s/%(dest_report)s" % locals() ),
                  (re.compile( '(%s)/_static' % curdir),
                   "http://www.cgat.org/downloads/%(project_id)s/%(dest_report)s/_static" % locals() )]
    
    _patterns.extend( patterns )

    for root, dirs, files in os.walk(os.path.join( web_dir, dest_report)):
        for f in files:
            fn = os.path.join( root, f )
            if fn.endswith(".html" ):
                with open( fn ) as inf:
                    data = inf.read()
                for rx, repl in _patterns:
                    data = rx.sub( repl, data )
                outf = open( fn, "w" )
                outf.write( data )
                outf.close()

    E.info( "report has been published at http://www.cgat.org/downloads/%(project_id)s/%(dest_report)s" % locals())
