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
Pipeline.py - Tools for ruffus pipelines
========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

The :mod:`Pipeline` module contains utility functions for ruffus pipelines.

The module brings together several functionalities into a single module:

* `Controlling ruffus from the command line`_
* `Parsing configuration values`_
* `Environment functions`_
* `Database upload`_
* `Running remote jobs`_

:py:mod:`Pipeline` configures ruffus_ to use threads instead of multiple
processes. The underlying motivation is that the pipeline only encodes the logic of successive
steps, while computation is farmed-out to a cluster. 

Typically, the pipeline would be started on an interactive submit host and run until all jobs have finished. 
Multiple pipelines might be running at the same time started by different users.
Activity on the submit host should thus be light-weight and not occupy multiple CPU. Python's threading 
implementation uses the `GIL <http://en.wikipedia.org/wiki/Global_Interpreter_Lock>`_, which naturally 
limits the resources the pipeline uses to a single CPU (for python code).

The basic usage for Pipeline.py is the following::

   import CGAT.Pipeline as P

   ## Read pipeline parameters
   PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini" ] )

   ## Your pipeline code here
   ##

   ## start the pipeline
   if __name__== "__main__":
       sys.exit( P.main(sys.argv) )

Controlling ruffus from the command line
-----------------------------------------

Pipeline.py implements a :py:func:`main` function. This function
parses options from the command line and depending on the command given,
performs a variety of actions:

make <task>
   run all tasks required to build :term:`task`

show <task>
   show tasks required to build :term:`task` without executing them

plot <task>
   plot image (requires `inkscape <http://inkscape.org/>`_) of pipeline state for :term:`task`

touch <task>
   touch files without running :term:`task` or its pre-requisites. This sets the 
   timestamps for files in :term:`task` and its pre-requisites such that they will 
   seem up-to-date to the pipeline.

config
   write a new configuration file :file:`pipeline.ini` with default values. An existing 
   configuration file will not be overwritten (:py:func:`writeConfigFiles`).

clone <srcdir>
   clone a pipeline from :file:`srcdir` into the current
   directory. Cloning attempts to conserve disk space by linking. (:py:func:`clonePipeline`).

dump
   output configuration options and exit.

The :py:func:`main` function also sets up logging for the pipeline.
Logging information is output to stdout, but also written to the 
file :file:`pipeline.log` in the current directory. The logging format
follows the :mod:`Experiment` module.

.. autofunction:: Pipeline.main

Parsing configuration values
-----------------------------

:py:mod:`Pipeline` provides a mechanism for working with configuration values.
Pipeline.py reads configuration values in the ConfigParser format and provides them
in a global dictionary called :py:data:`PARAMS`. 

The central function is :py:func:`getParameters`, which will read configuration
values from one or more configuration files. Values are automatically converted.

.. autofunction:: Pipeline.getParameters

The following functions convert a parameter string of comma separated values into a list or tuple:

.. autofunction:: Pipeline.asList

.. autofunction:: Pipeline.asTuple

Database upload
---------------

The :py:mod:`Pipeline` module provides several utility functions for uploading
:term:`tsv` formatted files into a database. These methods are designed
to work within ruffus_ tasks.

The basic method is the :py:func:`load` method. Variants such as 
:py:func:`mergeAndLoad` and :py:func:`concatenateAndLoad` will combine
multiple input files before loading.

Uploading is implemented via the :doc:`../scripts/csv2db` script.

.. autofunction:: Pipeline.load

.. autofunction:: Pipeline.mergeAndLoad

.. autofunction:: Pipeline.concatenateAndLoad

The module also provides some convenience functions for working with tables.

.. autofunction:: Pipeline.quote

.. autofunction:: Pipeline.toTable

Environment functions
--------------------------

The :py:mod:`Pipeline` module provides utility functions to work with
temporary files and to check for pipeline requirements.

.. autofunction:: Pipeline.getTempFilename

.. autofunction:: Pipeline.getTempFile

.. autofunction:: Pipeline.getTempDir

The :py:mod:`Pipeline` module provides functions to check for executables
and scripts in the :envvar:`PATH`.

.. autofunction:: Pipeline.which

.. autofunction:: Pipeline.checkExecutables

.. autofunction:: Pipeline.checkScripts

.. autofunction:: Pipeline.snip

Running remote jobs
---------------------------

:py:mod:`Pipeline` provides a mechanism to run jobs remotely within
pipelines. The :py:func:`run` method runs command line statements
remotely, while the :py:func:`submit` method runs python code
remotely.

.. autofunction:: Pipeline.run

The concept of :py:func:`.run` is to minimize the coding necessary
to run remote jobs within a ruffus task. Typically, not much more than
the statement and a call to :py:func:`run` are required. For example::

   @transform( "*.unsorted.gz", suffix(".unsorted.gz"), ".sorted.gz")
   def sortAToB( infile, outfile ):
       """sort input."""

       statement = """
       gunzip < %(infile)s | sort -k %(sort_key)s | gzip > %(outfile)s
       """

       P.run()

Here, the variable ``infile`` and ``outfile`` will be interpolated in the
context of ``sortAToB``, while ``sort_key`` will be set by a configuration
variable in :py:data:`PARAMS`.

See :ref:`PipelineCommands` for more information.

.. autofunction:: Pipeline.submit

The purpose of :py:func:`submit` is equivalent to :py:func:`run`. Instead
of running a command line statement, a function within a script can be specified
that will be run. 

Note that filenames only are passed through :py:func:`submit`, no variable state is
transmitted. For example::

   @transform( "*.a.gz", suffix(".a.gz"), ".b.gz")
   def countAToB( infile, outfile ):
       """sort input."""

       P.submit( "MyPipelineFunctions.py", "sortFile", 
                  infiles = infile,
                  outfiles = outfile )

The code above will call the method sortFile in MyPipelineFunctions.py. Note that
MyPipelineFunctions.py need to be importable by being in the current directory,
specified with a relative or absolute path or within :envvar:`PYTHONPATH`.

The call signature of ``MyPipelineFunction.sortFile`` needs to correspond to
a typical ruffus_ task ``(infile, outfile, ...)``, for example::

   def sortFile(infile, outfile ):
       """sort a file"""
       with open( infile, "r") as inf:
           lines = inf.readlines()
       with open( outfile, "w") as outf:
           outf.write( "".join( sorted( lines )))

For some tasks, which do not provide an output, it is sometimes useful to provide
an empty sentinel file. 

.. autofunction:: Pipeline.touch

.. _ruffus: http://www.ruffus.org.uk/
.. _sqlite: http://www.sqlite.org/
.. _sphinxreport: http://code.google.com/p/sphinx-report/

Complete reference
------------------

'''
import os, sys, re, subprocess, optparse, stat, tempfile
import time, random, inspect, types, multiprocessing
import logging, collections, shutil, glob, gzip
import ConfigParser
import Database

# talking to a cluster
import drmaa
# talking to mercurial
import hgapi

# CGAT specific options - later to be removed
from CGAT import *

from ruffus import *

# use threading instead of multiprocessing
from multiprocessing.pool import ThreadPool
# note that threading can cause problems with rpy.
task.Pool = ThreadPool

import logging as L
import Experiment as E
import IOTools

# global options and arguments
GLOBAL_OPTIONS, GLOBAL_ARGS = None, None

# global drmaa session
GLOBAL_SESSION = None

# local temporary directory to use
TMPDIR = '/scratch'

hostname = os.uname()[0]

class PipelineError( Exception ): pass

# possible to use defaultdict, but then statements will
# fail on execution if a parameter does not exists, and not
# while building the statement. Hence, use dict.
LIB_DIR=os.path.dirname( __file__ )
ROOT_DIR=os.path.dirname( LIB_DIR )
SCRIPTS_DIR=os.path.join( ROOT_DIR, "scripts" )

####################################################################
####################################################################
####################################################################
## Configuration variables
####################################################################

## global variable containing configuration values
PARAMS = { 
    'scriptsdir' : SCRIPTS_DIR,
    'toolsdir' : SCRIPTS_DIR,
    'cmd-farm' : """%s/farm.py 
                --method=drmaa 
                --cluster-priority=-10 
		--cluster-queue=all.q 
		--cluster-num-jobs=100 
                --bashrc=%s/bashrc.cgat
		--cluster-options="" """ % (SCRIPTS_DIR,SCRIPTS_DIR),
    'cmd-sql' : """sqlite3 -header -csv -separator $'\\t' """,
    'cmd-run' : """%s/run.py""" % SCRIPTS_DIR
    }

# patch until parameter sharing is resolved between CGAT module
# and the pipelines module.
import CGAT
CGAT.PARAMS = PARAMS

CONFIG = {}

def configToDictionary( config ):
    '''convert a :py:class:`ConfigParser.ConfigParser` object 
    into a dictionary.

    :param config: ConfigParser object
    :returns: dictionary of configuration values.

    Sections and keys are combined with an underscore. If
    a key without section does not exist, it will be added 
    plain.
    
    The section [DEFAULT] is equivalent to [general].
    '''
    
    p = {}
    for section in config.sections():
        for key,value in config.items( section ):
            v = IOTools.convertValue( value )
            p["%s_%s" % (section,key)] = v
            if section in ( "general", "DEFAULT" ):
                p["%s" % (key)] = v
               
    for key, value in config.defaults().iteritems():
        p["%s" % (key)] =  IOTools.convertValue( value )
        
    return p

def getParameters( filenames = ["pipeline.ini",],
                   defaults = None ):
    '''read a config file and return as a dictionary.

    :param filenames: list of filenames to read
    :param defaults: dictionary with default values.
    :returns: dictionary with parameters. 

    Sections and keys are combined with an underscore. If
    a key without section does not exist, it will be added 
    plain.

    For example::

       [general]
       input=input1.file

       [special]
       input=input2.file

    will be entered as ``{ 'general_input' : "input1.file",
    'input: "input1.file", 'special_input' : "input2.file" }``

    This function also updates the module-wide parameter map.
    
    Configuration values are automatically converted into
    integers and float values.

    The section [DEFAULT] is equivalent to [general].
    '''

    global CONFIG

    CONFIG = ConfigParser.ConfigParser()
    CONFIG.read( filenames )

    p = configToDictionary( CONFIG )
    if defaults: PARAMS.update( defaults )
    PARAMS.update( p )

    return PARAMS

def loadParameters( filenames  ):
    '''load parameters from a tuple of filenames.

    :param filenames: list of filenames
    :returns: dictionary with parameters. 

    Parameters are processed in the same way as
    :py:func:`getParameters`, but the global parameter 
    dictionary is not updated.
    '''
    config = ConfigParser.ConfigParser()
    config.read( filenames )

    p = configToDictionary( config )
    return p

def substituteParameters( **kwargs ):
    '''combine configuration values with local and task specific parameters

    :param kwargs: Entries in ``**kwargs`` substitute default
                   values in :py:data:`PARAMS`.
    :returns: Dictionary with configuration values.

    Task specific configuration values take highest priority. Task
    specific parameters start with the contents of the ``outfile``
    parameter. 

    For example, if :py:data:`PARAMS` contains the following::

       PARAMS = { 'outfile' : 'test.gtf', 
                  'option1' : 12,
                  'outfile_option1' : 10 }

    this method will return::

       PARAMS = { 'outfile' : 'test.gtf', 'option1' : 10 }

    '''

    # build parameter dictionary
    # note the order of addition to make sure that kwargs takes precedence
    local_params = dict(PARAMS.items() + kwargs.items())

    if "outfile" in local_params:
        # replace specific parameters with task (outfile) specific parameters
        outfile = local_params["outfile"]
        for k in local_params.keys():
            if k.startswith(outfile):
                p = k[len(outfile)+1:]
                if p not in local_params:
                    raise KeyError( "task specific parameter '%s' does not exist for '%s' " % (p,k))
                E.debug( "substituting task specific parameter for %s: %s = %s" % (outfile,p,local_params[k] ) )
                local_params[p] = local_params[k]

    return local_params

def checkParameter( key ):
    if key not in PARAMS:
        raise ValueError("need `%s` to be set" % key )

def asList( param ):
    '''return a param as a list'''
    if type(param) not in (types.ListType, types.TupleType):
        try:
            params = [x.strip() for x in param.strip().split(",")]
        except AttributeError:
            params = [param.strip()]
        return [ x for x in params if x != ""]
    else: return param

def asTuple( param ):
    '''return a param as a list'''
    return tuple(asList( param ))

def flatten(l, ltypes=(list, tuple)):
    '''flatten a nested list/tuple.'''
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)

def asDict( param ):
    '''return a section of configuration file as a dictionary.'''
    return dict(CONFIG.items(param))

##############################################################################################
##############################################################################################
##############################################################################################
def checkFiles( filenames ):
    """check for the presence/absence of files"""
    
    missing = []
    for filename in filenames:
        if not os.path.exists( filename ):
            missing.append( filename )

    if missing:
        raise ValueError( "missing files: %s" % ",".join(missing) )

##############################################################################################
##############################################################################################
##############################################################################################
def clone( infile, outfile ):
    '''create a clone of ``infile`` named ``outfile``
    by creating a soft-link.
    '''
    try:
        os.symlink( infile, outfile )
    except OSError:
        pass
    
##############################################################################################
##############################################################################################
##############################################################################################
def touch( filename, times = None ):
    '''update/create a sentinel file.

    Compressed files (ending in .gz) are created
    as empty 'gzip' files, i.e., with a header.
    '''
    existed = os.path.exists(filename)
    fhandle = file(filename, 'a')

    if filename.endswith (".gz") and not existed:
        # this will automatically add a gzip header
        fhandle = gzip.GzipFile( filename, fileobj = fhandle )
        
    try:
        os.utime(filename, times)
    finally:
        fhandle.close()

##############################################################################################
##############################################################################################
##############################################################################################
def isEmpty( filename ):
    '''return true if file *filename* is empty.'''
    return os.stat(filename)[6]==0

##############################################################################################
##############################################################################################
##############################################################################################
def snip( filename, extension = None, alt_extension = None, path = False):
    '''return filename without extension.

    :param filename: filename to be truncated
    :param extension: extension to take off. If given
           this method ensures that filename has the extension.
    :param alt_extension: alternative extension. If extension
           not present, try alternative extension. 
    :param path: If path is set to True, the path is stripped from the file name.
    :returns: the filename without the extension.
    '''

    if extension: 
        if filename.endswith( extension ):
            root = filename[:-len(extension)]
        elif alt_extension and filename.endswith( alt_extension ):
            root = filename[:-len(alt_extension)]
        else:
            raise ValueError("'%s' expected to end in '%s'" % (filename, extension))
    else:
        root, ext = os.path.splitext( filename )

    if path==True: snipped = os.path.basename(root)
    else: snipped = root

    return snipped
    
##############################################################################################
##############################################################################################
## Temporary files and directories
##############################################################################################
def getTempFilename( dir = TMPDIR ):
    '''return temporary filename.

    :param dir: directory in which to create temporary file.
    :returns: a filename of a temporary file.
    
    The temporary file is persistent. The caller needs to delete 
    the temporary file.
    '''
    tmpfile = tempfile.NamedTemporaryFile(dir=dir,delete=False, prefix = "ctmp" )
    tmpfile.close()
    return tmpfile.name

def getTempFile( dir = TMPDIR ):
    '''return a temporary file.

    :param dir: directory in which to create temporary file.
    :returns: a file object to a temporary file.
    
    The temporary file is persistent. The caller needs to delete 
    the temporary file.
    '''
    return tempfile.NamedTemporaryFile(dir=dir, delete=False, prefix = "ctmp" )

def getTempDir( dir = TMPDIR ):
    '''return a temporary directory.

    :param dir: directory in which to create temporary directory.
    :returns: name of the temporary directory.

    The temporary directory is persistent. The caller needs to delete 
    the temporary directory.
    '''
    return tempfile.mkdtemp( dir=dir, prefix = "ctmp" )

##############################################################################################
##############################################################################################
##############################################################################################
def checkExecutables( filenames ):
    """check for the presence/absence of executables in PATH.

    *param filenames* list of filenames to check check.
    *raise ValueError* if any excecutables are missing
    """
    
    missing = []

    for filename in filenames:
        if not which( filename ):
            missing.append( filename )

    if missing:
        raise KeyError( "missing executables: %s" % ",".join(missing) )

##############################################################################################
##############################################################################################
##############################################################################################
def which( filename ):
    '''return path to executable.

    *param filename* filename of executable.
    *returns* the path to the executable. Returns None if not found.

    This method looks up an executable in the :envvar:`PATH`
    path variable. The first location will be returned.
    '''

    if not os.environ.has_key('PATH') or os.environ['PATH'] == '':
        p = os.defpath
    else:
        p = os.environ['PATH']

    pathlist = p.split (os.pathsep)

    for path in pathlist:
        f = os.path.join(path, filename)
        if os.access(f, os.X_OK):
            return f
    return None


##############################################################################################
##############################################################################################
##############################################################################################
def checkScripts( filenames ):
    """check for the presence/absence of files.

    *param filenames*  list of filenames to check check. Filenames should include
                       any paths.
    *raise ValueError* if any excecutables are missing.
    
    """
    missing = []
    for filename in filenames:
        if not os.path.exists( filename ):
            missing.append( filename )

    if missing:
        raise ValueError( "missing scripts: %s" % ",".join(missing) )

##############################################################################################
##############################################################################################
##############################################################################################
## Logging functions
##############################################################################################
def log( loglevel, message ):
    """log message at loglevel."""
    E.log( loglevel, message )

def info( message ):
    L.info( message )

def warning( message ):
    E.warning( message )

def warn( message ):
    E.warning( message )

def debug( message ):
    L.debug( message )

def error( message ):
    E.error( message )
        
def critical( message):
    E.critical( message )

##############################################################################################
##############################################################################################
##############################################################################################
def quote( s ):
    '''quote string such that is applicable for a table name.

    :param s: string to be quoted
    :returns: quoted string
    '''
    return re.sub( "[-(),\[\].]", "_", s)

def toTable( outfile ):
    '''convert an outfile (filename) into
    a quoted table name.

    :param outfile: filename to be quoted
    :returns: quoted string.

    '''
    assert outfile.endswith( ".load" ) 
    name = os.path.basename( outfile[:-len(".load")] )
    return quote( name )

def load( infile, 
          outfile = None, 
          options = "", 
          transpose = None,
          tablename = None):
    '''import a :term:`tsv` formatted table into database.

    :param infile: data to be uploaded
    :param outfile: output file. The output file will contain
         the logging information of the loading. 
         Determines the table name if given and tablename not set.
         Needs to end in ".load".
    :param options: options for loading. See :py:mod:`CSV2DB`.
    :param transpose: transpose data before uploading.
    :param tablename: tablename to use. 

    The table name is given by outfile without the
    ".load" suffix.

    If *transpose* is set, the table will be transposed before loading.
    The first column in the first row will be set to the string
    within transpose.
    '''

    if not tablename:
        tablename = toTable( outfile )
    
    statement = []
    if infile.endswith(".gz"): statement.append( "zcat %(infile)s" )
    else: statement.append( "cat %(infile)s" )

    if transpose:
        statement.append( "python %(scriptsdir)s/table2table.py --transpose --set-transpose-field=%(transpose)s" )

    statement.append('''
    python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              %(options)s 
              --table=%(tablename)s 
    > %(outfile)s
    ''')

    statement = " | ".join( statement)

    run()

def concatenateAndLoad( infiles, outfile, regex_filename = None, header = None, cat = None, titles = False, options = "" ):
    '''concatenate categorical tables and load into a database.

    Concatenation assumes that the header is the same in all files.
    The first file will be taken in completion, headers
    in other files will be removed.
    '''
    
    infiles = " ".join(infiles)

    tablename = toTable( outfile )

    passed_options = options
    load_options,options = [], []

    if regex_filename:
        options.append( '--regex-filename="%s"' % regex_filename )

    if header:
        load_options.append( "--header=%s" % header )

    if not cat:
        cat = "track"
        
    if titles == False:
        no_titles = "--no-titles"
    else: no_titles = ""

    options = " ".join(options)
    load_options = " ".join(load_options) + " " + passed_options
    statement = '''python %(scriptsdir)s/combine_tables.py
                     --cat=%(cat)s
                     %(no_titles)s
                     %(options)s
                   %(infiles)s
                   | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
                      --index=track
                      --table=%(tablename)s 
                      %(load_options)s
                   > %(outfile)s'''
    run()

def mergeAndLoad( infiles, 
                  outfile, 
                  suffix = None, 
                  columns=(0,1), 
                  regex = None, 
                  row_wise = True ):
    '''merge categorical tables and load into a database.

    Columns denotes the columns to be taken.

    The tables are merged and entered row-wise, i.e each file is 
    a row unless *row_wise* is set to False. The latter is useful if 
    histograms are being merged.

    Filenames are stored in a ``track`` column. Directory names
    are chopped off.
    '''
    if len(infiles) == 0:
        raise ValueError( "no files for merging")
    if suffix:
        header = ",".join( [ os.path.basename( snip( x, suffix) ) for x in infiles] )
    elif regex:
        header = ",".join( [ "-".join(re.search( regex, x).groups()) for x in infiles] )        
    else:
        header = ",".join( [ os.path.basename( x ) for x in infiles] )

    columns = ",".join( map(str, [ x + 1 for x in columns ]))

    if infiles[0].endswith(".gz"):
        filenames = " ".join( [ "<( zcat %s | cut -f %s )" % (x,columns) for x in infiles ] )
    else:
        filenames = " ".join( [ "<( cat %s | cut -f %s )" % (x,columns) for x in infiles ] )

    tablename = toTable( outfile )

    if row_wise:
        transform = """| perl -p -e "s/bin/track/" | python %(scriptsdir)s/table2table.py --transpose""" % PARAMS
    else:
        transform = ""

    statement = """python %(scriptsdir)s/combine_tables.py
                      --headers=%(header)s
                      --skip-titles
                      --missing=0
                      --ignore-empty
                   %(filenames)s
                %(transform)s
                | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s
            """
    run()

def createView( dbhandle, tables, tablename, outfile, view_type = "TABLE",
                ignore_duplicates = True ):
    '''create a view in database for tables.

    Tables should be a list of tuples. Each tuple
    should contain the name of a table and the field
    to join with the first table.
    
    If ignore duplicates is set to False, duplicate column
    names will be added with the tablename as prefix. The
    default is to ignore.

    For example::
       tables = ("reads_summary", "track",
                 "bam_stats", "track",
                 "context_stats", "track",
                 "picard_stats_alignment_summary_metrics", "track" )

    view_type 
       type of view. If a view is to be created across multiple database,
       use "TABLE", otherwise, use "VIEW"

    '''

    Database.executewait( dbhandle, "DROP %(view_type)s IF EXISTS %(tablename)s" % locals() )

    tracks, columns = [], []
    tablenames = [ x[0] for x in tables ]
    for table, track in tables:
        d = Database.executewait( dbhandle,"SELECT COUNT(DISTINCT %s) FROM %s" % (track,table))
        tracks.append( d.fetchone()[0] )
        columns.append( [x.lower() for x in Database.getColumnNames( dbhandle, table ) if x != track ] )

    E.info( "creating %s from the following tables: %s" % (tablename, str(zip( tablenames, tracks ))))
    if min(tracks) != max(tracks):
        raise ValueError("number of rows not identical - will not create view" )

    from_statement = " , ".join( [ "%s as t%i" % (y[0],x) for x,y in enumerate(tables) ] )
    f = tables[0][1]
    where_statement = " AND ".join( ["t0.%s = t%i.%s" % (f,x+1,y[1]) for x,y in enumerate(tables[1:]) ] )
    
    all_columns, taken = [], set()
    for x, c in enumerate(columns):
        i = set(taken).intersection( set(c))
        if i:
            E.warn("duplicate column names: %s " % i )
            if not ignore_duplicates:
                table = tables[x][0]
                all_columns.extend( ["t%i.%s AS %s_%s" % (x,y,table,y) for y in i ] )
                c = [ y for y in c if y not in i ]
                
        all_columns.extend( ["t%i.%s" % (x,y) for y in c ] )
        taken.update( set(c) )

    all_columns = ",".join( all_columns)
    statement = '''CREATE %(view_type)s %(tablename)s AS SELECT t0.track, %(all_columns)s
                   FROM %(from_statement)s
                   WHERE %(where_statement)s
    ''' % locals()

    Database.executewait( dbhandle, statement )

    nrows = Database.executewait( dbhandle, "SELECT COUNT(*) FROM view_mapping" ).fetchone()[0]
    
    if nrows == 0:
        raise ValueError( "empty view mapping, check statement = %s" % (statement % locals()) )
    if nrows != min(tracks):
        E.warn( "view creates duplicate rows, got %i, expected %i" % (nrows, min(tracks)))

    E.info( "created view_mapping with %i rows" % nrows )

    touch( outfile )

##############################################################################################
##############################################################################################
##############################################################################################
## Remote job execution
##############################################################################################
def shellquote( statement ):
    '''shell quote a string to be used as a function argument.

    from http://stackoverflow.com/questions/967443/python-module-to-shellquote-unshellquote
    '''
    _quote_pos = re.compile('(?=[^-0-9a-zA-Z_./\n])')

    if statement:
        return _quote_pos.sub('\\\\', statement).replace('\n',"'\n'")
    else:
        return "''"

def getCallerLocals(decorators=0):
    '''returns locals of caller using frame.

    optional pass number of decorators
    
    from http://pylab.blogspot.com/2009/02/python-accessing-caller-locals-from.html
    '''
    f = sys._getframe(2+decorators)
    args = inspect.getargvalues(f)
    return args[3]

def getCaller(decorators=0):
    frm = inspect.stack()[2+decorators]
    mod = inspect.getmodule(frm[0])
    return mod

def execute( statement, **kwargs ):
    '''execute a statement locally.'''

    if not kwargs: kwargs = getCallerLocals()

    kwargs = dict( PARAMS.items() + kwargs.items() )    

    L.debug("running %s" % (statement % kwargs))

    process = subprocess.Popen(  statement % kwargs,
                                 cwd = os.getcwd(), 
                                 shell = True,
                                 stdin = subprocess.PIPE,
                                 stdout = subprocess.PIPE,
                                 stderr = subprocess.PIPE )

    # process.stdin.close()
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        raise PipelineError( "Child was terminated by signal %i: \nThe stderr was: \n%s\n%s\n" % (-process.returncode, stderr, statement ))

    return stdout, stderr

# Definition of helper functions for job scripts
# detect_pipe_error(): propagate error of programs not at the end of a pipe
# checkpoint(): exit a set of chained commands (via ;) if the previous command failed.
_exec_prefix = '''detect_pipe_error_helper() 
    {
    while [ "$#" != 0 ] ; do
        # there was an error in at least one program of the pipe
        if [ "$1" != 0 ] ; then return 1 ; fi
        shift 1
    done
    return 0 
    }
    detect_pipe_error() {
    detect_pipe_error_helper "${PIPESTATUS[@]}"
    return $?
    }
    checkpoint() {
        detect_pipe_error;
        if [ $? != 0 ]; then exit 1; fi;
    }
    '''

_exec_suffix = "; detect_pipe_error"

def buildStatement( **kwargs ):
    '''build statement from kwargs.

    Options in PARAMS are added, but kwargs take precedence.
    
    If outfile is in kwargs, 
    '''

    if "statement" not in kwargs:
        raise ValueError("'statement' not defined")

    local_params = substituteParameters( **kwargs )

    # build the statement
    try:
        statement = kwargs.get("statement") % local_params
    except KeyError, msg:
        raise KeyError( "Error when creating command: could not find %s in dictionaries" % msg)
    except ValueError, msg:
        raise ValueError( "Error when creating command: %s, statement = %s" % (msg, kwargs.get("statement") ) )

    # add bash as prefix to allow advanced shell syntax like 'wc -l <( gunzip < x.gz)'
    # executable option to call() does not work. Note that there will be an extra
    # indirection.
    statement = " ".join( re.sub( "\t+", " ", statement).split( "\n" ) ).strip()
    if statement.endswith(";"): statement = statement[:-1]

    L.debug( "running statement:\n%s" % statement )

    return statement

def expandStatement( statement ):
    '''add exec_prefix and exec_suffix to statement.'''
    
    return " ".join( (_exec_prefix, statement, _exec_suffix) )

def joinStatements( statements, infile ):
    '''join a chain of statements into a single statement.

    Each statement contains an @IN@ or a @OUT@ or both.
    These will be replaced by the names of successive temporary
    files.
    
    In the first statement, @IN@ is replaced with *infile*. 
    
    The last statement should move @IN@ to outfile.

    returns a single statement.
    '''
    
    prefix = getTempFilename()
    pattern = "%s_%%i" % prefix

    result = []
    for x, statement in enumerate(statements):
        if x == 0:
            s = re.sub( "@IN@", infile, statement )
        else:
            s = re.sub( "@IN@", pattern % x, statement )

        s = re.sub( "@OUT@", pattern % (x+1), s ).strip()

        if s.endswith(";"): s = s[:-1]
        result.append( s )

    assert prefix != ""
    result.append( "rm -f %s*" % prefix )
    
    result = "; checkpoint ; ".join( result )
    return result

def getStdoutStderr( stdout_path, stderr_path, tries=5 ):
    '''get stdout/stderr allowing for same lag.

    Try at most *tries* times. If unsuccessfull, throw PipelineError.

    Removes the files once they are read. 

    Returns tuple of stdout and stderr.
    '''
    x = tries
    while x >= 0:
        if os.path.exists( stdout_path ): break
        time.sleep(1)
        x -= 1
            
    x = tries
    while x >= 0:
        if os.path.exists( stderr_path ): break
        time.sleep(1)
        x -= 1

    try:
        stdout = open( stdout_path, "r" ).readlines()
    except IOError, msg:
        E.warn( "could not open stdout: %s" % msg )
        stdout = []
            
    try:
        stderr = open( stderr_path, "r" ).readlines()
    except IOError, msg:
        E.warn( "could not open stdout: %s" % msg )
        stderr = []

    try:
        os.unlink( stdout_path )
        os.unlink( stderr_path )
    except OSError, msg:
        pass

    return stdout, stderr

def run( **kwargs ):
    """run a command line statement (remotely).

    :param kwargs: dictionary with command line statement and options.

    The :py:func:`run` method runs one or more statements. The method
    will examine the contents of kwargs for the key ``statement``. If
    found, it will take the contents of statement and run it. If ``statement``
    is not found but ``statements`` containing a list of 
    command line statements, these will be executed instead.

    Several indicator variables permit setting job specific parameters:
    
    to_cluster
       send job to cluster. If not set or False, the command is executed 
       locally. Remote job execution is implemented using DRMAA.
      
    job_options
       job options. This is a string with job options that is passed to the 
       scheduler.
       
    job_queue
       queue to use for the job

    job_priority 
       priority to use for the job

    job_array
       run job as an array job. The value are three comma separated values
       denoting the range and increment. For example: ``job_array = 0,10,5``
       will execute jobs 0 to 9 in increments of 5.

    Before executing a statement, the statement is interpolated with the
    contents of the global :py:data`PARAMS` dictionary, any variables
    defined in the context of the calling function or those given in
    ``kwargs``.

    This method creates several files while running. The job script,
    stdout and stderr streams are temporary files in the current working
    directory. 

    The script writes a list of installed modules and environment variables
    to the file :file:`shell.log`.

    Troubleshooting:
       1. DRMAA creates sessions and there is a limited number
          of sessions available. If there are too many in use or sessions
          are not freed after failed jobs, use ``qconf -secl``
          to list sessions and ``qconf -kec #`` to delete sessions.
       2. Memory: 1G of free memory can be requested using 
          ``job_options="-l mem_free=1G"``.
          If there are error messages like "no available queue", then the
          problem could be that a particular complex attribute has
          not been defined (the code should be ``hc`` for ``host:complex``
          and not ``hl`` for ``host:local``. Note that qrsh/qsub directly
          still works.
    """

    # combine options using correct preference
    options = dict(PARAMS.items())
    options.update( getCallerLocals().items() )
    options.update( kwargs.items() )

    # setup a job/session with DRMAA
    def setupJob( session ):

        jt = session.createJobTemplate()
        jt.workingDirectory = os.getcwd()
        jt.jobEnvironment = { 'BASH_ENV' : '~/.bashrc' }
        jt.args = []
        jt.nativeSpecification = "-V -q %s -p %i -N %s %s" % \
            (options.get("job_queue", GLOBAL_OPTIONS.cluster_queue ),
             options.get("job_priority", GLOBAL_OPTIONS.cluster_priority ),
             "_" + re.sub( "[:]", "_", os.path.basename(options.get("outfile", "ruffus" ))),
             options.get("job_options", GLOBAL_OPTIONS.cluster_options))
            
        # keep stdout and stderr separate
        jt.joinFiles=False

        return jt

    shellfile = os.path.join( os.getcwd(), "shell.log" )

    pid = os.getpid()
    L.debug( 'task: pid = %i' % pid )

    # connect to global session 
    session = GLOBAL_SESSION
    L.debug( 'task: pid %i: sge session = %s' % (pid, str(session)))

    def buildJobScript( statement ):
        '''build job script from statement.

        returns (name_of_script, stdout_path, stderr_path)
        '''
        
        tmpfile = getTempFile( dir = os.getcwd() )
        tmpfile.write( "#!/bin/bash\n" ) #  -l -O expand_aliases\n" )
        tmpfile.write( 'echo "START--------------------------------" >> %s \n' % shellfile )
        # disabled - problems with quoting
        # tmpfile.write( '''echo 'statement=%s' >> %s\n''' % (shellquote(statement), shellfile) )
        tmpfile.write( "set &>> %s\n" % shellfile)
        tmpfile.write( "module list &>> %s\n" % shellfile )
        tmpfile.write( 'echo "END----------------------------------" >> %s \n' % shellfile )
        tmpfile.write( expandStatement( statement ) + "\n" )
        tmpfile.close()

        job_path = os.path.abspath( tmpfile.name )
        stdout_path = job_path + ".stdout" 
        stderr_path = job_path + ".stderr" 

        os.chmod( job_path, stat.S_IRWXG | stat.S_IRWXU )

        return (job_path, stdout_path, stderr_path)
    
    # run multiple jobs
    if options.get( "statements" ):

        statement_list = []
        for statement in options.get("statements"): 
            options["statement"] = statement
            statement_list.append(buildStatement( **options))
            
        if options.get( "dryrun", False ): return

        jt = setupJob( session )
        
        jobids, filenames = [], []
        for statement in statement_list:

            job_path, stdout_path, stderr_path = buildJobScript( statement )

            jt.remoteCommand = job_path
            jt.outputPath=":"+ stdout_path
            jt.errorPath=":" + stderr_path

            os.chmod( job_path, stat.S_IRWXG | stat.S_IRWXU )

            jobid = session.runJob(jt)
            jobids.append( jobid )
            filenames.append( (job_path, stdout_path, stderr_path) )

            L.debug( "job has been submitted with jobid %s" % str(jobid ))
        
        L.debug( "waiting for %i jobs to finish " % len(jobids) )
        session.synchronize(jobids, drmaa.Session.TIMEOUT_WAIT_FOREVER, False)
        
        # collect and clean up
        for jobid, statement, paths in zip( jobids, statement_list, filenames) :
            job_path, stdout_path, stderr_path = paths
            retval = session.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)

            stdout, stderr = getStdoutStderr( stdout_path, stderr_path )

            if retval.exitStatus != 0:
                raise PipelineError( "---------------------------------------\n"
                                     "Child was terminated by signal %i: \n"
                                     "The stderr was: \n%s\n%s\n" 
                                     "---------------------------------------\n" % \
                                         (retval.exitStatus, 
                                          "".join( stderr),
                                          statement ) )

            try:
                os.unlink( job_path )
            except OSError:
                L.warn( "temporary job file %s not present for clean-up - ignored" % job_path )
            
        session.deleteJobTemplate(jt)

    # run a single parallel job
    elif (options.get( "job_queue" ) or options.get( "to_cluster" )) \
            and not GLOBAL_OPTIONS.without_cluster:

        statement = buildStatement( **options )

        if options.get( "dryrun", False ): return

        job_path, stdout_path, stderr_path = buildJobScript( statement )

        jt = setupJob( session )
        
        jt.remoteCommand = job_path
        # later: allow redirection of stdout and stderr to files; can even be across hosts?
        jt.outputPath=":"+ stdout_path
        jt.errorPath=":" + stderr_path

        if "job_array" in options and options["job_array"] != None:
            # run an array job
            start, end, increment = options.get("job_array" )
            L.debug("starting an array job: %i-%i,%i" % (start, end, increment ))
            # sge works with 1-based, closed intervals
            jobids = session.runBulkJobs( jt, start+1, end, increment )
            L.debug( "%i array jobs have been submitted as jobid %s" % (len(jobids), jobids[0]) )
            retval = session.synchronize(jobids, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)
        else:
            jobid = session.runJob(jt)
            L.debug( "job has been submitted with jobid %s" % str(jobid ))
            try:
                retval = session.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)
            except Exception, msg:
                # ignore message 24 in PBS
                # code 24: drmaa: Job finished but resource usage information and/or termination status could not be provided.":
                if not msg.message.startswith("code 24"): raise
                retval = None

        stdout, stderr = getStdoutStderr( stdout_path, stderr_path )

        if "job_array" not in options:
            if retval and retval.exitStatus != 0:
                raise PipelineError( "---------------------------------------\n"
                                     "Child was terminated by signal %i: \n"
                                     "The stderr was: \n%s\n%s\n"
                                     "-----------------------------------------" % \
                                         (retval.exitStatus, 
                                          "".join( stderr), statement))
            
        session.deleteJobTemplate(jt)
        try:
            os.unlink( job_path )
        except OSError:
            L.warn( "temporary job file %s not present for clean-up - ignored" % job_path )
    else:
        statement = buildStatement( **options )

        if options.get( "dryrun", False ): return
 
        if "<(" in statement:
            if "'" in statement: raise ValueError( "advanced bash syntax combined with single quotes" )
            statement = """/bin/bash -c '%s'""" % statement

        process = subprocess.Popen(  expandStatement( statement ),
                                     cwd = os.getcwd(), 
                                     shell = True,
                                     stdin = subprocess.PIPE,
                                     stdout = subprocess.PIPE,
                                     stderr = subprocess.PIPE )

        # process.stdin.close()
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            raise PipelineError( "---------------------------------------\n"
                                 "Child was terminated by signal %i: \n"
                                 "The stderr was: \n%s\n%s\n"
                                 "-----------------------------------------" % \
                                     (-process.returncode, stderr, statement ))

def submit( module, 
            function, 
            params = None,
            infiles = None, 
            outfiles = None, 
            toCluster = True,
            logfile = None ):
    '''submit a python *function* as a job to the cluster.

    :param module: Filename of module
    :param function: function within module to execute
    :param params: list of parameters to function
    :param infiles: list of input filenames
    :param outfiles: list of output filenames
    :param toCluster: if true, execute remotely
    :param logfile: logging output goes to this file.

    The function should reside in *module*. If *module* is
    not part of the PYTHONPATH, an absolute path can be given.

    *params* and *infiles*/*outfiles* are orthogonal ways to
    pass parameters to *function*. *params* is the generic method,
    while the latter corresponds to ruffus syntax.
    
    *infiles* and *output* are either a single filename or a list of 
    input/output filenames. Neither options supports yet nested lists.

    This method calls :doc:`../scripts/run_function` via :py:func:`.run`.
    '''

    if type( infiles ) in (list, tuple):
        infiles = " ".join( ["--input=%s" % x for x in infiles ] )
    else:
        infiles = "--input=%s" % infiles

    if type( outfiles ) in (list, tuple):
        outfiles = " ".join( ["--output=%s" % x for x in outfiles ] )
    else:
        outfiles = "--output=%s" % outfiles

    if logfile:
        logfile = "--log=%s" % logfile
    else:
        logfile = ""

    if params:
        params = "--params=%s" % ",".join(params)
    else:
        params = ""

    to_cluster = toCluster

    statement = '''python %(scriptsdir)s/run_function.py
                          --module=%(module)s
                          --function=%(function)s
                          %(logfile)s
                          %(infiles)s
                          %(outfiles)s
                          %(params)s
                '''
    run()

##############################################################################################
##############################################################################################
##############################################################################################
def clonePipeline( srcdir ):
    '''clone a pipeline from *srdir* into the current directory.

    :param srcdir: directory with pipeline to clone.

    The purpose of cloning is to create a copy of a pipeline with
    a minimum amount of duplicating data. Usually one would clone
    a pipeline, remove some output files and then re-run the analysis 
    from a certain point with changed parameters or code.

    Cloning involves creating copies for all configuration files and
    the database while creating soft links for all other files.
    
    Report related and temporary files are ignored.
    '''
    
    destdir = os.path.curdir

    copy_files = ("sphinxreport.ini", "conf.py", "pipeline.ini", "csvdb" )
    ignore_prefix = ("report", "_cache", "export", "tmp", "ctmp", "_static", "_templates" )

    def _ignore( p ):
        for x in ignore_prefix:
            if p.startswith( x ): 
                return True
        return False

    for root, dirs, files in os.walk(srcdir):

        relpath = os.path.relpath( root, srcdir )
        if _ignore( relpath ): continue

        for d in dirs:
            if _ignore( d ): continue
            dest = os.path.join( os.path.join(destdir, relpath, d ) )
            os.mkdir( dest )
            # touch
            s = os.stat( os.path.join(root, d ) )
            os.utime( dest, (s.st_atime, s.st_mtime ))

        for f in files:
            if _ignore( f ): continue

            fn = os.path.join( root, f )
            dest_fn = os.path.join( destdir, relpath, f ) 
            if f in copy_files:
                shutil.copyfile( fn, dest_fn )
            else:
                # realpath resolves links - thus links will be linked to
                # the original target
                os.symlink( os.path.realpath( fn),
                            dest_fn )

##############################################################################################
##############################################################################################
##############################################################################################
def writeConfigFiles( path ):
    '''output configuration files from files in *path*.
    
    :param path: Directory with configuration files.
    '''

    for dest in ( "pipeline.ini", "sphinxreport.ini", "conf.py" ):
        src = os.path.join( path, dest)
        if os.path.exists(dest):
            L.warn( "file `%s` already exists - skipped" % dest )
            continue

        if not os.path.exists( src ):
            raise ValueError( "default config file `%s` not found"  % src )
        shutil.copyfile( src, dest )
        L.info( "created new configuration file `%s` " % dest )


##############################################################################################
##############################################################################################
##############################################################################################
def clean( patterns, dry_run = False ):
    '''clean up files given by glob *patterns*.

    returns list of files deleted together with their statinfo.
    '''

    cleaned = []

    for p in patterns:
        files = glob.glob( p )
        for x in files:
            statinfo = os.stat( x )
            cleaned.append( (x, statinfo) )
            if dry_run: continue
            os.unlink( x )
        L.info( "%i files: %s" % (len(files), p ))

    return cleaned

#############################################################################
#############################################################################
#############################################################################
def peekParameters( workingdir, pipeline ):
    '''peak configuration parameters from an external directory.
    '''
    
    dirname = os.path.dirname( __file__ )

    if dirname == "":
        # special case: pipeline called in source directory
        dirname = os.path.abspath(".")
    else:
        # point towards location of CGAT pipelines
        dirname = os.path.join( os.path.dirname( dirname ), "CGATPipelines")

    pipeline = os.path.join( dirname, pipeline )
    assert os.path.exists( pipeline ), "can't find pipeline source %s" % pipeline
    if workingdir == "": workingdir = os.path.abspath(".")

    assert os.path.exists( workingdir ), "can't find working dir %s" % workingdir

    statement = "python %s -f -v 0 dump" % pipeline
    process = subprocess.Popen(  statement,
                                 cwd = workingdir, 
                                 shell = True,
                                 stdin = subprocess.PIPE,
                                 stdout = subprocess.PIPE,
                                 stderr = subprocess.PIPE )

    # process.stdin.close()
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        raise PipelineError( "Child was terminated by signal %i: \nThe stderr was: \n%s\n" % (-process.returncode, stderr ))

    for line in stdout.split("\n"):
        if line.startswith("dump"):
            exec( line )

    return dump

#############################################################################
#############################################################################
#############################################################################
def run_report( clean = True):
    '''run sphinxreport.'''

    dirname, basename = os.path.split( getCaller().__file__ )

    docdir = os.path.join( dirname, "pipeline_docs", snip( basename, ".py" ) )
    themedir = os.path.join( dirname, "pipeline_docs", "themes")
    relpath = os.path.relpath( docdir )
    trackerdir = os.path.join( docdir, "trackers" )

    # Requires an X11 connection to cluster nodes
    # A solution is to run xvfb on the nodes.
    to_cluster = False
    
    job_options= "-pe dedicated %i -R y" % PARAMS["report_threads"]

    # use a fake X display in order to avoid windows popping up
    # from R plots.
    xvfb_command = which("xvfb-run" )

    # permit multiple servers using -a option
    if xvfb_command: xvfb_command += " -a "
    else: xvfb_command = ""

    # if there is no DISPLAY variable set, xvfb runs, but
    # exits with error when killing process. Thus, ignore return
    # value.
    if not os.getenv("DISPLAY"):
        erase_return = "|| true"
    else:
        erase_return = ""

    if clean: clean = """rm -rf report _cache _static;"""
    else: clean = ""

    statement = '''
    %(clean)s
    ( export SPHINX_DOCSDIR=%(docdir)s; 
      export SPHINX_THEMEDIR=%(themedir)s; 
    %(xvfb_command)s
    sphinxreport-build 
           --num-jobs=%(report_threads)s
           sphinx-build 
                    -b html 
                    -d %(report_doctrees)s
                    -c . 
           %(docdir)s %(report_html)s
    >& report.log %(erase_return)s )
    '''

    run()

#############################################################################
#############################################################################
#############################################################################
class MultiLineFormatter(logging.Formatter):
    '''logfile formatter: add identation for multi-line entries.'''

    def format(self, record):
        s = logging.Formatter.format(self, record)
        if record.message:
            header, footer = s.split(record.message)
            s = s.replace('\n', '\n' + ' '*len(header))
        return s


#############################################################################
#############################################################################
#############################################################################
USAGE = '''
usage: %prog [OPTIONS] [CMD] [target]

Execute pipeline %prog.

Commands can be any of the following

make <target>
   run all tasks required to build *target*

show <target>
   show tasks required to build *target* without executing them

plot <target>
   plot image (using inkscape) of pipeline state for *target*

config
   write new configuration files pipeline.ini, sphinxreport.ini and conf.py
   with default values

dump
   write pipeline configuration to stdout

touch
   touch files only, do not run

clone <source>
   create a clone of a pipeline in <source> in the current directory. The cloning
   process aims to use soft linking to files (not directories) as much as possible. 
   Time stamps are preserved. Cloning is useful if a pipeline needs to be re-run
   from a certain point but the original pipeline should be preserved.

'''

def main( args = sys.argv ):

    global GLOBAL_OPTIONS
    global GLOBAL_ARGS
    global GLOBAL_SESSION

    parser = E.OptionParser( version = "%prog version: $Id: Pipeline.py 2799 2009-10-22 13:40:13Z andreas $",
                                    usage = USAGE )
    
    parser.add_option( "--pipeline-action", dest="pipeline_action", type="choice",
                       choices=("make", "show", "plot", "dump", "config", "clone" ),
                       help="action to take [default=%default]." )

    parser.add_option( "--pipeline-format", dest="pipeline_format", type="choice",
                      choices=( "dot", "jpg", "svg", "ps", "png" ),
                      help="pipeline format [default=%default]." )

    parser.add_option( "-n", "--dry-run", dest="dry_run", action="store_true",
                      help="perform a dry run (do not execute any shell commands) [default=%default]." )

    parser.add_option( "-f", "--force", dest="force", action="store_true",
                      help="force running the pipeline even if there are uncommited changes "
                           "in the repository [default=%default]." )

    parser.add_option( "-l", "--local", dest="without_cluster", action="store_true",
                      help="execute all jobs locally [default=%default]." )

    parser.add_option( "-p", "--multiprocess", dest="multiprocess", type="int",
                       help="number of parallel processes to use (different from number of jobs to use for cluster jobs) [default=%default]." ) 

    parser.add_option( "-t", "--tempdir", dest="tempdir", type="string",
                       help="temporary directory to use [default=%default].")

    parser.set_defaults(
        pipeline_action = None,
        pipeline_format = "svg",
        pipeline_targets = [],
        multiprocess = 2,
        logfile = "pipeline.log",
        dry_run = False,
        without_cluster = False,
        force = False,
        )

    (options, args) = E.Start( parser, 
                               add_cluster_options = True )


    GLOBAL_OPTIONS, GLOBAL_ARGS = options, args
    PARAMS["dryrun"] = options.dry_run
    
    # get mercurial version
    repo = hgapi.Repo( PARAMS["scriptsdir"] )
    version = repo.hg_id()

    # check for uncommitted changes in the repository
    status = repo.hg_status()
    if status["M"] or status["A"]:
        if not options.force:
            raise ValueError( "uncommitted change in code repository at '%s'. Either commit or use --force" % PARAMS["scriptsdir"])
        else:
            E.warn( "uncommitted changes in code repository - ignored ")

    if args: 
        options.pipeline_action = args[0]
        if len(args) > 1:
            options.pipeline_targets.extend( args[1:] )

    if options.pipeline_action in ("make", "show", "svg", "plot", "touch" ):

        try:
            if options.pipeline_action == "make":

                # create the session proxy
                GLOBAL_SESSION = drmaa.Session()
                GLOBAL_SESSION.initialize()
                
                #
                #   make sure we are not logging at the same time in different processes
                #
                #session_mutex = manager.Lock()

                # set up extra file logger
                handler = logging.FileHandler( filename = options.logfile, 
                                               mode = "a" )
                handler.setFormatter( MultiLineFormatter( '%(asctime)s %(levelname)s %(module)s.%(funcName)s.%(lineno)d %(message)s' ) )
                logger = logging.getLogger()
                logger.addHandler( handler )
                
                L.info( E.GetHeader() )
                L.info( "code location: %s" % PARAMS["scriptsdir"] )
                L.info( "code version: %s" % version[:-1] )

                pipeline_run( options.pipeline_targets, 
                              multiprocess = options.multiprocess, 
                              logger = logger,
                              verbose = options.loglevel )

                L.info( E.GetFooter() )

                GLOBAL_SESSION.exit()

            elif options.pipeline_action == "show":
                pipeline_printout( options.stdout, options.pipeline_targets, verbose = options.loglevel )

            elif options.pipeline_action == "touch":
                pipeline_run( options.pipeline_targets, 
                              touch_files_only = True,
                              verbose = options.loglevel )

            elif options.pipeline_action == "svg":
                pipeline_printout_graph( options.stdout, 
                                         options.pipeline_format,
                                         options.pipeline_targets )

            elif options.pipeline_action == "plot":
                outf, filename = tempfile.mkstemp()
                pipeline_printout_graph( os.fdopen(outf,"w"),
                                         options.pipeline_format,
                                         options.pipeline_targets )
                execute( "inkscape %s" % filename ) 
                os.unlink( filename )

        except ruffus_exceptions.RethrownJobError, value:
            E.error("some tasks resulted in errors - error messages follow below" )
            # print value
            E.error( value )
            E.error( "end of error messages" )
            raise

    elif options.pipeline_action == "dump":
        # convert to normal dictionary (not defaultdict) for parsing purposes
        print "dump = %s" % str(dict(PARAMS))

    elif options.pipeline_action == "config":
        f = sys._getframe(1)
        caller = inspect.getargvalues(f).locals["__file__"]
        prefix = os.path.splitext(caller)[0]
        writeConfigFiles( prefix )

    elif options.pipeline_action == "clone":
        clonePipeline( options.pipeline_targets[0] )

    else:
        raise ValueError("unknown pipeline action %s" % options.pipeline_action )

    E.Stop()

if __name__ == "__main__":

    main()

    #parser = E.OptionParser( version = "%prog version: $Id: Pipeline.py 2799 2009-10-22 13:40:13Z andreas $")

    #(options, args) = E.Start( parser, add_cluster_options = True )

    #global GLOBAL_OPTIONS
    #global GLOBAL_ARGS
    #GLOBAL_OPTIONS, GLOBAL_ARGS = options, args

    #L.info("in main")

    #run( **{"statement" : "printenv > test.out", "job_queue" : "server_jobs.q", "job_priority" : -10 } )
    #run( **{"statement" : "printenv > test2.out", "job_queue" : "server_jobs.q", "job_priority" : -10 } )
    #run( **{"statement" : "printenv > test3.out", "job_queue" : "server_jobs.q", "job_priority" : -10 } )
    
