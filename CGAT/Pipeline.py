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

The :mod:`Pipeline` module contains various utility functions
for interfacing CGAT ruffus pipelines with databases and the 
cluster.

API
----

'''
import os, sys, re, subprocess, optparse, stat, tempfile
import time, random, inspect, types, multiprocessing
import logging, collections, shutil, glob, gzip
import ConfigParser
import Database

# talking to a cluster
try:
    import drmaa
    HAS_DRMAA = True
except RuntimeError:
    HAS_DRMAA = False

# talking to mercurial
import hgapi

# CGAT specific options - later to be removed
from Local import *

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

class PipelineError( Exception ): pass

# Sort out the important paths
LIB_DIR=os.path.dirname( __file__ )
ROOT_DIR=os.path.dirname( LIB_DIR )
SCRIPTS_DIR=os.path.join( ROOT_DIR, "scripts" )
# if Pipeline.py is from an installed version,
# scripts are located in the "bin" directory.
if not os.path.exists( SCRIPTS_DIR):
    SCRIPTS_DIR = os.path.join(sys.exec_prefix, "bin")
    
# possible to use defaultdict, but then statements will
# fail on execution if a parameter does not exists, and not
# while building the statement. Hence, use dict.
PARAMS = { 
    'scriptsdir' : SCRIPTS_DIR,
    'toolsdir' : SCRIPTS_DIR,
    'cmd-farm' : """python %s/farm.py 
                --method=drmaa 
                --cluster-priority=-10 
		--cluster-queue=all.q 
		--cluster-num-jobs=100 
                --bashrc=%s/bashrc.cgat
		--cluster-options="" """ % (SCRIPTS_DIR,SCRIPTS_DIR),
    'cmd-sql' : """sqlite3 -header -csv -separator $'\\t' """,
    'cmd-run' : """%s/run.py""" % SCRIPTS_DIR
    }

# path until parameter sharing is resolved between CGAT module
# and the pipelines module.
import Local
Local.PARAMS = PARAMS

hostname = os.uname()[0]

CONFIG = {}

# local temporary directory to use
TMPDIR = '/scratch'

def configToDictionary( config ):

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

    Sections and keys are combined with an underscore. If
    a key without section does not exist, it will be added 
    plain.

    For example::

       [general]
       input=input1.file

       [special]
       input=input2.file

    will be entered as { 'general_input' : "input1.file",
    'input: "input1.file", 'special_input' : "input2.file" }

    This function also updates the module-wide parameter map.
    
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
    
    Parameters are processed in the same way as
    getParameters, but the global parameter dictionary
    is not updated.
    '''
    config = ConfigParser.ConfigParser()
    config.read( filenames )

    p = configToDictionary( config )
    return p

def substituteParameters( **kwargs ):
    '''return a local PARAMS dictionary.

    Options in ``**kwargs`` substitute default
    values in PARAMS.

    Finally, task specific configuration values 
    are inserted.
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

def checkFiles( filenames ):
    """check for the presence/absence of files"""
    
    missing = []
    for filename in filenames:
        if not os.path.exists( filename ):
            missing.append( filename )

    if missing:
        raise ValueError( "missing files: %s" % ",".join(missing) )

def which( filename ):

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

def clone( infile, outfile ):
    '''create a clone of ``infile`` named ``outfile``
    by creating a soft-link.
    '''
    # link via relative paths, otherwise it 
    # fails if infile and outfile are in different
    # directories or in a subdirectory
    if os.path.dirname(infile) != os.path.dirname(outfile):
        relpath = os.path.relpath( os.path.dirname(infile), os.path.dirname(outfile) )
    else:
        relpath = "."
    target = os.path.join( relpath, os.path.basename( infile ) )

    try:
        os.symlink( target, outfile )
    except OSError:
        pass
    
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
    
def getTempFilename( dir = TMPDIR ):
    '''get a temporary filename.

    The caller needs to delete it.
    '''
    tmpfile = tempfile.NamedTemporaryFile(dir=dir,delete=False, prefix = "ctmp" )
    tmpfile.close()
    return tmpfile.name

def getTempFile( dir = TMPDIR ):
    '''get a temporary file.

    The caller needs to delete it.
    '''
    return tempfile.NamedTemporaryFile(dir=dir, delete=False, prefix = "ctmp" )

def getTempDir( dir = TMPDIR ):
    return tempfile.mkdtemp( dir=dir, prefix = "ctmp" )

def checkExecutables( filenames ):
    """check for the presence/absence of executables in the path"""
    
    missing = []

    for filename in filenames:
        if not which( filename ):
            missing.append( filename )

    if missing:
        raise ValueError( "missing executables: %s" % ",".join(missing) )

def checkScripts( filenames ):
    """check for the presence/absence of scripts"""
    missing = []
    for filename in filenames:
        if not os.path.exists( filename ):
            missing.append( filename )

    if missing:
        raise ValueError( "missing scripts: %s" % ",".join(missing) )

def checkParameter( key ):
    if key not in PARAMS:
        raise ValueError("need `%s` to be set" % key )

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

def isEmpty( filename ):
    '''return true if file *filename* is empty.
    
    If filename ends in .gz, it checks if the contents are empty
    by opening it.
    '''
    if filename.endswith(".gz"):
        n = 0
        with gzip.open(filename) as inf:
            return len(inf.read(10)) == 0
    else:
        return os.stat(filename)[6]==0

def asList( param ):
    '''return a param as a list'''
    if type(param) == str:
        try:
            params = [x.strip() for x in param.strip().split(",")]
        except AttributeError:
            params = [param.strip()]
        return [ x for x in params if x != ""]
    elif type(param) in (types.ListType, types.TupleType):
        return param
    else: 
        return [param]

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

def quote( track ):
    '''quote track such that is applicable for a table name.'''
    return re.sub( "[-(),\[\].]", "_", track)

def shellquote( statement ):
    '''shell quote a string to be used as a function argument.

    from http://stackoverflow.com/questions/967443/python-module-to-shellquote-unshellquote
    '''
    _quote_pos = re.compile('(?=[^-0-9a-zA-Z_./\n])')

    if statement:
        return _quote_pos.sub('\\\\', statement).replace('\n',"'\n'")
    else:
        return "''"

def toTable( outfile ):
    '''convert an outfile (filename) into
    a table name.

    The table name is quoted.
    '''
    assert outfile.endswith( ".load" ) 
    name = os.path.basename( outfile[:-len(".load")] )
    return quote( name )

def load( infile, 
          outfile = None, 
          options = "", 
          collapse = None,
          transpose = None,
          tablename = None):
    '''straight import from tab separated table.

    The table name is given by outfile without the
    ".load" suffix.

    If *collapse* is set, the table will be collapsed before loading.
    The value of collapse is the value used for missing values.

    If *transpose* is set, the table will be transposed before
    loading.  The first column in the first row will be set to the
    string within transpose.
    '''

    if not tablename:
        tablename = toTable( outfile )
    
    statement = []
    if infile.endswith(".gz"): statement.append( "zcat %(infile)s" )
    else: statement.append( "cat %(infile)s" )

    if collapse != None:
        statement.append( "python %(scriptsdir)s/table2table.py --collapse=%(collapse)s" )

    if transpose != None:
        statement.append( "python %(scriptsdir)s/table2table.py --transpose --set-transpose-field=%(transpose)s" )

    statement.append('''
    python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              %(options)s 
              --table=%(tablename)s 
    > %(outfile)s
    ''')

    statement = " | ".join( statement)

    run()

def concatenateAndLoad( infiles, 
                        outfile, 
                        regex_filename = None, 
                        header = None, 
                        cat = None, 
                        has_titles = True, 
                        missing_value = "na",
                        options = "" ):
    '''concatenate categorical tables and load into a database.

    If *has_titles* is False, the tables are assumed to have no titles.
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
        
    if has_titles == False: no_titles = "--no-titles"
    else: no_titles = ""

    options = " ".join(options)
    load_options = " ".join(load_options) + " " + passed_options
    statement = '''python %(scriptsdir)s/combine_tables.py
                     --cat=%(cat)s
                     --missing-value=%(missing_value)s
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
                  row_wise = True,
                  options = "",
                  prefixes = None ):
    '''merge categorical tables and load into a database.

    The tables are merged and entered row-wise, i.e each file is a
    row. If *row_wise* is set to False, each table will be a column in
    the resulting table. This is useful if histograms are being
    merged.

    *columns* denotes the columns to be taken. By default, the first
    two columns are taken with the first being the key. Filenames are 
    stored in a ``track`` column. Directory names are chopped off.

    If *columns* is set to None, all columns will be taken. Here,
    column names will receive a prefix (*prefixes*). If *prefixes* is
    None, the filename will be added as a prefix. 

    If *prefixes* is a list, the respective prefix will be added to
    each column. The length of *prefixes* and *infiles* need to be the
    same.

    The track names are given by the filenames without any paths. If
    *suffix* is given, the suffix will be removed. If *regex* is set,
    the full filename will be used to extract a track name via the
    supplied regular expression.

    *options* are passed on to ``csv2db.py``.
    '''
    if len(infiles) == 0:
        raise ValueError( "no files for merging")

    if suffix:
        header = ",".join( [ os.path.basename( snip( x, suffix) ) for x in infiles] )
    elif regex:
        header = ",".join( [ "-".join(re.search( regex, x).groups()) for x in infiles] )        
    else:
        header = ",".join( [ os.path.basename( x ) for x in infiles] )

    header_stmt = "--header=%s" % header

    if columns:
        column_filter = "| cut -f %s" % ",".join( map(str, [ x + 1 for x in columns ]))
    else:
        column_filter = ""
        if prefixes:
            assert len(prefixes) == len(infiles)
            header_stmt = "--prefixes=%s" % ",".join( prefixes)
        else:
            header_stmt = "--add-file-prefix"

    if infiles[0].endswith(".gz"):
        filenames = " ".join( [ "<( zcat %s %s )" % (x,column_filter) for x in infiles ] )
    else:
        filenames = " ".join( [ "<( cat %s %s )" % (x,column_filter) for x in infiles ] )

    tablename = toTable( outfile )

    if row_wise:
        transform = """| perl -p -e "s/bin/track/" | python %(scriptsdir)s/table2table.py --transpose""" % PARAMS
    else:
        transform = ""

    statement = """python %(scriptsdir)s/combine_tables.py
                      %(header_stmt)s
                      --skip-titles
                      --missing=0
                      --ignore-empty
                   %(filenames)s
                %(transform)s
                | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
                      --index=track
                      --table=%(tablename)s 
                      %(options)s
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


def snip( filename, extension = None, alt_extension = None, path = False):
    '''return prefix of filename.

    If extension is given, make sure that filename has the extension.

    If path is set to false, the path is stripped from the file name.
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
    
    if "cwd" not in kwargs:
        cwd = os.getcwd()
    else:
        cwd = kwargs["cwd"]

    process = subprocess.Popen(  statement % kwargs,
                                 cwd = cwd, 
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
    
    prefix = getTempFilename(".")
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
    """run a statement.

    runs it on the cluster using drmaa if to_cluster is set.

    Troubleshooting:
       1. DRMAA creates sessions and their is a limited number
          of sessions available. If there are two many or sessions
          become not available after failed jobs, use ``qconf -secl``
          to list sessions and ``qconf -kec #`` to delete sessions.
       2. Memory: 1G of free memory can be requested using ``-l mem_free=1G``.
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

    # Run a single parallel job if
    #   1. job_queue is set, or
    #   2. to_cluster is not defined or to_cluster is set to True.
    # If the cluster has not been disabled through the command line, do not
    #     run on cluster
    elif (options.get( "job_queue" ) or 
          ("to_cluster" not in options or options.get( "to_cluster" ))) \
          and (GLOBAL_OPTIONS and not GLOBAL_OPTIONS.without_cluster):

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

class MultiLineFormatter(logging.Formatter):
    '''logfile formatter: add identation for multi-line entries.'''

    def format(self, record):
        s = logging.Formatter.format(self, record)
        if record.message:
            header, footer = s.split(record.message)
            s = s.replace('\n', '\n' + ' '*len(header))
        return s

def submit( module, function, params = None,
            infiles = None, outfiles = None, 
            toCluster = True,
            logfile = None, 
            jobOptions = ""):
    '''Submit a python *function* as a job to the cluster.

    The function should reside in *module*. If *module* is
    not part of the PYTHONPATH, an absolute path can be given.

    *infiles* and *output* are either a single filename or a list of 
    input/output filenames. Neither options supports yet nested lists.
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

    job_options = jobOptions 

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

def clonePipeline( srcdir ):
    '''clone a pipeline from srcdir into the current directory.
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

def writeConfigFiles( path ):
    
    for dest in ( "pipeline.ini", "sphinxreport.ini", "conf.py" ):
        src = os.path.join( path, dest)
        if os.path.exists(dest):
            L.warn( "file `%s` already exists - skipped" % dest )
            continue

        if not os.path.exists( src ):
            raise ValueError( "default config file `%s` not found"  % src )
        shutil.copyfile( src, dest )
        L.info( "created new configuration file `%s` " % dest )



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

def peekParameters( workingdir, pipeline ):
    '''peak configuration parameters from an external directory.
    '''
    
    # attempt to locate directory with pipeline source code
    # This is a patch as pipelines might be called within the repository directory
    # or from an installed location
    dirname = os.path.dirname( __file__ )

    # called without a directory, use current directory
    if dirname == "":
        dirname = os.path.abspath(".")
    else:
        # else: use location of Pipeline.py
        # remove CGAT part, add CGATPipelines
        dirname = os.path.join( os.path.dirname(dirname), "CGATPipelines")
        # if not exists, assume we want version located
        # in directory of calling script.
        if not os.path.exists( dirname ):
            # directory is path of calling script
            v = getCallerLocals()
            dirname = os.path.dirname(v['__file__'])

    pipeline = os.path.join( dirname, pipeline )
    assert os.path.exists( pipeline ), "can't find pipeline source %s" % ( dirname, pipeline )
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
    print os.getenv("DISPLAY"), "command=", xvfb_command
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

    parser = optparse.OptionParser( version = "%prog version: $Id: Pipeline.py 2799 2009-10-22 13:40:13Z andreas $",
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

    parser.add_option( "-e", "--exceptions", dest="log_exceptions", action="store_true",
                      help="echo exceptions immediately as they occur [default=%default]." )

    parser.add_option( "-i", "--terminate", dest="terminate", action="store_true",
                      help="terminate immediately at the first exception [default=%default]." )

    parser.set_defaults(
        pipeline_action = None,
        pipeline_format = "svg",
        pipeline_targets = [],
        multiprocess = 2,
        logfile = "pipeline.log",
        dry_run = False,
        without_cluster = False,
        force = False,
        log_exceptions = False,
        exceptions_terminate_immediately = False,
        )

    (options, args) = E.Start( parser, 
                               add_cluster_options = True )

    GLOBAL_OPTIONS, GLOBAL_ARGS = options, args
    PARAMS["dryrun"] = options.dry_run

    version = None

    try:
        # this is for backwards compatibility
        # get mercurial version
        repo = hgapi.Repo( PARAMS["scriptsdir"] )
        version = repo.hg_id()

        status = repo.hg_status()
        if status["M"] or status["A"]:
            if not options.force:
                raise ValueError( "uncommitted change in code repository at '%s'. Either commit or use --force" % PARAMS["scriptsdir"])
            else:
                E.warn( "uncommitted changes in code repository - ignored ")
        version = version[:-1]
    except:
        # try git:
        try:
            stdout, stderr = execute( "git rev-parse HEAD", cwd = PARAMS["scriptsdir"] )
        except: 
            stdout = "NA"
        version = stdout

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
                L.info( "code version: %s" % version )

                pipeline_run( options.pipeline_targets, 
                              multiprocess = options.multiprocess, 
                              logger = logger,
                              verbose = options.loglevel,
                              log_exceptions = options.log_exceptions,
                              exceptions_terminate_immediately = options.exceptions_terminate_immediately,
                              )

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

    #parser = optparse.OptionParser( version = "%prog version: $Id: Pipeline.py 2799 2009-10-22 13:40:13Z andreas $")

    #(options, args) = E.Start( parser, add_cluster_options = True )

    #global GLOBAL_OPTIONS
    #global GLOBAL_ARGS
    #GLOBAL_OPTIONS, GLOBAL_ARGS = options, args

    #L.info("in main")

    #run( **{"statement" : "printenv > test.out", "job_queue" : "server_jobs.q", "job_priority" : -10 } )
    #run( **{"statement" : "printenv > test2.out", "job_queue" : "server_jobs.q", "job_priority" : -10 } )
    #run( **{"statement" : "printenv > test3.out", "job_queue" : "server_jobs.q", "job_priority" : -10 } )
    
