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
import os, sys, re, subprocess, optparse, stat, tempfile, time, random, inspect, types
import logging, collections, shutil, glob, gzip
import ConfigParser

import drmaa
from ruffus import *
import logging as L
import Experiment as E
import IOTools

global_options, global_args = None, None

# dictionary mapping process ids (pids) to drmaa sessions
global_sessions = {}

class PipelineError( Exception ): pass

# possible to use defaultdict, but then statements will
# fail on execution if a parameter does not exists, and not
# while building the statement. Hence, use dict.
PARAMS= { 
    'scriptsdir' : os.path.dirname( __file__ ),
    'toolsdir' : os.path.dirname( __file__ ),
    'cmd-farm' : """%s/farm.py 
                --method=drmaa 
                --cluster-priority=-10 
		--cluster-queue=all.q 
		--cluster-num-jobs=100 
		--cluster-options="" """ % os.path.dirname( __file__ ),
    'cmd-sql' : """sqlite3 -header -csv -separator $'\\t' """,
    'cmd-run' : """%s/run.py""" % os.path.dirname( __file__ ),
    }

CONFIG = {}

PROJECT_ROOT = '/ifs/projects'

# local temporary directory to use
TMPDIR = '/scratch'

def configToDictionary( config ):

    p = {}
    for section in config.sections():
        for key,value in config.items( section ):
            v = IOTools.convertValue( value )
            p["%s_%s" % (section,key)] = v
            if section == "general":
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

    Options in **kwargs substitute default
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
    try:
        os.symlink( infile, outfile )
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
    '''return true if file *filename* is empty.'''
    return os.stat(filename)[6]==0

def asList( param ):
    '''return a param as a list'''
    if type(param) not in (types.ListType, types.TupleType):
        try:
            return [x.strip() for x in param.split(",")]
        except AttributeError:
            return [param]
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

def quote( track ):
    '''quote track such that is applicable for a table name.'''
    return re.sub( "[-(),\[\].]", "_", track)

def toTable( outfile ):
    '''convert an outfile (filename) into
    a table name.

    The table name is quoted.
    '''
    assert outfile.endswith( ".load" ) 
    name = os.path.basename( outfile[:-len(".load")] )
    return quote( name )

def getProjectId():
    '''cgat specific method: get the (obfuscated) project id
    based on the current working directory.
    '''
    curdir = os.path.abspath(os.getcwd())
    if not curdir.startswith( PROJECT_ROOT ):
        raise ValueError( "method getProjectId no called within %s" % PROJECT_ROOT )
    prefixes = len(PROJECT_ROOT.split("/"))
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

def load( infile, outfile = None, 
          options = "", transpose = None,
          tablename = None):
    '''straight import from tab separated table.

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

def mergeAndLoad( infiles, outfile, suffix ):
    '''load categorical tables (two columns) into a database.

    The tables are merged and entered row-wise.
    '''
    header = ",".join( [ quote( snip( x, suffix)) for x in infiles] )
    if suffix.endswith(".gz"):
        filenames = " ".join( [ "<( zcat %s | cut -f 1,2 )" % x for x in infiles ] )
    else:
        filenames = " ".join( [ "<( cat %s | cut -f 1,2 )" % x for x in infiles ] )

    tablename = toTable( outfile )

    statement = """python %(scriptsdir)s/combine_tables.py
                      --headers=%(header)s
                      --missing=0
                      --ignore-empty
                   %(filenames)s
                | perl -p -e "s/bin/track/" 
                | python %(scriptsdir)s/table2table.py --transpose
                | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s
            """
    run()

def snip( filename, extension = None):
    '''return prefix of filename.

    If extension is given, make sure that filename has the extension.
    '''
    if extension: 
        assert filename.endswith( extension )        
        return filename[:-len(extension)]

    root, ext = os.path.splitext( filename )
    return root

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
            (options.get("job_queue", global_options.cluster_queue ),
             options.get("job_priority", global_options.cluster_priority ),
             "_" + os.path.basename(options.get("outfile", "ruffus" )),
             options.get("job_options", global_options.cluster_options))

        # keep stdout and stderr separate
        jt.joinFiles=False

        return jt
    
    # run multiple jobs
    if options.get( "statements" ):

        statement_list = []
        for statement in options.get("statements"): 
            options["statement"] = statement
            statement_list.append(buildStatement( **options))
            
        if options.get( "dryrun", False ): return

        # get session for process - only one is permitted
        pid = os.getpid()
        if pid not in global_sessions: 

            L.debug( "creating new drmaa session for pid %i" % pid )
            global_sessions[pid]=drmaa.Session()            
            global_sessions[pid].initialize()

        session = global_sessions[pid]
        
        jt = setupJob( session )
        
        jobids, filenames = [], []
        for statement in statement_list:
            # create job script
            tmpfile = tempfile.NamedTemporaryFile( dir = os.getcwd() , delete = False )
            tmpfile.write( "#!/bin/bash\n" ) #  -l -O expand_aliases\n" )
            tmpfile.write( expandStatement(statement) + "\n" )
            tmpfile.close()

            # build paths
            job_path = os.path.abspath( tmpfile.name )
            stdout_path = job_path + ".stdout" 
            stderr_path = job_path + ".stderr" 

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

            os.unlink( job_path )
            
        session.deleteJobTemplate(jt)

    # run a single parallel job
    elif (options.get( "job_queue" ) or options.get( "to_cluster" )) and not global_options.without_cluster:

        statement = buildStatement( **options )

        if options.get( "dryrun", False ): return

        tmpfile = tempfile.NamedTemporaryFile( dir = os.getcwd() , delete = False )
        tmpfile.write( "#!/bin/bash\n" ) #  -l -O expand_aliases\n" )
        tmpfile.write( expandStatement( statement ) + "\n" )
        tmpfile.close()

        job_path = os.path.abspath( tmpfile.name )
        stdout_path = job_path + ".stdout" 
        stderr_path = job_path + ".stderr" 

        os.chmod( job_path, stat.S_IRWXG | stat.S_IRWXU )

        # get session for process - only one is permitted
        pid = os.getpid()
        if pid not in global_sessions:
            L.debug( "creating new drmaa session for pid %i" % pid )
            global_sessions[pid]=drmaa.Session()            
            global_sessions[pid].initialize()

        session = global_sessions[pid]

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
        os.unlink( job_path )

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
   write a new configuration file pipeline.ini with default values

dump
   write pipeline configuration to stdout

touch
   touch files only, do not run

'''

def main( args = sys.argv ):

    parser = optparse.OptionParser( version = "%prog version: $Id: Pipeline.py 2799 2009-10-22 13:40:13Z andreas $",
                                    usage = USAGE )
    
    parser.add_option( "--pipeline-action", dest="pipeline_action", type="choice",
                       choices=("make", "show", "plot", "dump", "config" ),
                       help="action to take [default=%default]." )

    parser.add_option( "--pipeline-format", dest="pipeline_format", type="choice",
                      choices=( "dot", "jpg", "svg", "ps", "png" ),
                      help="pipeline format [default=%default]." )

    parser.add_option( "-n", "--dry-run", dest="dry_run", action="store_true",
                      help="perform a dry run (do not execute any shell commands) [default=%default]." )

    parser.add_option( "-l", "--local", dest="without_cluster", action="store_true",
                      help="execute all jobs locally [default=%default]." )

    parser.add_option( "-p", "--multiprocess", dest="multiprocess", type="int",
                       help="number of parallel processes to use (different from number of jobs to use for cluster jobs) [default=%default]." ) 

    parser.set_defaults(
        pipeline_action = None,
        pipeline_format = "svg",
        pipeline_target = "full",
        multiprocess = 2,
        logfile = "pipeline.log",
        dry_run = False,
        without_cluster = False,
        )

    (options, args) = E.Start( parser, 
                               add_cluster_options = True )

    global global_options
    global global_args
    global_options, global_args = options, args
    PARAMS["dryrun"] = options.dry_run
    
    version, _ = execute( "hg identify %s" % PARAMS["scriptsdir"] )

    if args: 
        options.pipeline_action = args[0]
        if len(args) > 1:
            options.pipeline_target = args[1]

    if options.pipeline_action in ("make", "show", "svg", "plot", "touch" ):

        try:
            if options.pipeline_action == "make":

                # set up extra file logger
                handler = logging.FileHandler( filename = options.logfile, 
                                               mode = "a" )
                handler.setFormatter( MultiLineFormatter( '%(asctime)s %(levelname)s %(module)s.%(funcName)s.%(lineno)d %(message)s' ) )
                logger = logging.getLogger()
                logger.addHandler( handler )
                
                L.info( E.GetHeader() )
                L.info( "code location: %s" % PARAMS["scriptsdir"] )
                L.info( "code version: %s" % version[:-1] )

                pipeline_run( [ options.pipeline_target ], 
                              multiprocess = options.multiprocess, 
                              logger = logger,
                              verbose = options.loglevel )

                L.info( E.GetFooter() )

            elif options.pipeline_action == "show":
                pipeline_printout( options.stdout, [ options.pipeline_target ], verbose = options.loglevel )

            elif options.pipeline_action == "touch":
                pipeline_run( [ options.pipeline_target ], 
                              touch_files_only = True,
                              verbose = options.loglevel )

            elif options.pipeline_action == "svg":
                pipeline_printout_graph( options.stdout, 
                                         options.pipeline_format,
                                         [ options.pipeline_target ] )
            elif options.pipeline_action == "plot":
                outf, filename = tempfile.mkstemp()
                pipeline_printout_graph( os.fdopen(outf,"w"),
                                         options.pipeline_format,
                                         [ options.pipeline_target ] )
                execute( "inkscape %s" % filename ) 
                os.unlink( filename )

        except ruffus_exceptions.RethrownJobError, value:
            E.error("some tasks resulted in errors - error messages follow below" )
            E.error( value )
            raise

    elif options.pipeline_action == "dump":
        # convert to normal dictionary (not defaultdict) for parsing purposes
        print "dump = %s" % str(dict(PARAMS))
    elif options.pipeline_action == "config":
        if os.path.exists("pipeline.ini"):
            raise ValueError( "file `pipeline.ini` already exists" )
        f = sys._getframe(1)
        caller = inspect.getargvalues(f).locals["__file__"]
        configfile = os.path.splitext(caller)[0] + ".ini" 
        if not os.path.exists( configfile ):
            raise ValueError( "default config file `%s` not found"  % configfile )
        shutil.copyfile( configfile, "pipeline.ini" )
        L.info( "created new configuration file `pipeline.ini` " )
    else:
        raise ValueError("unknown pipeline action %s" % options.pipeline_action )

    E.Stop()

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
    
    dirname = os.path.dirname( __file__ )
    pipeline = os.path.join( dirname, pipeline )
    assert os.path.exists( pipeline ), "can't find pipeline source %s" % pipeline
    assert os.path.exists( workingdir ), "can't find working dir %s" % workingdir
    
    process = subprocess.Popen(  "python %s -v 0 dump" % pipeline,
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
    if xvfb_command: xvfb_command+= " -a "

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

def publish_report( prefix = "", patterns = [], project_id = None):
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

    def _link( src, dest ):
        '''create links. 

        Only link to existing targets.
        '''
        dest = os.path.abspath( os.path.join( PARAMS["web_dir"], dest ) )
        if os.path.exists( dest ):
            os.remove(dest)
        if os.path.exists( src ):
            os.symlink( os.path.abspath(src), dest )

    def _copy( src, dest ):
        dest = os.path.abspath( os.path.join( PARAMS["web_dir"], dest ) )
        if os.path.exists( dest ): shutil.rmtree( dest )
        shutil.copytree( os.path.abspath(src), dest ) 

    # publish export dir via symlinking
    _link( src_export, dest_export )

    # publish web pages by copying
    _copy( os.path.abspath("report/html"), dest_report ) 

    # substitute links to export
    _patterns = [ (re.compile( src_export ), 
                   "http://www.cgat.org/downloads/%(project_id)s/%(dest_export)s" % locals() ), 
                  ]
    
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

if __name__ == "__main__":

    main()

    #parser = optparse.OptionParser( version = "%prog version: $Id: Pipeline.py 2799 2009-10-22 13:40:13Z andreas $")

    #(options, args) = E.Start( parser, add_cluster_options = True )

    #global global_options
    #global global_args
    #global_options, global_args = options, args

    #L.info("in main")

    #run( **{"statement" : "printenv > test.out", "job_queue" : "server_jobs.q", "job_priority" : -10 } )
    #run( **{"statement" : "printenv > test2.out", "job_queue" : "server_jobs.q", "job_priority" : -10 } )
    #run( **{"statement" : "printenv > test3.out", "job_queue" : "server_jobs.q", "job_priority" : -10 } )
    
