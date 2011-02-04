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

Code
----

'''
import os, sys, re, subprocess, optparse, stat, tempfile, time, random, inspect, types, glob, shutil, logging
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

PARAMS= { 
    'scriptsdir' : "/ifs/devel/cgat",
    'toolsdir' : "/ifs/devel/cgat",
    'cmd-farm' : """/ifs/devel/cgat/farm.py 
                --method=drmaa 
                --cluster-priority=-10 
		--cluster-queue=all.q 
		--cluster-num-jobs=100 
		--cluster-options="" """,
    'cmd-sql' : """sqlite3 -header -csv -separator $'\\t' """,
           }

CONFIG = {}

def configToDictionary( config ):

    p = {}

    for section in config.sections():
        for key,value in config.items( section ):
            v = IOTools.convertValue( value )
            p["%s_%s" % (section,key)] = v
            if section == "general":
                p["%s" % (key)] = v

    return p

def getParameters( filenames = ["pipeline.ini",] ):
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
    
    '''

    global CONFIG

    CONFIG = ConfigParser.ConfigParser()
    CONFIG.read( filenames )

    p = configToDictionary( CONFIG )
    PARAMS.update( p )

    return p

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

def touch( filename, times = None ):
    '''update/create a sentinel file.'''
    fhandle = file(filename, 'a')
    try:
        os.utime(filename, times)
    finally:
        fhandle.close()
    
def getTempFilename( dir = None ):
    '''get a temporary filename.

    The caller needs to delete it.
    '''
    tmpfile = tempfile.NamedTemporaryFile(dir=dir,delete=False)
    tmpfile.close()
    return tmpfile.name

def getTempFile( dir = None ):
    '''get a temporary file.

    The caller needs to delete it.
    '''
    return tempfile.NamedTemporaryFile(dir=dir,delete=False)

def getTempDir( dir = None ):
    return tempfile.mkdtemp( dir=dir )

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

def asList( param ):
    '''return a param as a list'''
    if type(param) not in (types.ListType, types.TupleType):
        return [param,]
    else: return param

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
    name = outfile[:-len(".load")]
    return quote( name )

def load( infile, outfile, options = "" ):
    '''straight import from tab separated table.

    The table name is given by outfile without the
    ".load" suffix.
    '''

    tablename = toTable( outfile )

    if infile.endswith(".gz"): cat = "zcat"
    else: cat = "cat"

    statement = '''%(cat)s %(infile)s
    | csv2db.py %(csv2db_options)s 
              %(options)s 
              --table=%(tablename)s 
    > %(outfile)s
    '''

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
    '''

_exec_suffix = "; detect_pipe_error"

def buildStatement( **kwargs ):
    '''build statement from kwargs.'''

    # the actual statement
    try:
        statement = kwargs.get("statement") % dict( PARAMS.items() + kwargs.items() )
    except KeyError, msg:
        raise KeyError( "Error when creating command: could not find %s in dictionaries" % msg)

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

    if not kwargs: kwargs = getCallerLocals()

    # run multiple jobs
    if kwargs.get( "statements" ):

        statement_list = []
        for statement in kwargs.get("statements"): 
            kwargs["statement"] = statement
            statement_list.append(buildStatement( **kwargs))
            
        if kwargs.get( "dryrun", False ): return

        # get session for process - only one is permitted
        pid = os.getpid()
        if pid not in global_sessions: 

            L.debug( "creating new drmaa session for pid %i" % pid )
            global_sessions[pid]=drmaa.Session()            
            global_sessions[pid].initialize()

        session = global_sessions[pid]

        jt = session.createJobTemplate()
        jt.workingDirectory = os.getcwd()
        jt.jobEnvironment = { 'BASH_ENV' : '~/.bashrc' }
        jt.args = []
        jt.nativeSpecification = "-q %s -p %i -N %s %s" % \
            (kwargs.get("job_queue", global_options.cluster_queue ),
             kwargs.get("job_priority", global_options.cluster_priority ),
             os.path.basename(kwargs.get("outfile", "ruffus" )),
             kwargs.get("job_options", global_options.cluster_options))
        
        # keep stdout and stderr separate
        jt.joinFiles=False

        jobids, filenames = [], []
        for statement in statement_list:
            # create job scrip
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
                raise PipelineError( "Child was terminated by signal %i: \nThe stderr was: \n%s\n%s\n" % \
                                         (retval.exitStatus, 
                                          "".join( stderr),
                                          statement ) )

            os.unlink( job_path )
            
        session.deleteJobTemplate(jt)

    # run a single parallel job
    elif (kwargs.get( "job_queue" ) or kwargs.get( "to_cluster" )) and not global_options.without_cluster:

        statement = buildStatement( **kwargs )

        if kwargs.get( "dryrun", False ): return

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

        jt = session.createJobTemplate()
        jt.workingDirectory = os.getcwd()
        jt.remoteCommand = job_path
        jt.jobEnvironment = { 'BASH_ENV' : '~/.bashrc' }
        jt.args = []
        jt.nativeSpecification = "-q %s -p %i -N %s %s" % \
            (kwargs.get("job_queue", global_options.cluster_queue ),
             kwargs.get("job_priority", global_options.cluster_priority ),
             os.path.basename(kwargs.get("outfile", "ruffus" )),
             kwargs.get("job_options", global_options.cluster_options))
        
        # keep stdout and stderr separate
        jt.joinFiles=False
        # later: allow redirection of stdout and stderr to files; can even be across hosts?
        jt.outputPath=":"+ stdout_path
        jt.errorPath=":" + stderr_path

        if "job_array" in kwargs and kwargs["job_array"] != None:
            # run an array job
            start, end, increment = kwargs.get("job_array" )
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

        if "job_array" not in kwargs:
            if retval and retval.exitStatus != 0:
                raise PipelineError( "Child was terminated by signal %i: \nThe stderr was: \n%s\n%s\n" % \
                                         (retval.exitStatus, 
                                          "".join( stderr), statement))
            
        session.deleteJobTemplate(jt)
        os.unlink( job_path )

    else:
        statement = buildStatement( **kwargs )

        if kwargs.get( "dryrun", False ): return
 
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
            raise PipelineError( "Child was terminated by signal %i: \nThe stderr was: \n%s\n%s\n" % (-process.returncode, stderr, statement ))


class MultiLineFormatter(logging.Formatter):
    '''logfile formatter: add identation for multi-line entries.'''

    def format(self, record):
        str = logging.Formatter.format(self, record)
        header, footer = str.split(record.message)
        str = str.replace('\n', '\n' + ' '*len(header))
        return str

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

    if options.pipeline_action in ("make", "show", "svg", "plot"):

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
            print "re-raising exception"
            raise

    elif options.pipeline_action == "dump":
        print "dump = %s" % str(PARAMS)
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

if __name__ == "__main__":

    parser = optparse.OptionParser( version = "%prog version: $Id: Pipeline.py 2799 2009-10-22 13:40:13Z andreas $")

    (options, args) = E.Start( parser, add_cluster_options = True )

    #global global_options
    #global global_args
    global_options, global_args = options, args

    run( **{"statement" : "printenv > test.out", "job_queue" : "server_jobs.q", "job_priority" : -10 } )
    run( **{"statement" : "printenv > test2.out", "job_queue" : "server_jobs.q", "job_priority" : -10 } )
    run( **{"statement" : "printenv > test3.out", "job_queue" : "server_jobs.q", "job_priority" : -10 } )
    
