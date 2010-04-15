import os, sys, re, subprocess, optparse, stat, tempfile, time, random, inspect, types
import ConfigParser

import drmaa
from ruffus import *
import Experiment as E
import IOTools

global_options, global_args = None, None

# dictionary mapping process ids (pids) to drmaa sessions
global_sessions = {}

class PipelineError( Exception ): pass

PARAMS= { 
    'scriptsdir' : "/home/andreas/cgat",
    'toolsdir' : "/home/andreas/cgat",
    'cmd-farm' : """farm.py 
                --method=drmaa 
                --cluster-priority=-10 
		--cluster-queue=medium_jobs.q 
		--cluster-num-jobs=100 
		--cluster-options="" """,
    'cmd-sql' : """sqlite3 -header -csv -separator $'\\t' """,
           }

def getParameters( filename = "pipeline.ini" ):
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
    p = {}
    
    config = ConfigParser.ConfigParser()
    config.readfp(open(filename),"r")

    for section in config.sections():
        for key,value in config.items( section ):
            v = IOTools.convertValue( value )
            p["%s_%s" % (section,key)] = v
            if section == "general":
                p["%s" % (key)] = v


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
    E.info( message )

def warning( message ):
    E.warning( message )

def warn( message ):
    E.warning( message )

def debug( message ):
    E.debug( message )

def error( message ):
    E.error( message )
        
def critical( message):
    E.critical( message )

def asList( param ):
    '''return a param as a list'''
    if type(param) not in (types.ListType, types.TupleType):
        return [param,]
    else: return param

def getCallerLocals(decorators=0):
    '''returns locals of caller using frame.

    optional pass number of decorators
    
    from http://pylab.blogspot.com/2009/02/python-accessing-caller-locals-from.html
    '''
    f = sys._getframe(2+decorators)
    args = inspect.getargvalues(f)
    return args[3]

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

    # prefix to detect errors within pipes
    prefix = '''detect_pipe_error_helper() 
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
    suffix = "; detect_pipe_error"

    if not kwargs: kwargs = getCallerLocals()

    # the actual statement
    statement = kwargs.get("statement") % dict( PARAMS.items() + kwargs.items() )

    # add bash as prefix to allow advanced shell syntax like 'wc -l <( gunzip < x.gz)'
    # executable option to call() does not work. Note that there will be an extra
    # indirection.
    statement = " ".join( re.sub( "\t+", " ", statement).split( "\n" ) )

    E.debug( "running statement:\n%s" % statement )

    if kwargs.get( "dryrun", False ):
        return

    statement = " ".join( (prefix, statement, suffix) )

    if (kwargs.get( "job_queue" ) or kwargs.get( "to_cluster" )) and not global_options.without_cluster:

        tmpfile = tempfile.NamedTemporaryFile( dir = os.getcwd() , delete = False )
        tmpfile.write( "#!/bin/bash\n" ) #  -l -O expand_aliases\n" )
        tmpfile.write( statement + "\n" )
        tmpfile.close()

        job_path = os.path.abspath( tmpfile.name )
        stdout_path = job_path + ".stdout" 
        stderr_path = job_path + ".stderr" 

        os.chmod( job_path, stat.S_IRWXG | stat.S_IRWXU )

        # get session for process - only one is permitted
        pid = os.getpid()
        if pid not in global_sessions:
            global_sessions[pid]=drmaa.Session()            
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

        jobid = session.runJob(jt)
        retval = session.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)

        stdout = open( stdout_path, "r" ).readlines()
        stderr = open( stderr_path, "r" ).readlines()

        if retval.exitStatus != 0:
            raise PipelineError( "Child was terminated by signal %i: \nThe stderr was: \n%s\n%s\n" % (retval.exitStatus, "".join( stderr), statement))
        
        session.deleteJobTemplate(jt)
        os.unlink( job_path )
        os.unlink( stdout_path )
        os.unlink( stderr_path )

    else:
        if "<(" in statement:
            if "'" in statement: raise ValueError( "advanced bash syntax combined with single quotes" )
            statement = """/bin/bash -c '%s'""" % statement

        process = subprocess.Popen(  statement,
                                     cwd = os.getcwd(), 
                                     shell = True,
                                     stdin = subprocess.PIPE,
                                     stdout = subprocess.PIPE,
                                     stderr = subprocess.PIPE )

        # process.stdin.close()
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            raise PipelineError( "Child was terminated by signal %i: \nThe stderr was: \n%s\n%s\n" % (-process.returncode, stderr, statement ))

def main( args = sys.argv ):

    parser = optparse.OptionParser( version = "%prog version: $Id: Pipeline.py 2799 2009-10-22 13:40:13Z andreas $")
    
    parser.add_option( "--pipeline-action", dest="pipeline_action", type="choice",
                       choices=("make", "show", "plot" ),
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
        dry_run = False,
        without_cluster = False,
        )

    (options, args) = E.Start( parser, add_cluster_options = True )

    global global_options
    global global_args
    global_options, global_args = options, args
    PARAMS["dryrun"] = options.dry_run

    if args: 
        options.pipeline_action = args[0]
        if len(args) > 1:
            options.pipeline_target = args[1]


    try:
        if options.pipeline_action == "make":
            pipeline_run( [ options.pipeline_target ], multiprocess = options.multiprocess, verbose = options.loglevel )
        elif options.pipeline_action == "show":
            pipeline_printout( options.stdout, [ options.pipeline_target ], verbose = options.loglevel )
        elif options.pipeline_action == "plot":
            pipeline_printout_graph( options.stdout, 
                                     options.pipeline_format,
                                     [ options.pipeline_target ] )
    except ruffus_exceptions.RethrownJobError, value:
        print "re-raising exception"
        raise

    E.Stop()

if __name__ == "__main__":
    parser = optparse.OptionParser( version = "%prog version: $Id: Pipeline.py 2799 2009-10-22 13:40:13Z andreas $")
    (options, args) = E.Start( parser, add_cluster_options = True )
    #global global_options
    #global global_args
    global_options, global_args = options, args

    run( **{"statement" : "printenv > test.out", "job_queue" : "server_jobs.q", "job_priority" : -10 } )
    run( **{"statement" : "printenv > test2.out", "job_queue" : "server_jobs.q", "job_priority" : -10 } )
    run( **{"statement" : "printenv > test3.out", "job_queue" : "server_jobs.q", "job_priority" : -10 } )
