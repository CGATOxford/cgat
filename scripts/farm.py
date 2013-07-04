#!/bin/env python
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
farm.py - execute a cmd on the cluster
======================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   cat go | farm.py --split-at-colum=1 perl -p -e "s/GO/gaga/"

Type::

   python farm.py --help

for command line help.

Documentation
-------------

The input on stdin is split for embarrasingly parallel jobs.
The --split-at-.. options describe how standard input is to
be split. A temporary directory is created in the current
directory. This directory has be visible on the cluster nodes
and accessible under the same name.

The output is written to stdout. Results are returned in the
same order as they are submitted.

On error, error messages are echoed and nothing is returned.
The temporary directory is not deleted to allow manual recovery.

Examples
--------

The following command will split the file "go" at the first column 
and execute the command perl -p -e "s/GO/gaga/"::

   cat go | farm.py --split-at-colum=1 perl -p -e "s/GO/gaga/"

The following command will split a fasta file at each entry and
compute an approximate sequence length::

   cat genome.fasta | farm.py --split-at-regex="^>(\S+)" "wc -c"

The following command will split a fasta file at every 10 sequences::

   cat genome.fasta | farm.py --split-at-regex="^>(\S+)" --chunksize=10 "wc -c"

.. todo::

   implement continuation of jobs
   implement better error messages
   use sge array jobs for job control

Code
----

'''

import os
import sys
import re
import string
import optparse
import time
import glob
import subprocess
import tempfile
import shutil
import stat

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import threadpool
import multiprocessing

try:
    import drmaa
    HASDRMAA = True
except ImportError:
    HASDRMAA = False


def chunk_iterator_lines( infile, args, prefix, use_header = False ):
    """split by lines."""

    chunk_size = args[0]
    n = 0
    filename = "%s/%010i.in" % (prefix, n)
    outfile = IOTools.openFile( filename, "w" )
    header = None

    for line in infile:
        if line[0] == "#": continue

        if not header and n == 0 and use_header:
            header = line
            outfile.write( header )
            continue

        n += 1

        if n % chunk_size == 0:
            outfile.close()
            yield filename
            filename = "%s/%010i.in" % (prefix, n)
            outfile = IOTools.openFile( filename, "w" )
            if header: outfile.write( header )

        outfile.write( line )
    outfile.close()
    yield filename

def chunk_iterator_column( infile, args, prefix, use_header = False ):
    """split at column. 

    The table need not be sorted by this column.
    If num_files is given, files will randomly created
    and tags according to column randomly assigned.


    """

    column, max_files = args
    files = IOTools.FilePool()
    header = False

    if max_files:
        map_tag2file = {}

    for line in infile:
        if line[0] == "#": continue

        if not header and use_header: 
            files.setHeader( line )
            header = True
            continue

        key = line[:-1].split("\t")[column]
        if max_files:
            if key in map_tag2file:
                key= map_tag2file[key]
            else:
                n= "%010i" % (len(map_tag2file) % max_files)
                map_tag2file[key] = n
                key = n

        files.write( "%s/%s.in" % (prefix, key), line )

    for filename, count in files.items():
        E.info( "created file %s with %i items" % (filename, count ) )
        yield filename

def chunk_iterator_regex_group( infile, args, prefix, use_header = False ):
    """group by regular expression is true.

    Entries need to be consecutive.
    """

    rex = args[0]
    column = args[1]
    chunk_size = args[2]
    last = None
    header = None
    n = chunk_size
    outfile = None

    for line in infile:

        if line[0] == "#": continue

        if not header and use_header: 
            header = line
            continue

        try:
            this = rex.search(line[:-1]).groups()[0]
        except IndexError:
            if outfile: outfile.write(line)
            continue
        except AttributeError:
            if outfile: outfile.write(line)
            continue

        if last != this and n >= chunk_size:
            if last:
                outfile.close()
                yield filename

            last = this

            filename = "%s/%s.in" % (prefix, this)
            outfile = IOTools.openFile( filename, "w" )
            if header: outfile.write( header )
            n = 0

        outfile.write( line )
        n += 1

    if outfile:
        outfile.close()
        yield filename

def chunk_iterator_regex_split( infile, args, prefix, use_header = False ):
    """split where regular expression is true.
    """

    rex = args[0]
    column = args[1]
    chunk_size = args[2]
    max_lines = args[3]
    last = None
    header = None

    nlines = 0
    n = 0
    filename = "%s/%010i.in" % (prefix, n)
    outfile = IOTools.openFile( filename, "w" )

    for line in infile:

        if line[0] == "#": continue

        if rex.search(line[:-1]):
            if n > 0 and (n % chunk_size == 0 or (max_lines and nlines > max_lines)):
                outfile.close()
                yield filename
                filename = "%s/%010i.in" % (prefix, n)
                outfile = IOTools.openFile( filename, "w" )
                nlines = 0

            n += 1
                
        outfile.write( line )
        nlines += 1

    outfile.close()
    yield filename

def chunk_iterator_psl_overlap( infile, args, prefix, use_header = False ):
    """iterate over overlapping entries in a psl file."""

    iterator = Blat.BlatIterator( sys.stdin )

    processed_contigs = set()

    merge_distance=args[0]

    while 1:

        match = iterator.next()
        
        if match == None: break
        
        if match.mSbjctId != last_sbjct_id or match.mSbjctFrom >= (sbjct_end + merge_distance):
            if last_sbjct_id:
                outfile.close()
                yield filename

            if last_sbjct_id != match.mSbjctId and match.mSbjctId in processed_contigs:
                raise ValueError("input not sorted correctly (contig,start): already encountered %s\n%s" % (match.mSbjctId, str(match)))

            last_sbjct_id = match.mSbjctId
            processed_contigs.add( last_sbjct_id )

            sbjct_start = match.mSbjctFrom
            sbjct_end = match.mSbjctTo
            
        if match.mSbjctFrom < sbjct_start:
            raise ValueError("input not sorted correctly (contig,start): %i < %i\n%s" % (match.mSbjctFrom, sbjct_start, str(match)))

        sbjct_end = max( match.mSbjctTo, sbjct_end )
        outfile.write( str(match)+"\n" )

    if outfile:
        outfile.close()
        yield filename

class MapperGlobal:
    def __init__(self, pattern="%06i"):
        self.mMap = {}
        assert "%" in pattern, "please supply a pattern"
        self.mPattern = pattern

    def __call__(self, fn, id ):
        if id not in self.mMap:
            self.mMap[id] = self.mPattern % (len(self.mMap) + 1)
        return self.mMap[id]

class MapperLocal:
    def __init__(self, pattern="%06i"):
        self.mMap = {}
        assert "%" in pattern, "please supply a pattern"
        self.mPattern = pattern

    def __call__(self, fn, id ):
        key = "%s-%s" % (fn, id)
        if key not in self.mMap:
            self.mMap[key] = self.mPattern % (len(self.mMap) + 1)
        return self.mMap[key]

class MapperEmpty:
    def __init__(self):
        pass
    def __call__(self, fn, id ):
        return id

class ResultBuilder:
    """the default result builder for table formatted output."""

    def __init__( self, mapper = None, field_index = None, field_name = None):
        self.mMapper = mapper
        self.mFieldIndex = field_index
        self.mFieldName = field_name
        self.mHeader = None

    def parseHeader( self, infile, outfile, options ):
        """parse header in infile."""
        # skip comments until header
        while 1:
            l = infile.readline() 
            if not l or l[0] != "#": break
            options.stdlog.write(l)

        # Print only the first header and check if
        # all the headers are the same.
        if self.mHeader:
            if self.mHeader != l:
                raise ValueError( "inconsistent header in file %s\ngot=%s\nexpected=%s" % (infile, l, self.mHeader) )
        else:
            outfile.write( l )
            self.mHeader = l

            if self.mFieldIndex == None and self.mFieldName:
                try:
                    self.mFieldIndex = self.mHeader.split("\t").index( self.mFieldName )
                except ValueError:
                    if options.loglevel >= 1:
                        options.stdlog.write("# Warning: no mapping, can not find field %s in %s" % (self.mFieldName, self.mHeader) )
                    self.mFieldName = None
            
                if options.loglevel >= 2:
                    options.stdlog.write("# substituting field: %s, %s\n" % (self.mFieldName, self.mFieldIndex))

    def __call__( self, filenames, outfile, options ):

        for fi, fn in filenames:
            E.debug( "# merging %s" % fn )
            infile = IOTools.openFile( fn, "r" )

            if options.output_header: 
                self.parseHeader( infile, outfile, options )

            for l in infile:
                if l[0] == "#":
                    options.stdlog.write(l)
                else:
                    if self.mFieldIndex != None:
                        data = l[:-1].split("\t")
                        try:
                            data[self.mFieldIndex] = self.mMapper( fi, data[self.mFieldIndex] )
                        except IndexError:
                            raise IndexError( "can not find field %i in %s" % (self.mFieldIndex, l ) )
                        l = "\t".join( data ) + "\n"
                        
                    outfile.write(l)
            infile.close()

class ResultBuilderPSL ( ResultBuilder ):
    """Result builder for psl tables. Here, column 9, the query id, is substituted."""

    def __init__( self, *args, **kwargs):
        ResultBuilder.__init__(self, *args, **kwargs)
        self.mFieldIndex = 9
        self.mFirst = True

    def parseHeader( self, infile, outfile, options ):
        """parse header in infile."""
        # skip comments until header
        while 1:
            l = infile.readline() 
            if not l or l[0] != "#": break
            options.stdlog.write(l)

        if l.startswith( "psLayout version 3"):
            if self.mFirst:
                outfile.write( l )
                for x in range(0,4):
                    l = infile.readline()
                    outfile.write( l )
                self.mFirst = False
            else:
                for x in range(0,4):
                    l = infile.readline()

class ResultBuilderFasta( ResultBuilder) :

    def __init__( self, *args, **kwargs ):
        ResultBuilder.__init__(self, *args, **kwargs)
    
    def __call__( self, filenames, outfile, options ):
        for fi, fn in filenames:
            infile = IOtools.openFile( fn, "r" )
            for l in infile:
                if l[0] == "#": 
                    options.stdlog.write( l )
                    continue
                elif l[0] == ">":
                    x = re.search(">(\S+)", l[:-1])
                    id = self.mMapper( fi, x.groups()[0] )
                    l = ">%s%s" % (id, l[x.end(0):])
                outfile.write(l)
            infile.close()

class ResultBuilderBinary( ResultBuilder) :
    '''simply concatenate output files (without any parsing).'''

    def __init__( self, *args, **kwargs ):
        ResultBuilder.__init__(self, *args, **kwargs)
    
    def __call__( self, filenames, outfile, options ):
        for fi, fn in filenames:
            shutil.copyfileobj( IOTools.openFile(fn, "r"), outfile )

class ResultBuilderCopies( ResultBuilder) :
    '''create indexed copiers.'''

    def __init__( self, *args, **kwargs ):
        ResultBuilder.__init__(self, *args, **kwargs)
    
    def __call__( self, filenames, outfile, options ):
        idx = 0
        base, ext = os.path.splitext( outfile.name )
        for fi, fn in filenames:
            idx += 1
            shutil.copyfile( fn, base + ".%i" % idx + ext )

class ResultBuilderLog( ResultBuilder) :
    """processor for log files."""

    def __init__( self, *args, **kwargs ):
        ResultBuilder.__init__(self, *args, **kwargs)
    
    def __call__( self, filenames, outfile, options ):
        for fi, fn in filenames:
            infile = IOTools.openFile( fn, "r" )
            outfile.write("######### logging output for %s ###################\n" % fi)
            for l in infile:
                outfile.write(l)
            infile.close()

##--------------------------------------------------------------------
def runCommand( data ):

    filename, cmd, options, tmpdir, subdirs = data

    if subdirs:
        outdir = "%s.dir/" % (filename)
        os.mkdir( outdir )
        cmd = re.sub( "%DIR%", outdir, cmd )

    x = re.search( "'--log=(\S+)'", cmd ) or re.search( "'--L\s+(\S+)'", cmd )
    if x:
        logfile = filename + ".log" 
        cmd = cmd[:x.start()] + "--log=%s" % logfile + cmd[x.end():]
    else:
        logfile = filename + ".out"

    # working directory - needs to be the one from which the
    # the script is called to resolve input files.
    cwd = os.getcwd()

    if "<(" in cmd or "|" in cmd:
        if "'" in cmd: raise ValueError( "advanced bash syntax `<()` combined with single quotes" )
        cmd = """/bin/bash -c '%s'""" % cmd

    if "|" in cmd:
        if r"\|" not in cmd:
            E.warn( "pipes (`|`) within command need to be escaped, otherwise jobs run on submit host" )

    c = '%s -v "BASH_ENV=%s" -q %s -p %i %s %s' % (options.cluster_cmd, 
                                                   options.bashrc,
                                                   options.cluster_queue,
                                                   options.cluster_priority,
                                                   options.cluster_options,
                                                   cmd )

    iteration = 0

    while 1:

        iteration += 1
        if iteration > 1:
            E.info("%s: re-submitting command (repeat=%i): %s" % (filename, iteration, c ))
        else:
            E.info("%s: submitting command: %s" % (filename, c))

        infile = IOTools.openFile( filename, "r" )
        outfile = IOTools.openFile( filename + ".out", "w")
        errfile = IOTools.openFile( filename + ".err", "a")

        retcode = subprocess.call( c,
                                   shell = True,
                                   stdin = infile,
                                   stdout = outfile,
                                   stderr = errfile,
                                   cwd = cwd,
                                   close_fds = True)                              

        infile.close()
        outfile.close()
        errfile.close()

        if hasFinished( retcode, filename, options.output_tag, logfile ):
            break

        if iteration > options.resubmit:
            E.warn( "%s: giving up executing command: retcode=%i" % (filename, retcode) )
            break

        E.warn( "%s: error while executing command: retcode=%i" % (filename, retcode) )


    return (retcode, filename, cmd, logfile, iteration)

##--------------------------------------------------------------------
def hasFinished( retcode, filename, output_tag, logfile ):
    """check if a run has finished."""

    E.info( "checking status of job %s with returncode %i" % (filename, retcode) )
    if retcode != 0: 
        try:
            if not output_tag or not re.search( output_tag, IOTools.getLastLine( logfile )):
                return False
        except IOError:
            E.warn( "could not read output_tag from files %s" % (logfile) )
            return False
    return True


##--------------------------------------------------------------------
def runDRMAA( data, environment ):
    '''run jobs in data using drmaa to connect to the cluster.'''
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
    
    # working directory - needs to be the one from which the
    # the script is called to resolve input files.
    cwd = os.getcwd()

    session = drmaa.Session()            
    session.initialize()
    jt = session.createJobTemplate()
    
    jobids = []
    kwargs = {}

    for filename, cmd, options, tmpdir, subdirs in data:

        from_stdin, to_stdout = True, True

        if subdirs:
            outdir = "%s.dir/" % (filename)
            os.mkdir( outdir )
            cmd = re.sub( "%DIR%", outdir, cmd )

        x = re.search( "'--log=(\S+)'", cmd ) or re.search( "'--L\s+(\S+)'", cmd )
        if x:
            logfile = filename + ".log" 
            cmd = cmd[:x.start()] + "--log=%s" % logfile + cmd[x.end():]
        else:
            logfile = filename + ".out"

        if "%STDIN%" in cmd:
            cmd = re.sub("%STDIN%", filename, cmd )
            from_stdin = False

        if "%STDOUT%" in cmd:
            cmd = re.sub("%STDOUT%", filename + ".out", cmd )
            to_stdout = False

        cmd = " ".join( re.sub( "\t+", " ", cmd).split( "\n" ) )
        E.debug( "running statement:\n%s" % cmd )

        tmpfile = tempfile.NamedTemporaryFile( dir = os.getcwd() , delete = False )
        tmpfile.write( "#!/bin/bash\n" ) #  -l -O expand_aliases\n" )
        tmpfile.write( " ".join( (prefix, cmd, suffix)) + "\n" )
        tmpfile.close()

        job_path = os.path.abspath( tmpfile.name )

        os.chmod( job_path, stat.S_IRWXG | stat.S_IRWXU )

        # get session for process - only one is permitted
        jt.workingDirectory = os.getcwd()
        jt.remoteCommand = job_path
        e = { 'BASH_ENV' : options.bashrc }
        if environment:
            for en in environment:
                try:
                    e[en] = os.environ[en]
                except KeyError:
                    raise KeyError("could not export environment variable '%s'" % en )
        jt.jobEnvironment = e

        jt.args = []
        jt.nativeSpecification = "-V -q %s -p %i -N %s %s" % \
            (kwargs.get("job_queue", options.cluster_queue ),
             kwargs.get("job_priority", options.cluster_priority ),
             os.path.basename(kwargs.get("outfile", "farm.py" )),
             kwargs.get("job_options", options.cluster_options))
        
        # keep stdout and stderr separate
        jt.joinFiles=False

        # use stdin for data
        if from_stdin: jt.inputPath=":" + filename

        # later: allow redirection of stdout and stderr to files
        # could this even be across hosts?
        if to_stdout: 
            jt.outputPath=":"+ filename + ".out"
        else:
            jt.outputPath=":"+ filename + ".stdout"

        jt.errorPath=":" + filename + ".err"


        jobid = session.runJob(jt)
        jobids.append( (jobid, job_path, filename, cmd, logfile) )

    E.debug( "%i jobs have been submitted" % len(jobids) )

    results = []

    for jobid, job_path, filename, cmd, logfile in jobids:

        try:
            retval = session.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)
        except Exception, msg:
            # ignore message 24 in PBS
            # code 24: drmaa: Job finished but resource usage information and/or termination status could not be provided.":
            if not msg.message.startswith("code 24"): raise
            retval = None

        if retval and retval.exitStatus != 0:
            raise OSError( "Child was terminated by signal %i: \n%s\n" % \
                               (retval.exitStatus, cmd))
        
        results.append( (retval, filename, cmd, logfile, 1 ) )

        os.unlink( job_path )

    session.deleteJobTemplate(jt)
    session.exit()

##--------------------------------------------------------------------
def getOptionParser():
    """create parser and add options."""

    parser = E.OptionParser( version = "%prog version: $Id: farm.py 2782 2009-09-10 11:40:29Z andreas $", 
                                    usage = globals()["__doc__"] )


    parser.add_option( "--split-at-lines", dest="split_at_lines", type="int",
                       help = "split jobs according to line number [default=%default]." )

    parser.add_option( "--split-at-column", dest="split_at_column", type="int",
                       help = "split jobs according to column. Columns start at number 1 and the input should be sorted by this column [default=%default]." )

    parser.add_option( "--group-by-regex", dest="group_by_regex", type="string",
                       help = "group jobs according to a regular expression [default=%default]." )

    parser.add_option( "--split-at-regex", dest="split_at_regex", type="string",
                       help = "split jobs according to a regular expression [default=%default]." )

    parser.add_option( "--split-at-tag", dest="split_at_tag", type="int",
                       help = "split a file at a tag [default=%default]." )

    parser.add_option( "--chunksize", dest="chunksize", type="int",
                       help = "when splitting at regex or tag, aggregate x entries [default=%default]." )

    parser.add_option( "--debug", dest="debug", action="store_true",
                       help = "debug mode. Do not delete temporary file [default=%default]." )

    parser.add_option( "--dry-run", dest="dry_run", action="store_true",
                       help = "dry run. Do not split input and simply forward stdin to stdout - " 
                              "useful for debugging the command [default=%default]." )

    parser.add_option( "--input-header", dest="input_header", action="store_true",
                       help = "The input stream contains a table header. " 
                              "This header is replicated for each job [default=%default]." )

    parser.add_option( "--output-header", dest="output_header", action="store_true",
                       help = "The output jobs contain a table header. " 
                              "The header is removed for each job except for the first [default=%default]." )

    parser.add_option( "--output-tag", dest="output_tag", type="string",
                       help = "The output jobs contain a tag in the last line denoting job completion. If the unix return value denotes an error, the presence of this tag is checked [default=%default]." )

    parser.add_option( "--subdirs", dest="subdirs", action="store_true",
                       help = "Run within separate subdirs for jobs. This permits multiple output streams. Use a placeholder %DIR% if you supply the ouput pattern as a command line option [default=%default]." )

    parser.add_option( "-T", "--tmpdir", dest="tmpdir", type="string",
                       help = "Temporary directory to be used. Default is the current directory [default=%default]." )
    
    parser.add_option( "--max-files", dest="max_files", type="int",
                       help = "create at most x files [default=%default]." )

    parser.add_option( "--max-lines", dest="max_lines", type="int",
                       help = "in addition to splitting into chunksize, also split if more than max-lines is reached [default=%default]." )

    parser.add_option( "--renumber", dest="renumber", type="string",
                       help = "renumber ids consecutively, supply a pattern [default=%default]." )

    parser.add_option( "--renumber-column", dest="renumber_column", type="string", action="append",
                       help = "specify column to renumber. The format is regex:column, for example csv:1 or csv:id [default=%default]." )
    
    parser.add_option( "-r", "--reduce", dest="reduce", type="string", action="append",
                       help = "Add reduce functions for specific files. The format is file:reducer. The default reducer is 'table' for all files [default=%default]." )
    
    parser.add_option( "-m", "--map", dest="map", type="string", action="append",
                       help = "Map specific columns in tables. The format is file:column:pattern, for example .table:1:%06i [default=%default]." )

    parser.add_option( "--resume", dest="resume", type ="string",
                       help = "resume aborted run from files in dir [%default]")

    parser.add_option( "--collect", dest="collect", type ="string",
                       help = "collect files in dir and process as normally [%default]")

    parser.add_option( "--binary", dest="binary", action="store_true",
                       help = "the output is binary - files are concatenated without parsing [%default]")

    parser.add_option( "--resubmit", dest="resubmit", type="int",
                       help = "if a job fails, automatically resubmit # times. Set to 0 in order to disable resubmission [%default]")

    parser.add_option( "--fail", dest="resubmit", action="store_false",
                       help = "if a job fails, do not resubmit [%default]")

    parser.add_option( "--bashrc", dest="bashrc", type="string",
                       help = "bashrc file to use [%default]")

    parser.add_option( "--method", dest="method", type="choice",
                       choices=("multiprocessing", "threads", "drmaa"),
                       help = "method to submit jobs [%default]")

    parser.add_option( "-e", "--env", dest="environment", type="string", action="append",
                       help = "environment variables to be passed to the jobs [%default]" )

    parser.add_option( "--output-pattern", dest="output_pattern", type="string",
                       help = "Pattern for secondary output filenames. Should contain a '%s' " 
                              "[default=%default]." )
    parser.set_defaults( 
        split_at_lines = None,
        split_at_column = None,
        split_at_regex = None,
        group_by_regex = None,
        split_at_tag = None,
        chunksize = None,
        cluster_cmd = 'qrsh -cwd -now n',
        bashrc = "~/.bashrc",
        input_header = False,
        output_header = False,
        debug = False,
        dry_run = False, 
        tmpdir = "./",
        subdirs = False,
        renumber = None,
        output_tag = "# job finished",
        map = [], 
        reduce = [],
        resume = None,
        renumber_column = [],
        resubmit = 5, 
        collect = None, 
        method = "drmaa",
        max_files = None,
        max_lines = None,
        binary = False,
        environment = [],
        output_pattern = "%s",
        )

    ## stop parsing options at the first argument
    parser.disable_interspersed_args()

    return parser
    
##--------------------------------------------------------------------
def main():

    parser = getOptionParser()

    (options, args) = E.Start( parser, 
                               add_cluster_options = True )

    if len(args) == 0:
        raise ValueError( "command line argument missing - see usage information" )

    options.renumber_column = [ x.split( ":" ) for x in options.renumber_column ]

        
    cmd = args[0]
    if len(args) > 1:
        cmd += " '" + "' '".join(args[1:]) + "'" 

    if options.dry_run:

        cmd = re.sub( "%DIR%", "", cmd )
        retcode = subprocess.call( cmd,
                                   shell = True,
                                   stdin = sys.stdin,
                                   stdout = sys.stdout,
                                   cwd = os.getcwd(),
                                   close_fds = True)                              
        E.Stop()
        sys.exit(0)

    failed_requests = []
    started_requests = []
    niterations = 0

    if not options.collect:
        tmpdir = os.path.abspath( tempfile.mkdtemp( dir = options.tmpdir ) )

        E.info( " working in directory %s" % tmpdir )

        if options.split_at_lines:
            chunk_iterator = chunk_iterator_lines
            args = (options.split_at_lines,)
        elif options.split_at_column:
            chunk_iterator = chunk_iterator_column
            args = (options.split_at_column - 1, options.max_files )
        elif options.split_at_regex:
            chunk_iterator = chunk_iterator_regex_split
            args = (re.compile( options.split_at_regex ),0,options.chunksize, options.max_lines )
        elif options.group_by_regex:
            chunk_iterator = chunk_iterator_regex_group
            args = (re.compile( options.group_by_regex ),0, options.chunksize)
        else:
            raise ValueError( "please specify a way to chunk input data" )

        data = [ (x, cmd, options, None, options.subdirs) for x in chunk_iterator( options.stdin, 
                                                                                   args, 
                                                                                   prefix = tmpdir,
                                                                                   use_header = options.input_header ) ]

        started_requests = [ (x[0], x[0]+ ".out") for x in data ]

        if len(data) == 0:
            E.warn("no data received")
            E.Stop()
            sys.exit(0)

        if options.method == "multiprocessing":
            pool = multiprocessing.Pool( options.cluster_num_jobs )
            results = pool.map( runCommand, data, chunksize = 1 )
        elif options.method == "drmaa":
            results = []
            runDRMAA( data, environment = options.environment )
        elif options.method == "threading":
            results = []
            def reportError( request, exc_info ):
                """report errors.
                """
                E.warn( "exception occured in request #%s: %s" % \
                            (request.requestID, exc_info[1]) )
                options.stdlog.flush()

            def saveResult( request, result ):
                """save result.
                """
                results.append( result )

            pool = threadpool.ThreadPool( options.cluster_num_jobs )
            for a in data:
                request = threadpool.WorkRequest( runCommand, 
                                                  args = (a,),
                                                  callback = saveResult,
                                                  exc_callback = reportError )

                pool.putRequest( request )
            pool.wait()
            
        niterations = 0
        for retcode, filename, cmd, logfile, iterations in results:
            niterations += iterations
            if not hasFinished( retcode, filename, options.output_tag, logfile ):
                failed_requests.append( (filename, cmd ))

    else:
        tmpdir = options.collect
        started_requests = [ (x[:-4], x) for x in glob.glob( tmpdir + "/*.out" ) ]

        E.info( "collecting %i files from %s" % (len(started_requests),
                                                 tmpdir ) )

    if failed_requests:
        for fn, cmd in failed_requests:
            E.error( "failed request: filename= %s, cmd= %s" % (fn,cmd) )
    else:
        E.info( "building result from %i parts" % len(started_requests) )

        if options.renumber:
            mapper = MapperLocal( pattern = options.renumber )
        else:
            mapper = MapperEmpty()

        ## deal with stdout
        name = None
        index = None

        for pattern, column in options.renumber_column:

            if re.search( pattern, "stdout"):
                try:
                    index = int(column) - 1
                except ValueError:
                    name = column
                    break

        if options.binary:
            ResultBuilderBinary()( started_requests, options.stdout, options )
        else:
            ResultBuilder( mapper = mapper, field_index = index, field_name = name )( started_requests, options.stdout, options )

        ## deal with logfiles : combine them into a single file
        x = re.search( "'--log=(\S+)'", cmd ) or re.search( "'--L\s+(\S+)'", cmd )
        if x:
            E.info( "logging output goes to %s" % x.groups()[0])
            logfile = IOtools.openFile(x.groups()[0], "a")
            ResultBuilderLog() ( [ (x[0], "%s.log" % x[0]) for x in started_requests ], logfile, options )
            logfile.close()

        ## deal with other files
        if options.subdirs:

            files = glob.glob( "%s/*.dir/*" % tmpdir )
            # remove directory
            filenames = set([ os.path.basename(x) for x in files ])
            xx = len(".out")

            for filename in filenames:

                _, filetype = os.path.splitext( filename )

                name = None
                index = None

                for pattern, column in options.renumber_column:
                    if re.search( pattern, filename):
                        try:
                            index = int(column) - 1
                        except ValueError:
                            name = column
                        break

                if options.binary:
                    builder = ResultBuilderBinary( mapper = mapper )
                elif filetype in ( ".fa", ".fasta" ):
                    builder = ResultBuilderFasta( mapper = mapper )
                elif filetype in ( ".mali", ):
                    builder = ResultBuilderFasta( mapper = MapperEmpty() )
                elif filetype in ( ".psl" ):
                    builder = ResultBuilderPSL( mapper = mapper )
                elif filetype in ( ".gtf", ".gff"):
                    builder = ResultBuilderGFF( mapper = mapper, field_index = index, field_name = name )
                elif filetype in ( ".png" ):
                    builder = ResultBuilderCopies( mapper = mapper )
                else:
                    builder = ResultBuilder( mapper = mapper, field_index = index, field_name = name )

                E.debug("chose the following builder for %s: %s: %s" % (filename, filetype, str(builder)))

                E.info( "collecting results for %s" % filename )

                input_filenames = []
                for fi, fn in started_requests:
                    fn = fn[:-xx] + ".dir/" + filename
                    if os.path.exists( fn ):
                        input_filenames.append( (fi, fn)  )

                E.info( "output of %i files goes to %s" % (len(filenames), filename ) )

                outfile = IOTools.openFile( options.output_pattern % filename, "w" )
                builder( input_filenames, outfile, options )
                outfile.close()

    if not options.debug and (not options.resume or not options.collect):
        if len(failed_requests) == 0:
            E.info("removing directory %s" % tmpdir )
            shutil.rmtree( tmpdir )
        else:
            E.info("directory %s not removed due to %i failed jobs" % (tmpdir, len(failed_requests) ) )
                
    E.info( "job control: nstarted=%i, nfinished=%i, nerrors=%i, nrepeats=%i" % \
                (len(started_requests),
                 len(started_requests) - len(failed_requests),
                 len(failed_requests),
                 niterations) )

    E.Stop()

if __name__ == '__main__':
    sys.exit(main())
