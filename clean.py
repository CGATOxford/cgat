####
####
##
## Copyright (C) 2007 Andreas Heger All rights reserved
##
## Author: Andreas Heger <andreas.heger@dpag.ox.ac.uk>
##
## $Id: clean.py 2782 2009-09-10 11:40:29Z andreas $
##
##
####
####

USAGE="""python clean.py [OPTIONS] 

clean up output files from aborted runs.
"""

import os, sys, re, string, optparse, glob, subprocess, os.path

import Experiment

def getLastLine( filename, read_size = 1024 ):
  """return last line of a file.
  """
  f = open(filename, 'rU')    # U is to open it with Universal newline support
  offset = read_size
  f.seek(0, 2)
  file_size = f.tell()
  if file_size == 0: return ""
  while 1:
    if file_size < offset:
      offset = file_size
    f.seek(-1*offset, 2)
    read_str = f.read(offset)
    # Remove newline at the end
    if read_str[offset - 1] == '\n':
      read_str = read_str[:-1]
    lines = read_str.split('\n')
    if len(lines) >= 2:
        return lines[-1]
    if offset == file_size:   # reached the beginning
      return read_str
    offset += read_size
  f.close()

def checkPythonRuns( filename ):
    """returns true if a python run is complete."""
    last_line = getLastLine( filename )
    return re.match( "# job finished", last_line )
    
def isNewer( a, b ):
    """return true if file a is newer than file b."""

    ## get times of most recent access
    at = os.stat( a )[7]
    bt = os.stat( b )[7]

    return at > bt

if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id: clean.py 2782 2009-09-10 11:40:29Z andreas $", usage=USAGE )

    parser.add_option( "-g", "--glob", dest="glob_pattern", type="string" ,
                       help="glob pattern to use for collecting files [%default].")

    parser.add_option( "-n", "--dry-run", dest="dry_run", action="store_true",
                       help="only print out actions, do not execute them [%default]." )

    parser.add_option( "-f", "--file-pattern", dest="file_pattern", type="string",
                       help="only check files matching this pattern [%default]." )

    parser.set_defaults( glob_pattern = "data.dir",
                         file_pattern = ".out",
                         check_completeness = "python",
                         skip_dirs = [],
                         dry_run = False,
                         )

    (options, args) = Experiment.Start( parser, 
                                        add_pipe_options = True )

    if args:
        starts = args
    elif options.glob_pattern:
        starts = glob.glob( options.glob_pattern )
    else:
        starts = "."

    ndirs, nfiles, ndeleted = 0, 0, 0

    if options.check_completeness == "python":
        isComplete = checkPythonRuns

    rx = re.compile( options.file_pattern )

    for start in starts:
        for root, dirs, files in os.walk(start):
            
            ndirs += 1
            ## exclude directories
            for dir in options.skip_dirs:
                if dir in dirs:
                    dirs.remove(dir) 

            for filename in files:
                p = os.path.join(root,filename ) 
                if rx.search( filename ) and not isComplete( p ):
                    if options.dry_run:
                        options.stdlog.write("# removing file %s\n" % p )
                    else:
                        os.remove( p )
                    ndeleted += 1

    if options.loglevel >= 1:
        options.stdlog.write("# ndirs=%i, nfiles=%i, ndeleted=%i\n" %\
                                 (ndirs, nfiles, ndeleted) )

    Experiment.Stop()
