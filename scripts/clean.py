'''
clean.py - clean up output files from aborted runs
==================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script checks one or more output files have they
have completed successfully. It will remove output files
for those jobs that are incomplete.

The script checks for the "job finished" tag at the
end of the file.

Usage
-----

Example::

   python clean.py --help

Type::

   python clean.py --help

for command line help.

Command line options
--------------------
'''

import os
import sys
import re
import glob
import os.path
import CGAT.Experiment as E


def getLastLine(filename, read_size=1024):
    """return last line of a file.
    """
    f = open(
        filename, 'rU')    # U is to open it with Universal newline support
    offset = read_size
    f.seek(0, 2)
    file_size = f.tell()
    if file_size == 0:
        return ""
    while 1:
        if file_size < offset:
            offset = file_size
        f.seek(-1 * offset, 2)
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


def checkPythonRuns(filename):
    """returns true if a python run is complete."""
    last_line = getLastLine(filename)
    return re.match("# job finished", last_line)


def isNewer(a, b):
    """return true if file a is newer than file b."""

    # get times of most recent access
    at = os.stat(a)[7]
    bt = os.stat(b)[7]

    return at > bt


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: $Id: clean.py 2782 2009-09-10 11:40:29Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-g", "--glob", dest="glob_pattern", type="string",
                      help="glob pattern to use for collecting files [%default].")

    parser.add_option("-n", "--dry-run", dest="dry_run", action="store_true",
                      help="only print out actions, do not execute them [%default].")

    parser.add_option("-f", "--file-pattern", dest="file_pattern", type="string",
                      help="only check files matching this pattern [%default].")

    parser.set_defaults(glob_pattern="data.dir",
                        file_pattern=".out",
                        check_completeness="python",
                        skip_dirs=[],
                        dry_run=False,
                        )

    (options, args) = E.Start(parser,
                              add_pipe_options=True)

    if args:
        starts = args
    elif options.glob_pattern:
        starts = glob.glob(options.glob_pattern)
    else:
        starts = "."

    ndirs, nfiles, ndeleted = 0, 0, 0

    if options.check_completeness == "python":
        isComplete = checkPythonRuns

    rx = re.compile(options.file_pattern)

    for start in starts:
        for root, dirs, files in os.walk(start):

            ndirs += 1
            # exclude directories
            for dir in options.skip_dirs:
                if dir in dirs:
                    dirs.remove(dir)

            for filename in files:
                p = os.path.join(root, filename)
                if rx.search(filename) and not isComplete(p):
                    if options.dry_run:
                        options.stdlog.write("# removing file %s\n" % p)
                    else:
                        os.remove(p)
                    ndeleted += 1

    if options.loglevel >= 1:
        options.stdlog.write("# ndirs=%i, nfiles=%i, ndeleted=%i\n" %
                             (ndirs, nfiles, ndeleted))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
