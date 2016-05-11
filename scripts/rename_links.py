'''
rename_links.py - rename all links in a directory
=================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script reads links in a directory that match a pattern
and changes the destination into a new pattern keeping 
name of the link unchanged.

Usage
-----

Example::

   python rename_links.py --help

Type::

   python rename_links.py --help

for command line help.

Command line options
--------------------

'''
import os
import sys
import re
import CGAT.Experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: $Id: rename_links.py 2782 2009-09-10 11:40:29Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-r", "--recursive", dest="recursive", action="store_true",
                      help="proceed recursively through directories.")

    parser.add_option("-d", "--dry-run", dest="dry_run", action="store_true",
                      help="dry run, do not change anything.")

    parser.add_option("-o", "--old-pattern", dest="old_pattern", type="string",
                      help="pattern to be replaced.")

    parser.add_option("-n", "--new-pattern", dest="new_pattern", type="string",
                      help="pattern to be inserted.")

    parser.set_defaults(
        recursive=False,
        dry_run=False,
        old_pattern=None,
        new_pattern=None,
    )

    (options, args) = E.Start(parser)

    if options.old_pattern is None or options.new_pattern is None:
        raise "please specify both an old and a new pattern."

    rx_old = re.compile(options.old_pattern)

    def changeLink(fn, options):
        old_target = os.path.realpath(fn)
        if not re.search(options.old_pattern, old_target):
            return 0
        new_target = re.sub(
            options.old_pattern, options.new_pattern, old_target)

        cmds = []
        cmds.append('rm -f %s' % (fn))

        cmds.append('ln -fs %s %s' % (new_target, fn))

        if options.loglevel >= 1 or options.dry_run:
            for cmd in cmds:
                options.stdlog.write("%s\n" % cmd)

        if not options.dry_run:
            for cmd in cmds:
                os.system(cmd)

        return 1

    ninput, nvisited, nchanged = 0, 0, 0
    for directory in args:
        ninput += 1
        E.info("processing directory %s" % directory)

        for root, directories, files in os.walk(directory):
            for f in files + directories:
                fn = os.path.join(root, f)
                if not os.path.islink(fn):
                    continue
                nchanged += changeLink(os.path.abspath(fn), options)
                nvisited += 1

    if options.stdlog >= 1:
        options.stdlog.write(
            "# ndirs=%i, nvisited=%i, nchanged=%i\n" % (ninput, nvisited, nchanged))

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
