'''cgat_zap.py - replace files with empty placeholders
===================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

The purpose of this script is to help in preserving disk space by
replacing certain files with empty files. The new files will have the
same file properties (modification time, etc) as the original file in
order to preserve dependencies based on time-stamps in computational
work flows.

If a file is a link to another file, the link will be removed and an
empty file created in its place. The empty file will receive the file
properties of the file that was linked to.

Options
-------

Usage
-----

Example::

   python cgat_zap.py *.bam

Type::

   python cgat_zap.py --help

for command line help.

Command line options
--------------------

'''

import sys
import os
import CGAT.Experiment as E


def main(argv=None):
    """script main.

parses command line options in sys.argv, unless *argv* is given.
"""

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-n", "--dry-run", dest="dry_run", action="store_true",
                      help="do dry run, do not kill [default=%default].")

    parser.set_defaults(
        dry_run=False,
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    outfile = options.stdout

    fields = ('st_atime', 'st_blksize', 'st_blocks',
              'st_ctime', 'st_dev', 'st_gid', 'st_ino',
              'st_mode', 'st_mtime', 'st_nlink',
              'st_rdev', 'st_size', 'st_uid')

    outfile.write("filename\tlinkdest\t%s\n" % "\t".join(fields))

    # remove any duplicates and sort
    args = sorted(set(args))

    for fn in args:

        # stat follows times to links
        original = os.stat(fn)

        if os.path.islink(fn):
            linkdest = os.readlink(fn)
            E.info('breaking link from %s to %s' % (fn, linkdest))
            if not options.dry_run:
                os.unlink(fn)
                f = open(fn, "w")
                f.close()
        else:
            E.info('truncating file %s' % fn)
            linkdest = ""
            if not options.dry_run:
                f = open(fn, "w")
                f.truncate()
                f.close()

        outfile.write("%s\t%s\t%s\n" % (
            fn,
            linkdest,
            "\t".join([str(getattr(original, x)) for x in fields])))

        if not options.dry_run:
            # Set original times
            os.utime(fn, (original.st_atime, original.st_mtime))
            os.chmod(fn, original.st_mode)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
