'''
malis2profiles.py - build profiles from malis
=============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

convert a set of plain alignments to profiles.

Usage
-----

Example::

   python malis2profiles.py --help

Type::

   python malis2profiles.py --help

for command line help.

Command line options
--------------------

'''
import sys
import CGAT.Experiment as E
import CGAT.Mali as Mali


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.set_defaults(
    )

    (options, args) = E.Start(parser)

    mali = Mali.SequenceCollection()
    last_id = None
    ninput, noutput, nskipped = 0, 0, 0

    for line in sys.stdin:
        if line[0] == "#":
            continue

        start, ali, end, id = line[:-1].split("\t")
        ninput += 1
        if id != last_id:
            if last_id:
                mali.setName(last_id)
                mali.writeToFile(sys.stdout, format="profile")
                noutput += 1
            mali = Mali.SequenceCollection()
            last_id = id

        mali.addSequence(id, start, end, ali)

    if last_id:
        mali.setName(last_id)
        mali.writeToFile(sys.stdout, format="profile")
        noutput += 1

    if options.loglevel >= 1:
        options.stdlog.write(
            "# ninput=%i, noutput=%i, nskipped=%i.\n" % (ninput, noutput, nskipped))

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
