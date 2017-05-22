'''
split_links.py - 
======================================================

:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python split_links.py --help

Type::

   python split_links.py --help

for command line help.

Command line options
--------------------

'''
import os
import sys
import re
import CGAT.Experiment as E

open_files = {}


def WriteLine(a, b, line, prefix="%s-%s"):

    key1 = prefix % (a, b)
    key2 = prefix % (b, a)
    if key1 in open_files or os.path.exists(key1):
        key = key1
    else:
        key = key2

    if key not in open_files:
        open_files[key] = IOTools.openFile(key, "a")

    f = open_files[key]
    f.write(line)
    f.flush()


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: split_links.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-m", "--method", dest="method", type="string",
                      help="method for splitting.")

    parser.add_option("-r", "--regex", dest="regex", type="string",
                      help="regex to find prefix.")

    parser.add_option("-o", "--output-section", dest="output", type="string",
                      help="output filename.")

    parser.add_option("-t", "--targets", dest="targets", type="string",
                      help="output filename.")

    parser.set_defaults()

    (options, args) = E.Start(parser)

    if options.targets:
        options.targets = options.targets.split(",")

    nsame = 0
    ndiff = 0

    if options.method == "prefix":

        for line in sys.stdin:
            if line[0] == "#":
                continue

            data = line[:-1].split("\t")

            g1 = re.search(options.regex, data[0]).groups()[0]
            g2 = re.search(options.regex, data[1]).groups()[0]

            if g1 == g2:
                for t in options.targets:
                    if g1 == t:
                        continue
                    WriteLine(g1, t, line, options.output)
                    nsame += 1
            else:
                WriteLine(g1, g2, line, options.output)
                ndiff += 1

    print("nsame=%i, ndiff=%i" % (nsame, ndiff))

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
