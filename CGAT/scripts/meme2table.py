'''
meme2table - summarize infomation about a meme motif file
=========================================================

:Author:
:Tags: Python

Purpose
-------

.. Overall purpose and function of the script>

Usage
-----

.. Example use case

Example::

   python cgat_script_template.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.MEME as MEME

import itertools


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    motifs = MEME.MemeMotifFile(options.stdin)
    headers = list(set(itertools.chain(*[motif.properties.keys()
                                         for motif in motifs])))

    outlines = []
    for motif in motifs:
        outline = [motif.primary_id, motif.secondary_id, motif.consensus()]
        outline.extend([motif.properties.get(col, "") for col in headers])
        outlines.append(outline)

    headers = "\t".join(["primary_id", "secondary_id", "consensus"] + headers)

    output = "\n".join(["\t".join(map(str, line))
                        for line in outlines])
    output = headers + "\n" + output
    options.stdout.write(output)
        
    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
