'''
add_cgat_script_template.py - add template to python scipts
===========================================================

:Tags: Python

Purpose
-------

adds a preamble to python scripts. The original scripts
are saved as ``.bak`` files.

Usage
-----

Example::

   python add_cgat_script_template.py script1.py script2.py

Type::

   python add_cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import shutil

import CGAT.Experiment as E

SCRIPT = """################################################################################
'''
%(filename)s - 
======================================================

:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python %(filename)s --help

Type::

   python %(filename)s --help

for command line help.

Documentation
-------------

Code
----

'''
"""

MODULE = """################################################################################
'''
%(filename)s - 
======================================================

:Tags: Python

Code
----

'''
"""


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--template", dest="template", type="choice",
                      choices=("script", "module"),
                      help="which template to choose [default=%default].")

    parser.set_defaults(
        template="script",
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.template == "script":
        template = SCRIPT
    elif options.template == "module":
        template = MODULE

    for filename in args:
        E.info("processing %s" % filename)
        lines = IOTools.openFile(filename).readlines()
        x = 0
        while x < len(lines) and lines[x][0] in ("#", "", "\n"):
            x += 1
        l = lines[x]
        if l.startswith( "'''" ) or l.startswith('"""'):
            continue
        E.info("moving %s to %s.bak" % (filename, filename))
        shutil.move(filename, "%s.bak" % filename)
        outfile = IOTools.openFile(filename, "w")
        outfile.write(template % locals())
        outfile.write("".join(lines[x:]))
        outfile.close()

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
