##########################################################################
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
##########################################################################
'''
optic/update_orthology.py - 
======================================================

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

   python optic/update_orthology.py --help

Type::

   python optic/update_orthology.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import re
import glob
import CGAT.Experiment as E

""" program $Id: optic/update_orthology.py 2781 2009-09-10 11:33:14Z andreas $

update orthology assignments.

"""


def createLink(old_fn, new_fn):
    """create symbolic link between old and new."""
    if not os.path.exists(new_fn):
        try:
            os.symlink(old_fn, new_fn)
        except OSError:
            print "# ERROR: can't create symling to %s" % old_fn
            return False

    return True


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/update_orthology.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-o", "--old-dir", dest="old_directory", type="string",
                      help="old directory.")
    parser.add_option("-n", "--new-dir", dest="new_directory", type="string",
                      help="new directory.")
    parser.add_option("-s", "--schemas=", dest="schemas", type="string",
                      help="schemas.")

    parser.set_defaults(
        old_directory=None,
        new_directory=".",
        schemas=""
    )

    (options, args) = E.Start(parser)

    options.schemas = options.schemas.split(",")

    options.old_directory = os.path.abspath(options.old_directory)
    options.new_directory = os.path.abspath(options.new_directory)

    # locate old data
    rx = re.compile("^([^0-9]+)(\d+)")
    if options.old_directory:
        data = glob.glob("%s/pair_*" % options.old_directory)
        species = {}
        old_data = {}
        for d in data:
            file = os.path.basename(d)
            species1, species2 = file[5:].split("-")
            species[species1] = True
            species[species2] = True
            # extract species name and version
            s1, v1 = rx.search(species1).groups()
            s2, v2 = rx.search(species2).groups()
            key = "%s-%s" % (s1, s2)
            if key not in old_data:
                old_data[key] = []
            old_data[key].append((species1, species2, d))
    else:
        old_data = {}

    # set up links
    # todo

    nerrors, nlinked, nupdate, nnew = 0, 0, 0, 0
    ninput = 0

    for x in range(len(options.schemas) - 1):
        for y in range(x + 1, len(options.schemas)):
            s1 = options.schemas[x]
            s2 = options.schemas[y]

            ninput += 1

            is_found = False
            is_identical = False

            if not os.path.exists("%s/pair_%s-%s" % (options.new_directory, s1, s2)):
                if not os.path.exists("%s/pair_%s-%s" % (options.new_directory, s2, s1)):
                    print "directory pair_%s-%s or pair_%s-%s does not exist in %s!" % (s1, s2, s2, s1, options.new_directory)
                else:
                    s1, s2 = s2, s1

            g1 = rx.search(s1).groups()[0]
            g2 = rx.search(s2).groups()[0]
            old_dir = None

            if ("%s-%s" % (g1, g2)) in old_data:
                is_found = True
                src = old_data["%s-%s" % (g1, g2)]
                for ss1, ss2, old_dir in src:
                    if ss1 == s1 and ss2 == s2:
                        is_identical = True

            elif ("%s-%s" % (g2, g1)) in old_data:
                is_found = True,
                src = old_data["%s-%s" % (g2, g1)]
                for ss1, ss2, old_dir in src:
                    if ss1 == s2 and ss2 == s1:
                        is_identical = True

            new_dir = "pair_%s-%s" % (s1, s2)

            if is_found:
                if is_identical:
                    if options.loglevel >= 1:
                        print "# substituting directory %s with %s" % (new_dir, old_dir)
                    if os.path.islink(new_dir):
                        print "# dir %s is already linked. Skipped, but you make sure if link is correct." % new_dir
                    else:
                        if os.path.exists(new_dir):
                            if options.loglevel >= 2:
                                print "# removing everything under %s" % new_dir
                            for root, dirs, files in os.walk(new_dir, topdown=False):
                                for name in files:
                                    os.remove(os.path.join(root, name))
                                for name in dirs:
                                    os.rmdir(os.path.join(root, name))
                            os.rmdir(new_dir)
                        os.symlink(old_dir, new_dir)
                    nlinked += 1

                else:
                    is_ok = True

                    if options.loglevel >= 1:
                        print "# updating %s from %s" % (new_dir, old_dir)

                    is_ok = is_ok and createLink("%s/step1.dir/blast.links.gz" % old_dir,
                                                 "%s/previous.links" % new_dir)

                    is_ok = is_ok and createLink("%s/step1.queries.fasta" % old_dir,
                                                 "%s/previous.fasta" % new_dir)

                    if not is_ok:
                        nerrors += 1
                        continue

                    nupdate += 1
            else:
                if options.loglevel >= 1:
                    print "# de-novo calculation for %s" % new_dir
                nnew += 1

    if options.loglevel >= 1:
        options.stdlog.write("# ninput=%i, nlinked=%i, nupdate=%i, nnew=%i, nerrors=%i\n" % (
            ninput, nlinked, nupdate, nnew, nerrors))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
