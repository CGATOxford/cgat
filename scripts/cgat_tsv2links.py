'''
cgat_tsv2links.py - create softlinks to a series files
========================================================

:Author: Mike Morgan
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Create links from a tab-separated table. The table should have two
columns, for example::

   source     dest
   abc.fq.gz  sample1-R1.fq.gz
   def.fq.gz  sample1-R2.fq.gz

If the options ``--source`` is given, filenams will be matched
by walking through ``--source``. Otherwise, filenames in the source column 
need to contain the paths relative to the current working directory.

To create such a table, use the unix ``find`` command, for example::

   find /ifs/projects/proj013/backup/ -name "*.sanfastq.gz" > input_file.tsv
   
and then manually add table headers and a second column with the sample name.

Further ways to develop the script:
 
   * paired files - files might be grouped (read pairs), make sure
            they are all there. Use pattern matching to identify a
            group and create all appropriate links.

    *  requires exception handling if no options of files are provided

Usage
-----

Example::

   python cgat_tsv2links.py --source=../backup/dataset1 < input_file.tsv

Type::

   python cgat_tsv2links.py --help

for command line help.

Command line options
--------------------


'''

import os
import sys
import re
import optparse
import CGAT.Experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-s", "--source", dest="source_directory",
                      type="string", default=False,
                      help="The directory in which data"
                      "files are held [%default]")

    parser.add_option("-d", "--dest", dest="dest_directory",
                      type="string", default=False,
                      help="The directory in which links"
                      "are created [%default]")

    parser.set_defaults(source_directory=None,
                        dest_directory=".")
    
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    # read a map of input files to links with sanity checks
    map_filename2link = {}
    links = set()
    for line in options.stdin:
        if line.startswith("#"):
            continue

        # ignore header
        if line.startswith("source"):
            continue

        filename, link = line[:-1].split()[:2]
        if filename in map_filename2link:
            raise ValueError("duplicate filename '%s' " % filename)
        if link in links:
            raise ValueError("duplicate link '%s' " % link)
        map_filename2link[filename] = link
        links.add(link)

    counter = E.Counter()
    counter.input = len(map_filename2link)

    def _createLink(src, dest, counter):
        
        src = os.path.abspath(src)
        dest = os.path.abspath(os.path.join(options.dest_directory, dest))
        if os.path.exists(dest):
            E.warn("existing symlink %s" % dest)
            counter.link_exists += 1
        elif not os.path.exists(src):
            counter.file_not_found += 1
            E.warn("did not find %s" % src)
        else:
            try:
                os.symlink(src, dest)
                counter.success += 1
            except OSError:
                pass

    if not options.source_directory:
        # no source directory given, filenames must have complete path
        for filename, link in map_filename2link.items():
            _createLink(filename, link, counter)
    else:
        # walk through directory hierchy and create links
        # for files matching filenames in map_filename2link
        found = set()
        for dirName, subdirList, fileList in os.walk(options.source_directory):
            for f in fileList:
                if f in map_filename2link:
                    if f in found:
                        E.warn("found multiple files with "
                               "the same name %s" % f)
                    else:
                        _createLink(os.path.join(dirName, f), 
                                    map_filename2link[f], counter)
                        found.add(f)
                else:
                    E.info("Filename %s not in map" % f)

        notfound = set(map_filename2link.keys()).difference(found)
        counter.notfound = len(notfound)
        if notfound:
            E.warn("did not find %i files: %s" % (len(notfound),
                                                  str(notfound)))

    E.info(counter)
    # write footer and output benchmark information
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
