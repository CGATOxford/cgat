'''
combine_files.py - group and merge multiple files
=================================================

:Author:
:Tags: Python

Purpose
-------

This script merges multiple files. The files are grouped
by some criteria and a concatenated file is output for
each group.

Usage
-----

For example::

   python combine_files.py \
         --map=chipseq.map \
         --suffix=_1.fastq.gz=.fastq.1.gz \
         --suffix=_2.fastq.gz=.fastq.2.gz \
         ../backup/ChipSeq-2014-03/*.fastq.gz

The above command will merge 8 files in the directory
:file:`../backup/ChipSeq-2014-03`::

      ../backup/Irf5ChipSeq-2014-03/WTCHG_112417_217_1.fastq.gz
      ../backup/Irf5ChipSeq-2014-03/WTCHG_112417_217_2.fastq.gz
      ../backup/Irf5ChipSeq-2014-03/WTCHG_112417_218_1.fastq.gz
      ../backup/Irf5ChipSeq-2014-03/WTCHG_112417_218_2.fastq.gz
      ../backup/Irf5ChipSeq-2014-03/WTCHG_112418_217_1.fastq.gz
      ../backup/Irf5ChipSeq-2014-03/WTCHG_112418_217_2.fastq.gz
      ../backup/Irf5ChipSeq-2014-03/WTCHG_112418_218_1.fastq.gz
      ../backup/Irf5ChipSeq-2014-03/WTCHG_112418_218_2.fastq.gz

The mapping is given by the file :file:`chipseq.map`::

    regex   destination
    WTCHG_112417_217        Replicate1
    WTCHG_112417_218        Replicate2
    WTCHG_112418_217        Replicate1
    WTCHG_112418_218        Replicate2

The script will output 4 files::

    Replicate1.fastq.1.gz
    Replicate1.fastq.2.gz
    Replicate2.fastq.1.gz
    Replicate2.fastq.2.gz

Type::

   python combine_files.py --help

for command line help.

Command line options
--------------------

'''

import sys
import re
import collections

import CGAT.Experiment as E
import CGAT.IOTools as IOTools


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-s", "--suffix", dest="suffixes",
                      action="append",
                      type="string",
                      help="mapping of suffix [%default]")

    parser.add_option("-m", "--map", dest="filename_map",
                      type="string",
                      help="filename containing a mapping of filenames "
                      "to destinations [%default]")

    parser.add_option("-u", "--no-sort", dest="do_sort",
                      action="store_false",
                      help="do not sort filenames before grouping "
                      "[%default]")

    parser.add_option("-n", "--dry-run", dest="dry_run",
                      action="store_true",
                      help="dry run - do not merge "
                      "[%default]")

    parser.set_defaults(
        filename_map=None,
        do_sort=True,
        dry_run=False,
        suffixes=[])

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.filename_map:
        map_regex2dest = IOTools.readMap(
            IOTools.openFile(options.filename_map))
        map_regex2dest = dict([(re.compile(x), y) for x, y in
                               list(map_regex2dest.items())])

    map_suffix2dest = {}
    for suffix in options.suffixes:
        src, dest = suffix.split("=")
        map_suffix2dest[src.strip()] = dest.strip()

    filenames = args

    if options.do_sort:
        filenames.sort()

    dest2src = collections.defaultdict(list)
    for filename in filenames:
        dests = []
        for regex, dest in list(map_regex2dest.items()):
            if regex.search(filename):
                dests.append(dest)
        if len(dests) == 0:
            raise ValueError("no destination found for %s" % filename)
        elif len(dests) > 1:
            raise ValueError(
                "multiple destinations found for %s: %s " %
                (filename, dests))

        dest = dests[0]
        # implement suffix mapping, note that
        # suffixes can extend beyond an extension
        for suffix, new_suffix in list(map_suffix2dest.items()):
            if filename.endswith(suffix):
                if suffix in map_suffix2dest:
                    dest = dest + map_suffix2dest[suffix]
                    break

        dest2src[dest].append(filename)

    for dest, srcs in sorted(dest2src.items()):
        E.info("merging: %s <- %s" % (dest, srcs))
        if options.dry_run:
            continue
        E.run('cat %s > %s' % (" ".join(srcs),
                               dest))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
