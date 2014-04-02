'''
graph_map_links.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Map blast links.

If ``--map-identity`` is set, only identifiers will be updated.
Otherwise, the alignments will be combined.

Usage
-----

Example::

   python graph_map_links.py --help

Type::

   python graph_map_links.py --help

for command line help.

Command line options
--------------------

'''
import os
import sys
import string
import re
import getopt
import tempfile
import time
import popen2
import optparse
import hashlib


import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.BlastAlignments as BlastAlignments


def readIdentityMap(infile):
    """read identity map.

    multiple entries can either be separated over several lines
    or concatenated by semicolon in the same line.
    """

    m = {}
    for line in infile:

        if line[0] == "#":
            continue
        olds, news = map(lambda x: x.split(";"), line[:-1].split("\t")[:2])

        for old in olds:
            if old not in m:
                m[old] = []
            for new in news:
                m[old].append(new)
    return m

# ------------------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: graph_map_links.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-q", "--map-query", dest="filename_map_query", type="string",
                      help="filename with queries to map.")

    parser.add_option("-s", "--map-sbjct", dest="filename_map_sbjct", type="string",
                      help="filename with sbjcts to map.")

    parser.add_option("-m", "--multiple", dest="multiple", action="store_true",
                      help="map multiple options [%default].")

    parser.add_option("-k", "--keep-unmapped", dest="keep_unmapped", action="store_true",
                      help="keep unmapped entries [%default].")

    parser.add_option("-i", "--map-identity", dest="map_identity", action="store_true",
                      help="map by identifier [%default].")

    parser.add_option("-n", "--non-redundant", dest="non_redundant", action="store_true",
                      help="write only unique links (requires a lot of memory for large graphs) [%default]")

    parser.set_defaults(
        filename_map_query=None,
        filename_map_sbjct=None,
        multiple=False,
        keep_unmapped=False,
        map_identity=False,
        report_step=1000000,
        non_redundant=False)

    (options, args) = E.Start(parser)

    if options.filename_map_query:
        infile = IOTools.openFile(options.filename_map_query, "r")
        if options.map_identity:
            map_query = readIdentityMap(infile)
        else:
            map_query = BlastAlignments.ReadMap(infile, options.multiple)
        infile.close()
        E.info('read maps for %i queries' % len(map_query))
    else:
        map_query = None

    if options.filename_map_sbjct:
        if options.filename_map_sbjct == options.filename_map_query:
            map_sbjct = map_query
        else:
            infile = IOTools.openFile(options.filename_map_sbjct, "r")
            if options.map_identity:
                map_sbjct = readIdentityMap(infile)
            else:
                map_sbjct = BlastAlignments.ReadMap(infile, options.multiple)
            infile.close()
        E.info('read maps for %i sbjcts' % len(map_sbjct))
    else:
        map_sbjct = None

    nfailed = 0
    ninput = 0
    nskipped = 0
    noutput = 0

    # number of identical/mapped links
    nsame, nmapped = 0, 0

    printed = {}

    alignment = BlastAlignments.Map()

    for line in options.stdin:

        if line[0] == "#":
            continue

        data = line[:-1].split("\t")

        alignment.Read(line)
        skip = False
        ninput += 1

        E.debug(str(map))

        if options.loglevel >= 2 and ninput % options.report_step == 0:
            options.stderr.write(
                "# progress: ninput=%i, noutput=%i, nhash=%i\n" % (ninput, noutput, len(printed)))

        if options.multiple:
            skip = False
            if map_query is not None:
                if alignment.mQueryToken in map_query:
                    mq = map_query[alignment.mQueryToken]
                else:
                    skip = True
            else:
                mq = [None]

            if map_sbjct is not None:
                if alignment.mSbjctToken in map_sbjct:
                    ms = map_sbjct[alignment.mSbjctToken]
                else:
                    skip = True
            else:
                ms = [None]

            if skip:
                nskipped += 1
                continue

            if options.map_identity:

                # only if non_redundant is set, do global comparison
                if not options.non_redundant:
                    printed = {}

                new_map = alignment.GetClone()
                do_redundant = len(mq) > 1 or len(ms) > 1
                for q in mq:
                    for s in ms:

                        new_map.mQueryToken = q
                        new_map.mSbjctToken = s

                        # check for non-redundant links for 1:many or many:many
                        # mappings
                        if do_redundant:
                            key = "%s-%i-%i-%s-%i-%i" % (new_map.mQueryToken, new_map.mQueryFrom, new_map.mQueryTo,
                                                         new_map.mSbjctToken, new_map.mSbjctFrom, new_map.mSbjctTo)

                            # hash key to save space
                            hkey = hashlib.md5(key).digest()

                            if hkey in printed:
                                continue

                            printed[hkey] = 1

                        options.stdout.write(
                            '\t'.join([str(new_map)] + data[9:]) + '\n')
                        noutput += 1
                        if new_map.mQueryToken == alignment.mQueryToken and \
                                new_map.mSbjctToken == alignment.mSbjctToken:
                            nsame += 1
                        else:
                            nmapped += 1

            else:
                for q in mq:
                    for s in ms:
                        new_map = alignment.GetClone()

                        E.debug(str(q))
                        E.debug(str(s))

                        is_ok = new_map.MapAlignment(q, s)

                        if not is_ok:
                            nfailed += 1
                        else:
                            options.stdout.write(
                                '\t'.join([str(new_map)] + data[9:]) + '\n')
                            noutput += 1

        # options.multiple is False
        else:

            if map_query is not None:
                if alignment.mQueryToken in map_query:
                    mq = map_query[alignment.mQueryToken]
                else:
                    mq = None
                    skip = True
            else:
                mq = None

            if map_sbjct is not None:
                if alignment.mSbjctToken in map_sbjct:
                    ms = map_sbjct[alignment.mSbjctToken]
                else:
                    ms = None
                    skip = True
            else:
                ms = None

            if skip and not options.keep_unmapped:
                nskipped += 1
                continue

            E.debug(str(mq))
            E.debug(str(ms))

            if mq or ms:
                is_ok = alignment.MapAlignment(mq, ms)
            else:
                is_ok = True

            if not is_ok:
                nfailed += 1
            else:
                options.stdout.write(
                    '\t'.join([str(alignment)] + data[9:]) + '\n')
                noutput += 1

    E.info('ninput=%i, noutput=%i, nskipped=%i, nfailed=%i, nsame=%i, nmapped=%i' %
           (ninput, noutput, nskipped, nfailed, nsame, nmapped))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
