'''
graph_blast2besthits.py - 
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

   python graph_blast2besthits.py --help

Type::

   python graph_blast2besthits.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re
import getopt
import CGAT.Experiment as E
import CGAT.BlastAlignments as BlastAlignments


USAGE = """python %s [OPTIONS] < graph.in > graph.out

Version: $Id: graph_blast2besthits.py 2782 2009-09-10 11:40:29Z andreas $

Parse blast graph and only output best hits for each query per
genome. Note that the blast graph has been sorted by query_token and
sbjct_token.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-p, --pattern-genome=           pattern to extract genome
-m, --method=                   method [blast|score|pid]
-s, --score-factor=             thresholding on score (multiplicative)
-i, --pide-factor=              thresholding on pide (additive)
""" % sys.argv[0]

param_long_options = ["help", "verbose=", "pattern-genome=", "method=",
                      "score-factor=", "pide-factor=",
                      "version"]

param_short_options = "hv:p:f:i"

param_loglevel = 1
param_method = "evalue"
param_filename_self_scores = None
param_pattern_genome = "^([^|]+)|"
param_method = "evalue"

param_score_threshold_factor = 1.0
param_pide_threshold_factor = 0.0

# ------------------------------------------------------------------------


def PrintMatches(matches, method=None, filter=None,
                 score_threshold_factor=1.0,
                 pide_threshold_factor=0.0):
    """print best match per organism (sorted by E-Value).

    If there are several matches with the same E-Value, all are printed.
    """

    orgs = matches.keys()
    orgs.sort()

    if method == "evalue":
        f1 = lambda x, y: cmp(x.mEvalue, y.mEvalue)
        f2 = lambda x: x.mEvalue
    elif method == "score":
        f1 = lambda x, y: -cmp(x.score, y.score)
        f2 = lambda x: x.score
    elif method == "pid":
        f1 = lambda x, y: -cmp(x.mPercentIdentity, y.mPercentIdentity)
        f2 = lambda x: x.mPercentIdentity
    else:
        raise "unknown method"

    noutput = 0

    for org in orgs:

        if filter and filter == org:
            continue

        m = matches[org]
        score_threshold = int(m[0].score * score_threshold_factor)
        pide_threshold = int(m[0].mPercentIdentity + pide_threshold_factor)

        if param_loglevel >= 2:
            print "# %i matches for %s, min_score=%i" % (len(m), org, score_threshold)

        if param_loglevel >= 3:
            for mm in m:
                print str(mm)

        m.sort(f1)

        s = f2(m[0])
        x = 0

        if score_threshold_factor != 1.0 and pide_threshold_factor >= 0.0:
            # take best match and all within score threshold factor
            for x in range(len(m)):
                if m[x].score > score_threshold or \
                        m[x].mPercentIdentity >= pide_threshold:
                    print str(m[x])
                    noutput += 1
        else:
            # write all with same primary score
            while x < len(m) and s == f2(m[x]) and m[x].score >= score_threshold:
                print str(m[x])
                x += 1
                noutput += 1

    return noutput

# ------------------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    try:
        optlist, args = getopt.getopt(
            sys.argv[1:], param_short_options, param_long_options)
    except getopt.error, msg:
        print USAGE, msg
        sys.exit(2)

    for o, a in optlist:
        if o in ("-v", "--verbose"):
            param_loglevel = int(a)
        elif o in ("--version", ):
            print "version="
            sys.exit(0)
        elif o in ("-h", "--help"):
            print USAGE
            sys.exit(0)
        elif o in ("-p", "--pattern-genome"):
            param_pattern_genome = a
        elif o in ("-m", "--method"):
            param_method = a
        elif o in ("-s", "--score-factor"):
            param_score_threshold_factor = float(a)
        elif o in ("-i", "--pide-factor"):
            param_pide_threshold_factor = float(a)

    print E.GetHeader()
    print E.GetParams()

    ninput, noutput, nskipped, nfailed = 0, 0, 0, 0

    last_query_token = None
    rx = re.compile(param_pattern_genome)

    matches = {}

    for line in sys.stdin:

        if line[0] == "#":
            continue

        link = BlastAlignments.Link()
        try:
            link.Read(line)
        except ValueError:
            sys.stderr.write("parsing error in line %s\n" % line[:-1])
            continue

        ninput += 1
        if link.mQueryToken != last_query_token:
            if last_query_token:
                noutput += PrintMatches(matches,
                                        method=param_method,
                                        filter=rx.search(
                                            last_query_token).groups()[0],
                                        score_threshold_factor=param_score_threshold_factor,
                                        pide_threshold_factor=param_pide_threshold_factor)
            matches = {}
            last_query_token = link.mQueryToken

        org = rx.search(link.mSbjctToken).groups()[0]

        if org not in matches:
            matches[org] = []
        matches[org].append(link)

    noutput += PrintMatches(matches,
                            method=param_method,
                            filter=rx.search(last_query_token).groups()[0],
                            score_threshold_factor=param_score_threshold_factor,
                            pide_threshold_factor=param_pide_threshold_factor)

    print "# ninput=%i, noutput=%i" % (ninput, noutput)

    print E.GetFooter()

    if noutput == 0:
        if ninput == 0:
            raise "no output, because no input."
        else:
            raise "no output."


if __name__ == "__main__":
    sys.exit(main(sys.argv))
