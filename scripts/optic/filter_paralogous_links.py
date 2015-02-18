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
optic/filter_paralogous_links.py - 
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

   python optic/filter_paralogous_links.py --help

Type::

   python optic/filter_paralogous_links.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import getopt
import CGAT.Experiment as E
import CGAT.Exons as Exons
import alignlib_lite
import CGAT.BlastAlignments as BlastAlignments


USAGE = """python %s [OPTIONS] < orthologs > genes

Version: $Id: optic/filter_paralogous_links.py 1799 2008-03-28 11:44:19Z andreas $

Remove from a list of alignments those between paralogous sequences.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
--cds-gtf-file                           files with coding sequences (exon boundaries)
--write-exons=                  write exons (unaligned/aligned)
--write-introns=                write introns (unaligned/aligned)
--extend-introns=               extend intronic sequences with x residues from adjacent exons
--only-best                     take only best ortholog transcripts
--mask                          mask softmasked sequences
--min-intron-length=            minimum length of intron to align
--max-intron-length=            maximum length of intron to align
--missing-max-missing=          maximum number of missing introns
--missing-min-present=          minimum of introns present
--boundaries-max-slippage=      maximum of intron boundaries that they can deviate from reference
--boundaries-max-missed=        maximum number of intron boundaries that can be missed
--boundaries-allow-missed=      maximum number of intron boundaries that can be missed
--compress                      write alignment in compressed form
--disable-check-exon-number     do not check exon number
--report-step                   dump progress at each umpth step.
""" % sys.argv[0]

param_long_options = ["verbose=", "help",
                      "cds=",
                      "write-exons=", "write-introns=",
                      "extend-introns=", "only-best", "mask",
                      "min-intron-length=", "max-intron-length=",
                      "boundaries-max-missed=",
                      "boundaries-max-slippage=",
                      "boundaries-allow-missed=",
                      "missing-max-missing=", "missing-min-present=",
                      "disable-check-exon-number",
                      "report-step=",
                      "compress", "version"]


param_short_options = "v:hg:"

param_loglevel = 2

param_genome_file = "genome_%s.fasta"

param_filename_peptides1 = None
param_filename_peptides2 = None
param_filename_map1 = None
param_filename_map2 = None
param_filename_transcripts1 = None
param_filename_transcripts2 = None

# maximum number of gaps for ortholog pairs
param_max_gaps = 10

param_report_step = 100000

# whether or not to check number of introns
param_check_exon_number = 1

# whether or not to check number of introns
param_check_exon_boundaries = 1

# allowable margin for exon boundaries to vary.
param_boundaries_max_slippage = 3

# allowable number of exon boundaries that are wrong
param_boundaries_max_missed = 0

# above this number of exon boundiares, some boundaries may be missed:
param_boundaries_allow_missed = 0

##
param_boundaries_require_found = 1

param_compress = 0

# allowable

# flags of output
param_write_introns = ("")
param_write_exons = ("")

param_extend_introns = 0

param_min_coverage = 50.0

param_only_best = 0

param_mask = None

param_min_intron_length = None
param_max_intron_length = None

param_missing_min_present = 0
param_missing_max_missing = 0

# keep smatrix global, so that data does not get deleted when object goes
# out of scope!
global_substitution_matrix = None

param_min_score_sw = 50


# ------------------------------------------------------------
def UniquifyList(o):
    o.sort()
    n = []
    l = None
    for x in o:
        if x != l:
            n.append(x)
        l = x
    return n

# ------------------------------------------------------------


def IsParalogLink(link, cds1, cds2):
    """sort out ortholog relationships between
    transcripts of orthologous genes.

    """

    map_a2b = alignlib_lite.makeAlignmentVector()
    alignlib_lite.AlignmentFormatEmissions(
        link.mQueryFrom, link.mQueryAli,
        link.mSbjctFrom, link.mSbjctAli).copy(map_a2b)

    if link.mQueryLength < (map_a2b.getRowTo() - map_a2b.getRowFrom() + 1) or \
       link.mSbjctLength < (map_a2b.getColTo() - map_a2b.getColFrom() + 1):
        print "ERRONEOUS LINK: %s" % str(link)
        raise "length discrepancy"

    coverage_a = 100.0 * \
        (map_a2b.getRowTo() - map_a2b.getRowFrom() + 1) / link.mQueryLength
    coverage_b = 100.0 * \
        (map_a2b.getColTo() - map_a2b.getColFrom() + 1) / link.mSbjctLength

    # check exon boundaries, look at starts, skip first exon
    def MyMap(a, x):
        if x < a.getRowFrom():
            return 0
        while x <= a.getRowTo():
            c = a.mapRowToCol(x)
            if c:
                return c
            x += 1
        else:
            return 0

    mapped_boundaries = UniquifyList(
        map(lambda x: MyMap(map_a2b, x.mPeptideFrom / 3 + 1), cds1[1:]))
    reference_boundaries = UniquifyList(
        map(lambda x: x.mPeptideFrom / 3 + 1, cds2[1:]))

    nmissed = 0
    nfound = 0
    nmin = min(len(mapped_boundaries), len(reference_boundaries))
    nmax = max(len(mapped_boundaries), len(reference_boundaries))
    both_single_exon = len(cds1) == 1 and len(cds2) == 1
    one_single_exon = len(cds1) == 1 or len(cds2) == 1
    if len(mapped_boundaries) < len(reference_boundaries):
        mless = mapped_boundaries
        mmore = reference_boundaries
    else:
        mmore = mapped_boundaries
        mless = reference_boundaries

    # check if exon boundaries are ok
    for x in mless:
        is_ok = 0
        for c in mmore:
            if abs(x - c) < param_boundaries_max_slippage:
                is_ok = 1
                break
        if is_ok:
            nfound += 1
        else:
            nmissed += 1

    # set is_ok for dependent on exon boundaries
    # in single exon cases, require a check of coverage
    is_ok = False
    check_coverage = False
    if both_single_exon or one_single_exon:
        is_ok = True
        check_coverage = True
    else:
        if nmin == 1:
            is_ok = nmissed == 0
        elif nmin == 2:
            is_ok = nmissed <= 1
        elif nmin > 2:
            is_ok = nfound >= 2

    cc = min(coverage_a, coverage_b)

    if param_loglevel >= 3:
        print "# nquery=", len(cds1), "nsbjct=", len(cds2), "nmin=", nmin, "nmissed=", nmissed, "nfound=", nfound, \
              "is_ok=", is_ok, "check_cov=", check_coverage, \
              "min_cov=", cc, coverage_a, coverage_b, \
              "mapped=", mapped_boundaries, "reference=", reference_boundaries

    if not is_ok:
        return True, "different exon boundaries"

    if check_coverage and cc < param_min_coverage:
        return True, "low coverage"

    return False, None

# ------------------------------------------------------------


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
        elif o == "--cds-gtf-file":
            param_filename_cds = a
        elif o == "--peptides1":
            param_filename_peptides1 = a
        elif o == "--peptides2":
            param_filename_peptides2 = a
        elif o == "--transcripts1":
            param_filename_transcripts1 = a
        elif o == "--transcripts2":
            param_filename_transcripts2 = a
        elif o == "--write-exons":
            param_write_exons = a.split(",")
        elif o == "--write-introns":
            param_write_introns = a.split(",")
        elif o == "--extend-introns":
            param_extend_introns = int(a)
        elif o == "--min-intron-length":
            param_min_intron_length = int(a)
        elif o == "--max-intron-length":
            param_max_intron_length = int(a)
        elif o == "--only-best":
            param_only_best = 1
        elif o == "--boundaries-max-slippage":
            param_boundaries_max_slippage = int(a)
        elif o == "--boundaries-max-missed":
            param_boundaries_max_missed = int(a)
        elif o == "--boundaries-allow-missed":
            param_boundaries_allow_missed = int(a)
        elif o == "--missing-max-missing":
            param_missing_max_missing = int(a)
        elif o == "--missing-min-present":
            param_missing_min_present = int(a)
        elif o == "--disable-check-exon-number":
            param_check_exon_number = 0
        elif o == "--report-step":
            param_report_step = int(a)

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(2)

    print E.GetHeader()
    print E.GetParams()
    sys.stdout.flush()

    if param_loglevel >= 1:
        print "# reading exon boundaries."
        sys.stdout.flush()

    cds = Exons.ReadExonBoundaries(open(param_filename_cds, "r"))

    if param_loglevel >= 1:
        print "# read %i cds" % (len(cds))
        sys.stdout.flush()

    ninput, npairs, nskipped = 0, 0, 0

    for line in sys.stdin:
        if line[0] == "#":
            continue
        if line[0] == ">":
            print line[:-1]
            continue

        ninput += 1
        link = BlastAlignments.Link()

        link.Read(line)

        if link.mQueryToken == link.mSbjctToken:
            continue

        keep = 1
        if link.mQueryToken in cds and link.mSbjctToken in cds:
            is_paralog, reason = IsParalogLink(
                link, cds[link.mQueryToken], cds[link.mSbjctToken])
            if is_paralog:
                keep = 0
                if param_loglevel >= 2:
                    print "# DISCARDED because %s: %s" % (reason, str(link))
        else:
            if param_loglevel >= 2:
                print "# SKIPPED: %s" % str(link)
            nskipped += 1

        if keep:
            print str(link)
            npairs += 1

        if (ninput % param_report_step) == 0:
            if param_loglevel >= 1:
                print "# ninput=%i, noutput=%i, nskipped=%i" % (ninput, npairs, nskipped)
            sys.stdout.flush()

    if param_loglevel >= 1:
        print "# ninput=%i, noutput=%i, nskipped=%i" % (ninput, npairs, nskipped)

    print E.GetFooter()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
