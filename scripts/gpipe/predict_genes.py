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
gpipe/predict_genes.py - 
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

   python gpipe/predict_genes.py --help

Type::

   python gpipe/predict_genes.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import tempfile
import time
import subprocess
import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.IndexedFasta as IndexedFasta
import alignlib_lite
import CGAT.PredictionParser as PredictionParser
import CGAT.Exons as Exons

# import all, sort out names later
from CGAT.Predictor2 import *


USAGE = """python %s [OPTIONS] peptide genome

Version: $Id: gpipe/predict_genes.py 2462 2009-01-28 10:18:22Z andreas $

Wrapper for running gene predictions.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-i, --bracket-increment=        residues by how much to increase the genomic region to scan.
-b, --query-border=             residues allowed to be missing from query at either end.
-g, --genome-file=           pattern for filenames with the genomic DNA (FASTA).
-e, --exit-identical            exit, if there was no change from previous run
-m, --min-score=                minimum score
-p, --method=                   prediction method to use {genewise|exonerate}
-c, --incremental               do incremental search for genes.
-r, --recursive                 do recursive search for genes
-f, --refinement                do refinement of final prediction
-o, --probe                     do probing for gene region
--probe-options                 set probing options (for the expert)
-x, --exons-file                     filename with exon boundaries of queries.
-a, --mask-probe=               mask in probing step. Possible maskers are [seg,bias]
-f, --format                    input format
--keep-temp                     do not delete temporary files (for debugging purposes).
--graph-cutoff=                 in graph format, stop processing after this.
--peptides-fasta-file=                     filename with peptide sequences
""" % sys.argv[0]

HEADER = """# QUERY:        1  query id
# SBJCT:        2  sbjct id
# SCORE:        3  genewise alignment score
# QFROM:        4  query first residue
# QTO:          5  query last residue
# QALI:         6  query alignment
# SBJCT:        7  sbjct
# SFROM:        8  sbjct first residue
# STO:          9  sbjct last residue
# SALI:         10 sbjct alignment
# QLEN:         11 length of query
# CQUERY:       12 coverage of query (in percent)
# NGAPS:        13 number of gaps in alignment 
# NFR:          14 number of frame-shifts
# NINTRON:      15 number of introns
# NPHASE0:      16 number of phase 0 introns
# NPHASE1:      17 number of phase 1 introns
# NPHASE2:      18 number of phase 2 introns
# NSTOP:        19 number of stop codons in exons
# PIDE:         20 percent identity
# PSIM:         21 percent similarity
# PEP:          22 predicted peptide sequence
# SGFROM:       23 sbjct: genomic region first residue
# SGTO:         24 sbjct: genomic region last residue
# GALI:         25 peptide to query aligment
# NERROR:       26 number of errors in genewise parsing"""

SHORT_HEADER = """# QUERY\tSBJCT\tSCORE\tQFROM\tQTO\tSBJCT\tSFROM\tSTO\tSALI\tQLEN\tCQUERY\tNGAPS\tNFR\tNINTRON\tNPHASE0\tNPHAS1\tNPHASE2\tNSTOP\tPIDE\tPSIM\tPEP\tSGFROM\tSGTO\tGALI\tNERROR"""

global_options = {}


class Masker:

    mLogLevel = 3
    mExecutable = "biasdb.pl"
    mOptions = ""

    def __init__(self):
        pass

    def __call__(self, peptide_sequence):
        """mask peptide sequence
        """
        Masker.__init__(self)

        outfile, filename_peptide = tempfile.mkstemp()
        os.write(outfile, ">test\n%s\n" % (peptide_sequence))
        os.close(outfile)

        outfile, filename_output = tempfile.mkstemp()
        os.close(outfile)

        statement = string.join(map(str, (
            self.mExecutable,
            filename_peptide,
            self.mOptions
        )), " ")

        if self.mLogLevel >= 3:
            print "# statement: %s" % statement
            sys.stdout.flush()

        p = subprocess.Popen(statement,
                             shell=True,
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             close_fds=True)

        (file_stdout, file_stdin, file_stderr) = (p.stdout, p.stdin, p.stderr)
        lines = file_stdout.readlines()
        lines_stderr = file_stderr.readlines()
        exit_code = p.returncode

        if exit_code:
            print "# ERROR: %s returned with non-zero exit status %s" % (self.mExecutable, str(exit_code))
            for line in lines_stderr:
                print "#", line[:-1]
            sys.stdout.flush()
            return None

        os.remove(filename_peptide)

        if self.mLogLevel >= 3:
            print"# received %i lines from %s" % (len(lines), self.mExecutable)
            print lines

        masked_sequence = re.sub("\s", "", string.join(lines[1:], ""))
        return masked_sequence


class MaskerBias (Masker):

    mLogLevel = 0
    mExecutable = "biasdb.pl"
    mOptions = ""


class MaskerSeg (Masker):

    mLogLevel = 0
    mExecutable = "seg"
    mOptions = "12 2.2 2.5 -x"


# --------------------------------------------------------------------------
class TranscriptPredictor(Experiment):

    mSensitivityLevelStart = 0

    def __init__(self, filename_peptides, filename_genome):

        # border to be added at the refinement step
        self.mBorderRefinement = 50000

        # border to be added at the refinement step, if
        # Terminus has been reached
        self.mBorderRefinementSmall = 300

        # loglevel of module
        self.mLogLevel = 0

        # exit, if no change found in successive runs
        self.mExitIfIdentical = 0

        # refine predictions
        self.mDoRefinement = 0

        # run recursive refinement
        self.mDoRecursive = 0

        # do probing of alignments
        self.mDoProbe = 0

        # run incremental scanning for genes
        self.mDoIncremental = 0

        # residues which are permitted to be missing
        # at termini of the query.
        self.mQueryBorder = 0

        # mask in probing step
        self.mMaskProbe = []

        self.mRefinementMinOverlapResidues = 20

        # distance between neighbouring predictions to be combined
        self.mRefinementMinDistanceNucleotides = 100

        self.mRefinementMaxPermissiveOverlap = 10

        # availability of exon boundaries
        self.mExons = None

        # bracket in which prediction to run
        self.mBracketFrostart = 0
        self.mBracketFroend = 0

        # at least 100 residues for prediction
        self.mMinRegionLength = 100

    def SetFilenamePeptides(self, filename_peptides):
        self.mFilenamePeptides = filename_peptides

    def SetFilenameGenome(self, filename_genome):
        self.mFilenameGenome = filename_genome

    # ------------------------------------------------------------------------
    def RunIncremental(self, bracket_from, bracket_to, bracket_from_end, bracket_to_end):
        """Run predictions."""

        t0 = time.time()

        if self.mLogLevel >= 1:
            print "# INCREMENTAL: checking region: %i-%i (%i-%i)" % (self.mSbjctFrom + bracket_from,
                                                                     self.mSbjctFrom +
                                                                     bracket_to,
                                                                     self.mSbjctFrom +
                                                                     bracket_from_end,
                                                                     self.mSbjctFrom + bracket_to_end)

        last_result = None

        niterations = 0
        self.mSensitivityLevel = self.mSensitivityLevelStart
        nresults = 0

        query_coverage = 0.0

        # current increment to use
        num_increment = 0

        key = "%s_vs_%s_%s_%i_%i" % (self.mQueryToken,
                                     self.mSbjctToken, self.mSbjctStrand,
                                     bracket_from + self.mSbjctFrom, bracket_to + self.mSbjctFrom)

        last_mode = None
        last_coverage = 0

        while 1:

            t1 = time.time()

            # check, whether we have reached right/left border
            left_ok = bracket_from == bracket_from_end
            right_ok = bracket_to == bracket_to_end

            niterations += 1

            if self.mLogLevel >= 2:
                print "# INCREMENTAL: started iteration %i at %s" % (niterations, time.asctime(time.localtime(time.time())))

            result = self.RunSinglePrediction(bracket_from, bracket_to)

            if not result:
                return None, bracket_from, bracket_to

            if result.GetNumMatches() == 0:
                print "# WARNING: received empty result."

                # increase sensitivity, if there are more levels available and
                # the threshold permits it (is 0)
                if (self.mSensitivityLevel < len(self.mLevelsSensitivity) - 1) and\
                        (self.mLevelsSensitivity[self.mSensitivityLevel + 1][0] > 0):
                    if self.mLogLevel >= 2:
                        print "# Increasing sensitivity to %i." % (self.mSensitivityLevel)
                    last_mode = "sensitivity"
                    continue
                else:
                    return None, bracket_from, bracket_to

            best_entry = result.GetBestMatch()

            if self.mLogLevel >= 1:
                print "# INCREMENTAL: key=%s, iteration=%i, sensivity=%i, hits=%i, score=%5.2f, coverage=%5.2f, time=%i" % \
                      (self.mKey, niterations, self.mSensitivityLevel,
                       result.GetNumMatches(
                       ), best_entry.score, best_entry.mQueryCoverage,
                       time.time() - t1)

            ###################################################################
            # dump new results
            if last_result != result:
                last_coverage = best_entry.mQueryCoverage

                result.Write()
                last_result = result
                nresults += 1

                if last_coverage > best_entry.mQueryCoverage:
                    print "# WARNING: coverage has decreased."

            else:
                if self.mLogLevel >= 2:
                    print "# WARNING: received identical result."

            ###################################################################
            # decide, if to exit.
            if best_entry.score < self.mMinScore:
                if self.mLogLevel >= 2:
                    print "# EXITING: score below minimum"
                return None, bracket_from, bracket_to

            # stop run, if query has been found.
            if best_entry.mQueryFrom == 1 and best_entry.mQueryTo == best_entry.mQueryLength:
                if self.mLogLevel >= 2:
                    print "# EXITING: 100% coverage"
                break

            if best_entry.mQueryFrom < 1 + self.mQueryBorder and \
                    best_entry.mQueryTo > best_entry.mQueryLength - self.mQueryBorder:
                if self.mLogLevel >= 2:
                    print "# EXITING: almost complete entry."
                break

            # when bracket increment did produce no results, exit.
            if last_mode == "brackets" and \
               self.mSensitivityLevel == len(self.mLevelsSensitivity) - 1 and \
               self.mExitIfIdentical:
                if self.mLogLevel >= 2:
                    print "# EXITING: no change after bracket enlargement."
                return result, bracket_from, bracket_to

            ###############################################################
            # decide what to change from previous run.

            # increase sensitivity, if there are more levels available and
            # the threshold permits it (is larger than coverage)
            if last_mode != "sensitvity":
                last_mode = "sensitivity"
                if self.mSensitivityLevel < len(self.mLevelsSensitivity) - 1 and \
                        (self.mLevelsSensitivity[self.mSensitivityLevel + 1][0] <= best_entry.mQueryCoverage):
                    self.mSensitivityLevel += 1
                    if self.mLogLevel >= 2:
                        print "# Increasing sensitivity to %i" % (self.mSensitivityLevel)
                    continue

            # increase region to be searched.
            if last_mode != "brackets" and \
                    self.mBracketIncrements:

                if left_ok and right_ok:
                    if self.mLogLevel >= 2:
                        print "# EXITING: maximal region searched."
                    break

                # stop extension, if peptide is (almost) complete at either
                # terminus
                if best_entry.mQueryTo > best_entry.mQueryLength - self.mQueryBorder:
                    bracket_to = min(bracket_to_end,
                                     best_entry.mSbjctGenomeTo - self.mSbjctFrom + 10)
                    bracket_to_end = bracket_to

                if best_entry.mQueryFrom < 1 + self.mQueryBorder:
                    bracket_from = max(0,
                                       best_entry.mSbjctGenomeFrom - self.mSbjctFrom - 10)
                    bracket_from_end = bracket_from

                bracket_from = max(bracket_from - self.mBracketIncrements[num_increment],
                                   bracket_from_end)
                bracket_to = min(bracket_to + self.mBracketIncrements[num_increment],
                                 bracket_to_end)

                if self.mLogLevel >= 2:
                    print "# Increasing brackets by %i" % self.mBracketIncrements[num_increment]

                # change increments
                num_increment += 1
                if num_increment >= len(self.mBracketIncrements):
                    num_increment = 0

                last_mode = "brackets"
                continue

            if self.mLogLevel >= 2:
                print " EXITING: all possibilities exhausted."

            break

        if result:

            if self.mLogLevel >= 1:
                if result.GetNumMatches() > 0:
                    best_entry = result.GetBestMatch()
                    coverage = best_entry.mQueryCoverage
                    pide = best_entry.mPercentIdentity
                else:
                    coverage = 0
                    pide = 0
                print "# INCREMENTAL: key=%s, iteration=%i, sensitivity=%i, results=%i, hits=%i, coverage=%5.2f, pid=%5.2f, time=%i" % \
                      (key, niterations, self.mSensitivityLevel, nresults, result.GetNumMatches(), coverage, pide,
                       time.time() - t0)

        else:
            print "# INCREMENTAL: key=%s, iteration=%i, results=%i, hits=%i, coverage=%5.2f, pide=%5.2f, time=%i" % \
                  (key, niterations, 0, 0, 0, 0, time.time() - t0)

        return result, bracket_from, bracket_to

    # ------------------------------------------------------------------------
    def RunRecursive(self,
                     bracket_from, bracket_to,
                     bracket_from_end, bracket_to_end, level):

        if self.mLogLevel >= 1:
            print "# RECURSE: level %i: checking region: %i-%i (%i-%i)" % \
                  (level,
                   self.mSbjctFrom + bracket_from,
                   self.mSbjctFrom + bracket_to,
                   self.mSbjctFrom + bracket_from_end,
                   self.mSbjctFrom + bracket_to_end)

        if bracket_to - bracket_from <= 0:
            if self.mLogLevel >= 2:
                print "# RECURSE: exiting: bracket too small."
            return None

        result = None

        # probe for genomic region and change alignments if necessary:
        if self.mDoProbe:
            result = self.RunProbe(bracket_from, bracket_to,
                                   bracket_from_end, bracket_to_end)

            if result and result.GetNumMatches() > 0:
                match_from = result[0].mSbjctGenomeFrom
                match_to = result[0].mSbjctGenomeTo
                for x in range(1, result.GetNumMatches()):
                    match_from = min(match_from, result[x].mSbjctGenomeFrom)
                    match_to = max(match_to, result[x].mSbjctGenomeTo)

                if match_from - self.mSbjctFrom < bracket_from:
                    if self.mLogLevel >= 2:
                        print "# RECURSE: changed left border %i->%i" % (bracket_from + self.mSbjctFrom, match_from)
                    bracket_from = match_from - self.mSbjctFrom

                if match_to - self.mSbjctFrom > bracket_to:
                    if self.mLogLevel >= 2:
                        print "# RECURSE: changed right border: %i->%i" % (bracket_to + self.mSbjctFrom, match_to)
                    bracket_to = match_to - self.mSbjctFrom
            elif level == 0:
                if self.mLogLevel >= 2:
                    print "# RECURSE: continuing with dummy: first level received empty result from PROBE."
            else:
                if self.mLogLevel >= 2:
                    print "# RECURSE: exiting: received empty result from PROBE."
                return None

        if self.mDoIncremental:
            result, dummy_bracket_from, dummy_bracket_to = self.RunIncremental(bracket_from,
                                                                               bracket_to,
                                                                               bracket_from_end,
                                                                               bracket_to_end)

        if self.mDoRefinement:
            # create a dummy result, if nothing created previously.
            if not result:
                result = PredictionParser.Predictions()
                e = PredictionParser.PredictionParserEntry()
                e.mQueryToken = self.mQueryToken
                e.mSbjctToken = self.mSbjctToken
                e.mSbjctStrand = self.mSbjctStrand
                e.mSbjctGenomeFrom = bracket_from_end + self.mSbjctFrom
                e.mSbjctGenomeTo = bracket_to_end + self.mSbjctFrom
                result.Add(e)
            result = self.RunRefinement(
                result, bracket_from_end, bracket_to_end)

        if result and result.GetNumMatches() > 0 and self.mDoRecursive:

            match_from = result[0].mSbjctGenomeFrom
            match_to = result[0].mSbjctGenomeTo
            for x in range(1, result.GetNumMatches()):
                match_from = min(match_from, result[x].mSbjctGenomeFrom)
                match_to = max(match_to, result[x].mSbjctGenomeTo)
            match_from -= self.mSbjctFrom
            match_to -= self.mSbjctFrom

            # run left prediction
            if match_from > bracket_from:
                new_result = self.RunRecursive(
                    bracket_from, match_from, bracket_from_end, match_from, level + 1)
                if new_result:
                    result.Combine(new_result)

            # run right prediction
            if match_to < bracket_to:
                new_result = self.RunRecursive(
                    match_to, bracket_to, match_to, bracket_to_end, level + 1)
                if new_result:
                    result.Combine(new_result)

        return result

    # ------------------------------------------------------------------------
    def RunInitialize(self):
        """setup variables for run."""

        # read genomic sequence
        self.mForwardSequences, self.mReverseSequences = Genomics.ReadGenomicSequences(
            open(self.mFilenameGenome, "r"))

        # read peptide sequence from peptide sequences
        self.mPeptideSequences = Genomics.ReadPeptideSequences(
            open(self.mFilenamePeptides, "r"))

        # select the peptide sequence to work on.
        self.mPeptideSequence = self.mPeptideSequences[self.mQueryToken]

        self.mGenomicSequence = self.mForwardSequences[self.mSbjctToken]

        if self.mBracketToEnd == 0:
            self.mBracketToEnd = len(self.mGenomicSequence)

        if self.mBracketTo == 0:
            self.mBracketTo = self.mBracketToEnd

        # length of query
        self.mQueryLength = len(self.mPeptideSequence)

        # length of sbjct
        self.mSbjctLength = self.mBracketToEnd - self.mBracketFroend

        # key
        self.mKey = "%s_vs_%s_%s_%i_%i" % (self.mQueryToken,
                                           self.mSbjctToken, self.mSbjctStrand,
                                           self.mSbjctFrom, self.mSbjctTo)

    # ------------------------------------------------------------------------
    def Run(self):
        """Main startup routine.

        Calls iterative matching routine.
        """

        t0 = time.time()

        self.RunInitialize()

        if self.mLogLevel >= 1:
            print "# START: key=%s lquery=%i lsbjct=%i time=%s " %\
                  (self.mKey,
                   len(self.mPeptideSequences[self.mQueryToken]),
                   self.mBracketToEnd - self.mBracketFroend,
                   time.asctime(time.localtime(time.time())))

        result = self.RunRecursive(self.mBracketFrom, self.mBracketTo,
                                   self.mBracketFroend, self.mBracketToEnd,
                                   0)

        if result:

            result.Sort(lambda x, y: cmp(-x.score, -y.score))

            if len(result) > 0:
                best_entry = result.GetBestMatch()
                coverage = best_entry.mQueryCoverage
                pide = best_entry.mPercentIdentity

                if self.mLogLevel >= 2:
                    row_seq = alignlib_lite.makeSequence(
                        self.mPeptideSequences[best_entry.mQueryToken])
                    col_seq = alignlib_lite.makeSequence(
                        best_entry.mTranslation)

                    f = alignlib_lite.AlignmentFormatExplicit(
                        best_entry.mMapPeptide2Translation, row_seq, col_seq)
                    print "# TOPALI:", f.mRowAlignment
                    print "# TOPALI:", f.mColAlignment
            else:
                coverage = 0
                pide = 0

            if self.mLogLevel >= 1:
                print "# RESULT: key=%s, hits=%i, coverage=%5.2f, pide=%5.2f, time=%i" % \
                      (self.mKey, len(result), coverage, pide,
                       time.time() - t0)
        else:
            if self.mLogLevel >= 1:
                print "# RESULT: key=%s, time=%i no prediction possible." % \
                    (self.mKey, time.time() - t0)

        return result

    # ------------------------------------------------------------------------
    def RunRefinement(self, result, bracket_from_end, bracket_to_end):
        """refine result.

        The objective of this method is to return a best prediction from a list
        of results. The method stops, if a satisfactory assignment has been found.

        This is done the following way:

        Iterate over list of regions. Test predictions against exon structure.
        If the exon structure matches, return.

        Otherwise:
                1. Remove predictions adjacent in the genome but overlapping in the query.
                These are usually due to repeats. Take the one with best coverage.

                2. Combine predictions collinear in both genome and query.

                3. Combine predictions overlapping on sbjct.

        """

        t0 = time.time()

        if result.GetNumMatches() == 0:
            return result

        # sort results by score
        result.Sort(lambda x, y: cmp(-x.score, -y.score))

        boundaries_genome = []
        boundaries_peptide = []

        tmp_result = PredictionParser.Predictions()
        new_result = PredictionParser.Predictions()

        for p in result:

            is_compatible = True

            for xp in tmp_result:

                overlap = min(xp.mQueryTo, p.mQueryTo) - \
                    max(xp.mQueryFrom, p.mQueryFrom)
                goverlap = min(xp.mSbjctGenomeTo, p.mSbjctGenomeTo) - \
                    max(xp.mSbjctGenomeFrom, p.mSbjctGenomeFrom)

                if overlap > 0 or goverlap > 0:
                    is_compatible = False
                    break

            if is_compatible:
                tmp_result.append(p)

        if self.mLogLevel >= 2:
            print "# REFINE: running combination of %i predictions" % len(tmp_result)

        new_p = tmp_result.GetRange()

        if global_options.loglevel >= 2:
            print "# REFINE: predictions combined:"
            for p in tmp_result:
                print "#", str(p)

            print "#", str(new_p)

        new_prediction = self.RefinePrediction(
            new_p, bracket_from_end, bracket_to_end)
        if new_prediction:
            new_result.Add(new_prediction)

        new_result.Sort(lambda x, y: cmp((-x.score,),
                                         (-y.score,)))

        if new_result and new_result.GetNumMatches() > 0:
            if self.mLogLevel >= 1:
                if result.GetNumMatches() > 0:
                    best_entry = new_result.GetBestMatch()
                    coverage = best_entry.mQueryCoverage
                    pide = best_entry.mPercentIdentity
                else:
                    coverage = 0
                    pide = 0
                print "# REFINE: key=%s, hits=%i, coverage=%5.2f, pid=%5.2f, time=%i" % \
                      (self.mKey, new_result.GetNumMatches(), coverage, pide,
                       time.time() - t0)
                print str(new_result)
        else:
            print "# REFINE: key=%s, hits=0, time=%i" % (self.mKey, time.time() - t0)

        return new_result

    # ------------------------------------------------------------------------
    def RunProbeStep(self, level,
                     bracket_from, bracket_to, bracket_from_end, bracket_to_end,
                     peptide_from, peptide_to, peptide_from_end, peptide_to_end):
        """Run a prediction for bracket range and a peptide range."""

        if self.mLogLevel >= 1:
            print "# RunProbeStep: level %i: checking region: %i-%i (%i-%i) %i-%i (%i-%i)" % \
                  (level,
                   self.mSbjctFrom + bracket_from,
                   self.mSbjctFrom + bracket_to,
                   self.mSbjctFrom + bracket_from_end,
                   self.mSbjctFrom + bracket_to_end,
                   peptide_from, peptide_to, peptide_from_end, peptide_to_end)

        if bracket_to - bracket_from < 0:
            if self.mLogLevel >= 2:
                print "# exiting recursion bracket_to - bracket_from"
            return None

        if bracket_from < bracket_from_end < 0:
            if self.mLogLevel >= 2:
                print "# exiting recursion: bracket_from end"
            return None

        if bracket_to > bracket_to_end < 0:
            if self.mLogLevel >= 2:
                print "# exiting recursion: bracket_to end"
            return None

        if peptide_to - peptide_from < self.mQueryBorder:
            if self.mLogLevel >= 2:
                print "# exiting recursion: peptide_to - peptide_from"
            return None

        self.mSensitivityLevel = 0

        t1 = time.time()

        niterations = 0

        last_result = None

        # iterate over results.
        # rescan a region with higher sensitivity,
        # 1. if there has been no match
        while 1:

            niterations += 1
            result = self.RunSinglePrediction(bracket_from,
                                              bracket_to,
                                              peptide_from,
                                              peptide_to)

            if not result or result.GetNumMatches() == 0:

                if self.mLogLevel >= 2:
                    print "# PROBE: received empty result."

                # increase sensitivity, if there are more levels available and
                # the threshold permits it (is 0)
                if (self.mSensitivityLevel < self.mProbeMaxSensitivityLevel) and\
                        (self.mLevelsSensitivity[self.mSensitivityLevel + 1][0] >= 0):
                    self.mSensitivityLevel += 1
                    if self.mLogLevel >= 2:
                        print "# PROBE: increasing sensitivity to %i." % (self.mSensitivityLevel)
                    continue
                else:
                    break
            else:
                break

        if result and result.GetNumMatches() > 0:

            best_entry = result.GetBestMatch()
            self.mSensitivityLevelStart = max(
                self.mSensitivityLevelStart, self.mSensitivityLevel)

            if self.mLogLevel >= 1:
                print "# key=%s, iteration=%i, sensivity=%i, hits=%i, score=%5.2f, coverage=%5.2f, time=%i" % \
                      (self.mKey, niterations, self.mSensitivityLevel,
                       result.GetNumMatches(
                       ), best_entry.score, best_entry.mQueryCoverage,
                       time.time() - t1)

            best_entry = result.GetBestMatch()
            match_genome_from = best_entry.mSbjctGenomeFrom - self.mSbjctFrom
            match_genome_to = best_entry.mSbjctGenomeTo - self.mSbjctFrom
            match_peptide_from = best_entry.mQueryFrom - 1
            match_peptide_to = best_entry.mQueryTo

            if global_options.loglevel >= 2:
                print str(result)

            increment = self.mBracketIncrements[
                min(level, len(self.mBracketIncrements) - 1)]

            # run left prediction, if no internal match

            self.mMaxPeptideTo = max(self.mMaxPeptideTo, match_peptide_to)
            self.mMinPeptideFrom = min(
                self.mMinPeptideFrom, match_peptide_from)

            if match_peptide_from <= self.mMinPeptideFrom:
                new_result = self.RunProbeStep(level + 1,
                                               bracket_from -
                                               increment, match_genome_from,
                                               bracket_from_end, match_genome_from,
                                               peptide_from, match_peptide_from,
                                               peptide_from_end, match_peptide_from)
                if new_result:
                    result.Combine(new_result)

            if match_peptide_to >= self.mMaxPeptideTo:
                new_result = self.RunProbeStep(level + 1,
                                               match_genome_to, bracket_to +
                                               increment, match_genome_to, bracket_to_end,
                                               match_peptide_to, peptide_to, match_peptide_to, peptide_to_end)
                if new_result:
                    result.Combine(new_result)

        return result

    # -------------------------------------------------------------------------
    def RunProbe(self, bracket_from, bracket_to, bracket_from_end, bracket_to_end):
        """Run a probe
        This tries to run the predictor in probe mode in order to pin down the location
        of the gene more accurately."""

        t0 = time.time()

        if self.mLogLevel >= 1:
            print "# PROBE: checking region: %i-%i (%i-%i)" % (self.mSbjctFrom + bracket_from,
                                                               self.mSbjctFrom +
                                                               bracket_to,
                                                               self.mSbjctFrom +
                                                               bracket_from_end,
                                                               self.mSbjctFrom + bracket_to_end)

        # maximal range of peptide found
        self.mMinPeptideFrom = self.mQueryLength
        self.mMaxPeptideTo = 0

        result = self.RunProbeStep(0, bracket_from, bracket_to, bracket_from_end, bracket_to_end,
                                   0, query_length, 0, query_length)

        if result:
            query_length = result[0].mQueryLength
            coverage = 100 * \
                (self.mMaxPeptideTo - self.mMinPeptideFrom) / self.mQueryLength
        else:
            coverage = 0

        if self.mLogLevel >= 2:
            print "# PROBE: key=%s, sensitivity=%i, hits=%i, coverage=%i, from=%i, to=%i, time=%i" % \
                  (self.mKey, self.mSensitivityLevelStart, result.GetNumMatches(),
                   coverage,
                   self.mMinPeptideFrom + 1, self.mMaxPeptideTo,
                   time.time() - t0)

        return result

# --------------------------------------------------------------------------


class TranscriptPredictorTwoStep(TranscriptPredictor):

    """Transcript predictor that
    1. probes using Method 1
    2. predicts using Method 2
    """

    mProbePredictor = None
    mProbeOptions = None

    mRefinementPredictor = None
    mRefinementOptions = None

    def __init__(self, filename_genome, filename_peptides):

        TranscriptPredictor.__init__(self, filename_genome, filename_peptides)

    # -------------------------------------------------------------------------
    def DumpParameters(self):
        TranscriptPredictor.DumpParameters(self)
        self.mProbePredictor.DumpParameters()
        self.mRefinementPredictor.DumpParameters()

    # -------------------------------------------------------------------------
    def RunInitialize(self):
        """init."""

        TranscriptPredictor.RunInitialize(self)
        self.mProbeSequence = self.mPeptideSequence

        self.mProbePredictor.SetLogLevel(self.mLogLevel - 2)
        self.mRefinementPredictor.SetLogLevel(self.mLogLevel - 2)

        if self.mMaskProbe:
            for masker in self.mMaskProbe:
                self.mProbeSequence = masker(self.mProbeSequence)
            if not self.mProbeSequence:
                print "# WARNING: empty sequence after masking with %s, query was %s" % (str(masker), self.mPeptideSequence)
            if self.mLogLevel >= 4:
                print "# PROBE: sequence=%s" % self.mProbeSequence

    # -------------------------------------------------------------------------
    def RunProbe(self, bracket_from, bracket_to, bracket_from_end, bracket_to_end):
        """Run a probe
        This tries to run the predictor in probe mode in order to pin down the location
        of the gene more accurately."""

        t0 = time.time()

        if self.mLogLevel >= 1:
            print "# PROBE: checking region: %i-%i (size=%i)" % (self.mSbjctFrom + bracket_from_end,
                                                                 self.mSbjctFrom +
                                                                 bracket_to_end,
                                                                 bracket_to_end - bracket_from_end)

        result = self.mProbePredictor(self.mQueryToken, self.mProbeSequence,
                                      self.mSbjctToken, self.mGenomicSequence,
                                      self.mProbeOptions,
                                      bracket_from_end, bracket_to_end)

        if not result:
            print "# PROBE: received empty result."
            return None

        result.ShiftGenomicRegion(self.mSbjctFrom)
        result.SetStrand(self.mSbjctStrand)

        if result and result.GetNumMatches() > 0:
            match_from = result[0].mQueryFrom
            match_to = result[0].mQueryTo
            for x in range(1, result.GetNumMatches()):
                match_from = min(match_from, result[x].mQueryFrom)
                match_to = max(match_to, result[x].mQueryTo)

            coverage = 100 * (match_to - match_from + 1) / self.mQueryLength

            if self.mLogLevel >= 2:
                print "# PROBE: key=%s, hits=%i, coverage=%i, from=%i, to=%i, time=%i" % \
                      (self.mKey,
                       result.GetNumMatches(),
                       coverage,
                       match_from, match_to,
                       time.time() - t0)
                print str(result)
        else:
            print "# PROBE: key=%s, hits=0" % \
                  (self.mKey)

        return result

    # ------------------------------------------------------------------------
    def RefinePrediction(self, prediction, bracket_from_end, bracket_to_end):
        """refine a prediction.
        """

        # set region to check based on scanning result
        #
        # if not touching Terminus:
        #
        # 1. add fixed width - large boundary
        #
        # 2. increase refinement boundary according to query
        # gene size
        #
        #              ------           match
        # +++++++++++++           query size
        # ++++++++++++     query size
        # ....................    region to test
        #
        # if touching Terminus:
        #
        # 1. just add small range.

        bracket_from = prediction.mSbjctGenomeFrom - self.mSbjctFrom
        bracket_to = prediction.mSbjctGenomeTo - self.mSbjctFrom

        prediction_length = prediction.mSbjctGenomeTo - \
            prediction.mSbjctGenomeFrom

        if self.mExons:
            genesize = self.mExons[-1].mGenomeTo - self.mExons[0].mGenomeFrom
            gene_increment = genesize - prediction_length
        else:
            genesize = 0
            gene_increment = 0

        if prediction.mQueryFrom > 1 + self.mQueryBorder:
            left_increment = max(gene_increment, self.mBorderRefinement)
        else:
            left_increment = self.mBorderRefinementSmall

        if prediction.mQueryTo < prediction.mQueryLength - self.mQueryBorder:
            right_increment = max(gene_increment, self.mBorderRefinement)
        else:
            right_increment = self.mBorderRefinementSmall

        if self.mLogLevel >= 2:
            print "# REFINE: increments are: left=%i, right=%i, lpred=%i, lquery=%i" % (left_increment,
                                                                                        right_increment,
                                                                                        prediction_length,
                                                                                        genesize)

        bracket_from -= left_increment
        bracket_to += right_increment

        # make sure that range is within bracket
        bracket_from = max(bracket_from, bracket_from_end)
        bracket_to = min(bracket_to, bracket_to_end)

        if self.mLogLevel >= 2:
            print "# REFINE: searching region %i to %i, size=%i" % (bracket_from + self.mSbjctFrom,
                                                                    bracket_to +
                                                                    self.mSbjctFrom,
                                                                    bracket_to - bracket_from)
        if bracket_to - bracket_from < self.mMinRegionLength:
            print "# REFINE: region with %i residues to small, no refinement done" % (bracket_to - bracket_from)
            return None

        result = self.mRefinementPredictor(self.mQueryToken, self.mPeptideSequences[self.mQueryToken],
                                           self.mSbjctToken, self.mForwardSequences[
                                               self.mSbjctToken],
                                           self.mRefinementOptions,
                                           bracket_from, bracket_to)

        if not result:
            print "# WARNING: received no result, no refinement done."
            return None

        result.ShiftGenomicRegion(self.mSbjctFrom)
        result.SetStrand(self.mSbjctStrand)

        if result.GetNumMatches() == 0:
            print "# REFINE: received empty result."
            return None

        return result.GetBestMatch()

# --------------------------------------------------------------------------


class TranscriptPredictorTwoStepEG(TranscriptPredictorTwoStep):

    """Transcript predictor that
    1. probes using exonerate
    2. predicts using genewise
    """

    def __init__(self, filename_genome, filename_peptides):
        TranscriptPredictorTwoStep.__init__(
            self, filename_genome, filename_peptides)

        self.mProbeOptions = "--proteinwordlimit 5 --proteinhspdropoff 5  --proteinwordlen 3 --subopt TRUE"
        self.mProbePredictor = PredictorExonerate()

        self.mRefinementOptions = ""
        self.mRefinementPredictor = PredictorGenewise()

# --------------------------------------------------------------------------


class TranscriptPredictorTwoStepEE(TranscriptPredictorTwoStep):

    """Transcript predictor that
    1. probes using exonerate
    2. predicts using exonerate
    """

    def __init__(self, filename_genome, filename_peptides,
                 min_score=80):
        TranscriptPredictorTwoStep.__init__(
            self, filename_genome, filename_peptides)

        self.mProbeOptions = "--proteinwordlimit 5 --proteinhspdropoff 5  --proteinwordlen 3 --subopt TRUE --score '%s' " %\
                             str(min_score)
        self.mProbePredictor = PredictorExonerate()

        self.mRefinementOptions = "--exhaustive --subopt FALSE --score '%s' " % str(
            min_score)
        self.mRefinementPredictor = PredictorExonerate()

# ------------------------------------------------------------------------


def EvaluatePrediction(prediction, query_exons, query_sequence):
    """Evaluate the result of a gene prediction.

    If it is a successfull gene prediction, return 1, otherwise
    return 0.

    A gene prediction is successful, if the exon boundaries in the
    query are similar to the exon boundaries in the prediction.
    """

    if global_options.loglevel >= 2:
        print "# EVAL: Exons in query:"
        i = 0
        for e in query_exons:
            print "# EVAL: %i" % i, str(e)
            i += 1

    exons = Exons.Alignment2Exons(prediction.mMapPeptide2Genome,
                                  query_from=prediction.mQueryFrom - 1,
                                  sbjct_from=prediction.mSbjctGenomeFrom,
                                  add_stop_codon=1)

    for e in exons:
        e.mQueryToken = "prediction"

    if global_options.loglevel >= 2:
        print "# EVAL: Exons in prediction:"
        i = 0
        for e in exons:
            print "# EVAL: %i" % i, str(e)
            i += 1

    comparison = Exons.CompareGeneStructures(
        exons, query_exons,
        map_ref2cmp=prediction.mMapPeptide2Translation,
        cmp_sequence=prediction.mTranslation,
        ref_sequence=query_sequence,
        threshold_min_pide=prediction.mPercentIdentity *
        global_options.quality_threshold_pide / 100,
        threshold_slipping_exon_boundary=global_options.quality_slipping_exon_boundary)

    if global_options.loglevel >= 2:
        print comparison.Pretty(prefix="# EVAL: ")

    is_ok = False
    status = "unknown"

    max_nexons = max(len(exons), len(query_exons))

    # more than two exons in result: check number of identical exons
    if len(exons) > 2:
        if comparison.mNumIdenticalExons >= (len(exons) - 2) and \
                abs(comparison.mNumDifferenceExons) <= 1:
            status = "conserved multi exon"
            is_ok = True
        elif max_nexons > 10 and \
                (100 * comparison.mNumIdenticalExons) / max_nexons > global_options.evaluate_min_percent_exon_identity:
            status = "semi-conserved multi exon"
            is_ok = True
    # two exons in result:
    elif len(exons) == 2:
        # accept if two exons in query and at least one exon is identical
        if len(query_exons) == 2:
            if comparison.mNumIdenticalExons >= 1:
                status = "conserved double exon"
                is_ok = True
        else:
            # accept if both exons are identical
            if comparison.mNumIdenticalExons == 2:
                status = "semi-conserved double exon"
                is_ok = True
    # single exon in result: accept, if query is single exon and coverage
    # above threshold
    else:
        if len(query_exons) == 1 and \
           prediction.mQueryCoverage >= global_options.evaluate_single_exons_min_coverage:
            status = "conserved single exon"
            is_ok = True

    if global_options.loglevel >= 1:
        print "# EVAL: status=%s is_ok=%s" % (status, str(is_ok))

    return is_ok

# ------------------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gpipe/predict_genes.py 2462 2009-01-28 10:18:22Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-b", "--query-border", dest="query_border", type="int")
    parser.add_option(
        "-i", "--bracket-increment", dest="bracket_increments", type="string")
    parser.add_option(
        "-e", "--exit-identical", dest="exit_if_identical", action="store_true")
    parser.add_option("-m", "--min-score", dest="min_score", type="float")
    parser.add_option("-p", "--method", dest="method", type="string")
    parser.add_option(
        "-r", "--recursive", dest="recursive", action="store_true")
    parser.add_option("--refinement", dest="refinement", action="store_true")
    parser.add_option("--probe", dest="probe", action="store_true")
    parser.add_option("--incremental", dest="incremental", action="store_true")
    parser.add_option(
        "--border-refinement", dest="border_refinement", type="int")
    parser.add_option("-x", "--exons-file", dest="filename_exons", type="string")
    parser.add_option("-a", "--mask-probe", dest="mask_probe", type="string")
    parser.add_option("-f", "--format", dest="input_format", type="string")
    parser.add_option("--probe-options", dest="probe_options", type="string")
    parser.add_option("-g", "--genome-file", dest="genome_file", type="string")
    parser.add_option("--peptides-fasta-file", dest="filename_peptides", type="string")
    parser.add_option("--keep-temp", dest="keep_temp", action="store_true")

    parser.set_defaults(
        exit_if_identical=None,
        loglevel=2,
        query_border=0,
        short_options="v:hi:g:b:em:procx:af:",
        filename_peptides=None,
        filename_genome=None,
        filename_exons=None,
        genome_file="genome_%s.fasta",
        # bracket extension parameters
        bracket_increments=None,
        # method to use for prediction
        method="genewise",
        use_header=1,
        recursive=0,
        probe=0,
        incremental=0,
        border_refinement=50000,
        min_score=0,
        refinement=0,
        mask_probe="",
        evaluate_single_exons_min_coverage=80,
        evaluate_min_percent_exon_identity=70,
        quality_threshold_pide=75,
        quality_slipping_exon_boundary=9,
        graph_cutoff=5,
        input_format=None,
        probe_options=None,
        keep_temp=False,
    )

    (global_options, args) = E.Start(parser)

    if global_options.mask_probe:
        global_options.mask_probe = global_options.mask_probe.split(",")
    if global_options.bracket_increments:
        global_options.bracket_increments = map(
            int, global_options.bracket_increments.split(","))

    if global_options.method == "twostep_eg":
        predictor = TranscriptPredictorTwoStepEG(
            global_options.filename_peptides, global_options.filename_genome)
    elif global_options.method == "twostep_ee":
        predictor = TranscriptPredictorTwoStepEE(
            global_options.filename_peptides, global_options.filename_genome, global_options.min_score)
    else:
        print "unknown method", global_options.method
        sys.exit(1)

    if global_options.probe_options:
        predictor.mProbeOptions = global_options.probe_options

    predictor.mQueryBorder = global_options.query_border
    predictor.mBorderRefinement = global_options.border_refinement
    predictor.mBracketIncrements = global_options.bracket_increments
    predictor.mMinScore = global_options.min_score
    predictor.mExitIfIdentical = global_options.exit_if_identical
    predictor.mLogLevel = global_options.loglevel
    predictor.mDoRecursive = global_options.recursive
    predictor.mDoRefinement = global_options.refinement
    predictor.mDoProbe = global_options.probe
    predictor.mDoIncremental = global_options.incremental
    predictor.mKeepTemp = global_options.keep_temp

    for m in global_options.mask_probe:
        if m == "seg":
            predictor.mMaskProbe.append(MaskerSeg())
        elif m == "bias":
            predictor.mMaskProbe.append(MaskerBias())
        else:
            raise "unknown masker %s" % m

    exit_code = 0

    print E.GetHeader()
    print E.GetParams()
    predictor.DumpParameters()

    if global_options.filename_exons:
        file = open(global_options.filename_exons, "r")
        exons = Exons.ReadExonBoundaries(file, do_invert=1, remove_utr=1)
        file.close()
    else:
        exons = {}

    if global_options.loglevel >= 2:
        print "# read exon boundaries for %i queries" % len(exons)

    fasta = IndexedFasta.IndexedFasta(global_options.genome_file)

    # --------------------------------------------------------------------
    # Format "pairs": process peptide and genomic sequences given in fasta format.
    # Query and genomic sequence are linked by common identifier
    if global_options.input_format == "pairs":
        if len(args) == 2:
            global_options.filename_peptides, global_options.filename_genome = args

        peptide_sequences = Genomics.ReadPeptideSequences(
            open(global_options.filename_peptides, "r"))
        forward_sequences = Genomics.ReadGenomicSequences(
            open(global_options.filename_genome, "r"), do_reverse=0)

        for query_token in peptide_sequences:
            if query_token not in peptide_sequences:
                print "# WARNING: no genomic sequence found "\
                    "for query %s" % query_token
                continue

            query_sequence = peptide_sequences[query_token]
            sbjct_sequence = forward_sequences[query_token]
            predictor.mBracketFrom = 0
            predictor.mBracketTo = len(sbjct_sequence)
            predictor.mQueryToken = query_token
            predictor.mSbjctToken = query_token
            predictor.mSbjctStrand = "+"
            predictor.mSbjctFrom = 0
            predictor.mSbjctTo = len(sbjct_sequence)
            predictor.mBracketFroend = 0
            predictor.mBracketToEnd = len(sbjct_sequence)

            # create temporary files
            query_outfile, query_filename = tempfile.mkstemp()
            os.write(query_outfile, ">%s\n%s\n" %
                     (query_token, query_sequence))
            os.close(query_outfile)
            predictor.SetFilenamePeptides(query_filename)

            sbjct_outfile, sbjct_filename = tempfile.mkstemp()
            os.write(sbjct_outfile, ">%s\n%s\n" %
                     (query_token, sbjct_sequence))
            os.close(sbjct_outfile)
            predictor.SetFilenameGenome(sbjct_filename)

            predictor.Run()

            os.remove(query_filename)
            os.remove(sbjct_filename)

    elif global_options.input_format == "single":

        query_token, sbjct_token, sbjct_strand, sbjct_from, sbjct_to = args[:5]

        sbjct_from, sbjct_to = int(sbjct_from), int(sbjct_to)
        peptide_sequences = Genomics.ReadPeptideSequences(
            open(global_options.filename_peptides, "r"))
        query_sequence = peptide_sequences[query_token]
        sbjct_sequence = fasta.getGenomicSequence(sbjct_token, sbjct_strand,
                                                  sbjct_from, sbjct_to)

        predictor.mBracketFrom = 0
        predictor.mBracketTo = len(sbjct_sequence)
        predictor.mQueryToken = query_token
        predictor.mSbjctToken = sbjct_token
        predictor.mSbjctStrand = sbjct_strand
        predictor.mSbjctFrom = sbjct_from
        predictor.mSbjctTo = sbjct_to
        predictor.mBracketFroend = 0
        predictor.mBracketToEnd = len(sbjct_sequence)

        print ">%s\n%s" % (query_token, query_sequence)
        print ">%s\n%s" % (sbjct_token, sbjct_sequence)

        if query_token in exons:
            predictor.mExons = exons[query_token]
        else:
            predictor.mExons = []

        # create temporary files
        query_outfile, query_filename = tempfile.mkstemp()
        os.write(query_outfile, ">%s\n%s\n" % (query_token, query_sequence))
        os.close(query_outfile)
        predictor.SetFilenamePeptides(query_filename)

        sbjct_outfile, sbjct_filename = tempfile.mkstemp()
        os.write(sbjct_outfile, ">%s\n%s\n" % (sbjct_token, sbjct_sequence))
        os.close(sbjct_outfile)
        predictor.SetFilenameGenome(sbjct_filename)

        result = predictor.Run()

        # dump result and check, if it is satisfactory
        if result:
            result.Write()

        os.remove(query_filename)
        os.remove(sbjct_filename)

    # --------------------------------------------------------------------
    # Format "lists": given are chunks of priority lists.
    # After very prediction, the prediction is evaluated. If it is sufficient,
    # the remaining entries are skipped and the next list is processed.
    elif global_options.input_format == "graph":

        # get genomes in input set (so that for single genome files,
        # unnecessary entries can be skipped).
        data = map(lambda x: x[:-1].split("\t"),
                   filter(lambda x: x[0] != "#", sys.stdin.readlines()))
        sbjct_tokens = {}
        for g in map(lambda x: x[2], data):
            sbjct_tokens[g] = True

        if len(args) == 1 and args[0] == "-":
            # process several chunks of data
            skip = False
            last_time = None
            last_region_id = None
            last_sbjct_sequence = None
            forward_sequences = None
            reverse_sequences = None

            for d in data:
                (query_token,
                 query_sequence,
                 sbjct_token,
                 sbjct_strand,
                 sbjct_sequence,
                 sbjct_from,
                 sbjct_to,
                 min_bracket_from, min_bracket_to,
                 region_id, region_nr, region_max_nr) = d

                (sbjct_from, sbjct_to, min_bracket_from, min_bracket_to, region_id, region_nr, region_max_nr) = \
                    map(int, (sbjct_from, sbjct_to, min_bracket_from,
                        min_bracket_to, region_id, region_nr, region_max_nr))

                if sbjct_sequence == "":
                    if last_sbjct_sequence is None:
                        try:
                            sbjct_sequence = fasta.getSequence(
                                sbjct_token, sbjct_strand, sbjct_from, sbjct_to)
                        except AssertionError:
                            global_options.stderr.write("# WARNING: could not retrieve sequence for in region %i-%i: %s:%s:%i:%i - skipped\n" %
                                                        (region_id, region_nr,
                                                         sbjct_token, sbjct_strand, sbjct_from, sbjct_to))
                            global_options.stdlog.write("# WARNING: could not retrieve sequence for in region %i-%i: %s:%s:%i:%i - skipped\n" %
                                                        (region_id, region_nr,
                                                         sbjct_token, sbjct_strand, sbjct_from, sbjct_to))
                            continue
                    else:
                        sbjct_sequence = last_sbjct_sequence
                else:
                    last_sbjct_sequence = sbjct_sequence

                # do not test on region_nr, as first region_nr might not
                # be 1 due to duplicated key removal in
                # gpipe/assignments2pairs.py
                if region_id != last_region_id:
                    this_time = time.time()
                    if global_options.loglevel >= 1 and last_time:
                        print "## GRAPH: region %i: finished in %i seconds" % (last_region_id, this_time - last_time)
                        print "####################################################################"
                    last_time = this_time
                    last_region_id = region_id
                    if global_options.loglevel >= 1:
                        print "####################################################################"
                        print "## GRAPH: region %i: starting with %i members" % (region_id, region_max_nr)
                    skip = False

                if global_options.loglevel >= 2:
                    print "## GRAPH: region %i: processing %i of %i members" % (region_id, region_nr, region_max_nr)

                if skip and region_nr <= region_max_nr:
                    if global_options.loglevel >= 2:
                        print "## GRAPH: skipping entry %i/%i" % (region_nr, region_max_nr)
                    continue

                if global_options.graph_cutoff and region_nr > global_options.graph_cutoff:
                    if global_options.loglevel >= 2:
                        print "## GRAPH: omitting entry %i/%i" % (region_nr, region_max_nr)
                    continue

                (bracket_from, bracket_to,
                 bracket_from_end, bracket_to_end) = map(lambda x: x - sbjct_from,
                                                         (min_bracket_from,
                                                          min_bracket_to,
                                                          sbjct_from,
                                                          sbjct_to))

                predictor.mBracketFrom = bracket_from
                predictor.mBracketTo = bracket_to
                predictor.mQueryToken = query_token
                predictor.mSbjctToken = sbjct_token
                predictor.mSbjctStrand = sbjct_strand
                predictor.mSbjctFrom = sbjct_from
                predictor.mSbjctTo = sbjct_to
                predictor.mBracketFroend = bracket_from_end
                predictor.mBracketToEnd = bracket_to_end

                if query_token in exons:
                    predictor.mExons = exons[query_token]
                else:
                    predictor.mExons = []

                # create temporary files
                query_outfile, query_filename = tempfile.mkstemp()
                os.write(query_outfile, ">%s\n%s\n" %
                         (query_token, query_sequence))
                os.close(query_outfile)
                predictor.SetFilenamePeptides(query_filename)

                sbjct_outfile, sbjct_filename = tempfile.mkstemp()
                os.write(sbjct_outfile, ">%s\n%s\n" %
                         (sbjct_token, sbjct_sequence))
                os.close(sbjct_outfile)
                predictor.SetFilenameGenome(sbjct_filename)

                result = predictor.Run()

                # dump result and check, if it is satisfactory
                if result:
                    result.Write()

                if result and query_token in exons:
                    skip = False
                    if global_options.loglevel >= 1:
                        print "# EVAL: evaluating %i predictions." % len(result)
                    for r in result:
                        skip = skip or EvaluatePrediction(
                            r, exons[query_token], query_sequence)

                os.remove(query_filename)
                os.remove(sbjct_filename)

            this_time = time.time()
            if global_options.loglevel >= 1 and last_time:
                print "## GRAPH: region %i: finished in %i seconds" % (last_region_id, this_time - last_time)
                print "####################################################################"

    # --------------------------------------------------------------------
    # process default format
    else:
        # --------------------------------------------------------------------
        # two arguments: a single prediction with one peptide file and one
        # genome file.
        if len(args) == 2:
            global_options.filename_peptides, global_options.filename_genome = args

            sbjct_token = ""
            sbjct_strand = ""
            sbjct_from = 0
            sbjct_to = 0
            query_token = ""
            bracket_from = 0
            bracket_to = 0
            bracket_from_end = 0
            bracket_to_end = 0

            if global_options.use_header:
                # read sbjct token and get bracket information
                line = open(global_options.filename_genome, "r").readline()

                data = re.split("\s+", line[1:-1])

                if len(data) == 1:
                    sbjct_token = data[0]

                elif len(data) >= 6:
                    (sbjct_token, sbjct_strand,
                     sbjct_from, sbjct_to,
                     min_bracket_from, min_bracket_to) = data[:6]

                    (sbjct_from, sbjct_to, min_bracket_from, min_bracket_to) = \
                        map(int, (
                            sbjct_from, sbjct_to, min_bracket_from, min_bracket_to))

                    # convert to relative coordinates on genomic sequence
                    (bracket_from, bracket_to,
                     bracket_from_end, bracket_to_end) = map(lambda x: x - sbjct_from,
                                                             (min_bracket_from,
                                                              min_bracket_to,
                                                              sbjct_from,
                                                              sbjct_to))

                # read query token
                line = open(global_options.filename_peptides, "r").readline()
                (query_token,) = re.split("\s+", line[1:-1])[:1]

            predictor.mBracketFrom = bracket_from
            predictor.mBracketTo = bracket_to
            predictor.mQueryToken = query_token
            predictor.mSbjctToken = sbjct_token
            predictor.mSbjctStrand = sbjct_strand
            predictor.mSbjctFrom = sbjct_from
            predictor.mSbjctTo = sbjct_to
            predictor.mBracketFroend = bracket_from_end
            predictor.mBracketToEnd = bracket_to_end
            predictor.SetFilenamePeptides(global_options.filename_peptides)
            predictor.SetFilenameGenome(global_options.filename_genome)

            if query_token in exons:
                predictor.mExons = exons[query_token]
            else:
                predictor.mExons = []

            predictor.Run()
        # --------------------------------------------------------------------
        # one argument, which is -: read input as chunks from stdin
        elif len(args) == 1 and args[0] == "-":

            # get genomes in input set (so that for single genome files,
            # unnecessary entries can be skipped).
            data = map(
                lambda x: x[:-1].split("\t"), filter(lambda x: x[0] != "#", sys.stdin.readlines()))
            sbjct_tokens = {}
            for g in map(lambda x: x[2], data):
                sbjct_tokens[g] = True

            last_sbjct_sequence = None

            # process several chunks of data
            for d in data:

                (query_token,
                 query_sequence,
                 sbjct_token,
                 sbjct_strand,
                 sbjct_sequence,
                 sbjct_from,
                 sbjct_to,
                 min_bracket_from,
                 min_bracket_to) = d

                (sbjct_from, sbjct_to, min_bracket_from, min_bracket_to) = \
                    map(int, (
                        sbjct_from, sbjct_to,
                        min_bracket_from, min_bracket_to))

                (bracket_from, bracket_to,
                 bracket_from_end, bracket_to_end) = map(
                     lambda x: x - sbjct_from,
                     (min_bracket_from,
                      min_bracket_to,
                      sbjct_from,
                      sbjct_to))

                if sbjct_sequence == "":
                    if last_sbjct_sequence is None:
                        sbjct_sequence = fasta.getSequence(sbjct_token,
                                                           sbjct_strand,
                                                           sbjct_from,
                                                           sbjct_to)
                    else:
                        sbjct_sequence = last_sbjct_sequence
                else:
                    last_sbjct_sequence = sbjct_sequence

                predictor.mBracketFrom = bracket_from
                predictor.mBracketTo = bracket_to
                predictor.mQueryToken = query_token
                predictor.mSbjctToken = sbjct_token
                predictor.mSbjctStrand = sbjct_strand
                predictor.mSbjctFrom = sbjct_from
                predictor.mSbjctTo = sbjct_to
                predictor.mBracketFroend = bracket_from_end
                predictor.mBracketToEnd = bracket_to_end

                if query_token in exons:
                    predictor.mExons = exons[query_token]
                else:
                    predictor.mExons = []

                # create temporary files
                query_outfile, query_filename = tempfile.mkstemp()
                os.write(query_outfile, ">%s\n%s\n" %
                         (query_token, query_sequence))
                os.close(query_outfile)
                predictor.SetFilenamePeptides(query_filename)

                sbjct_outfile, sbjct_filename = tempfile.mkstemp()
                os.write(sbjct_outfile, ">%s\n%s\n" %
                         (sbjct_token, sbjct_sequence))
                os.close(sbjct_outfile)
                predictor.SetFilenameGenome(sbjct_filename)

                result = predictor.Run()

                if result:
                    result.Write()

                os.remove(query_filename)
                os.remove(sbjct_filename)

    E.Stop()
    sys.exit(exit_code)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
