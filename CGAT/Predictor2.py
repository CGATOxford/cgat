'''Predictor2 -
====================

Old Predictor classes. Needs to be sorted out
with Predictor.py

Used by predict_genes.py and gff2predictions.py
'''


import os
import sys
import string
import re
import tempfile
import subprocess

import CGAT.Experiment as E
import CGAT.PredictionParser as PredictionParser


class Experiment:

    mShortOptions = ""
    mLongOptions = []

    mLogLevel = 0
    mTest = 0
    mDatabaseName = None

    mName = sys.argv[0]

    def __init__(self):

        # process command-line arguments
        (self.mOptlist, self.mArgs) = self.ParseCommandLine()

        # set options now
        self.ProcessOptions(self.mOptlist)

    def DumpParameters(self):
        """dump parameters of this object. All parameters start with a
        lower-case m."""

        members = self.__dict__

        print "#" + "-" * 50
        print "#" + string.join(sys.argv)
        print "# pid: %i, system:" % os.getpid(), string.join(os.uname(), ",")
        print "#" + "-" * 50
        print "# Parameters for instance of <" + self.mName + \
            "> on " + time.asctime(time.localtime(time.time()))

        member_keys = list(members.keys())
        member_keys.sort()
        for member in member_keys:
            if member[0] == 'm':
                print "# %-40s:" % member, members[member]

        print "#" + "-" * 50
        sys.stdout.flush()

    def ProcessOptions(self, optlist):
        """Sets options in this module. Please overload as necessary."""

        for o, a in optlist:
            if o in ("-V", "--Verbose"):
                self.mLogLevel = string.atoi(a)
            elif o in ("-T", "--test"):
                self.mTest = 1

    def ProcessArguments(self, args):
        """Perform actions as given in command line arguments."""

        if self.mLogLevel >= 1:
            self.DumpParameters()

        for arg in args:
            if arg[-1] == ")":
                statement = "self.%s" % arg
            else:
                statement = "self.%s()" % arg
            exec statement

            if self.mLogLevel >= 1:
                print "-" * 50
                print statement + " finished at " + \
                    time.asctime(time.localtime(time.time()))
                print "-" * 50

    def ParseCommandLine(self):
        """Call subroutine with command line arguments."""

        self.mShortOptions = self.mShortOptions + "V:D:T"
        self.mLongOptions.append("Verbose=")
        self.mLongOptions.append("Database=")
        self.mLongOptions.append("Test")

        try:
            optlist, args = getopt.getopt(sys.argv[1:],
                                          self.mShortOptions,
                                          self.mLongOptions)
        except getopt.error, msg:
            self.PrintUsage()
            print msg
            sys.exit(2)

        return optlist, args

    def Process(self):
        self.ProcessArguments(self.mArgs)

    def PrintUsage(self):
        """print usage information."""

        print "# valid short options are:", self.mShortOptions
        print "# valid long options are:", str(self.mLongOptions)


class Predictor(Experiment):

    mLogLevel = 0
    mKeepTemp = False

    def __init__(self):
        pass

    def __call__(self,
                 query_token, peptide_sequence,
                 sbjct_token, genomic_sequence,
                 options="",
                 bracket_from=0, bracket_to=0,
                 peptide_from=0, peptide_to=0):
        """run a single prediction using exonerate.
        """

        if self.mLogLevel >= 1:
            print "# Predictor:  %i-%i" % (bracket_from, bracket_to)

        if bracket_from or bracket_to:
            genomic_sequence = genomic_sequence[bracket_from:bracket_to]
        if peptide_from or peptide_to:
            peptide_sequence = peptide_sequence[peptide_from:peptide_to]

        peptide_sequence = self.ProcessPeptideSequence(peptide_sequence)

        if self.mLogLevel >= 2:
            print "# Predictor: using peptide sequence:", peptide_sequence

        if self.mLogLevel >= 10:
            print "# Predictor using genomic sequence:", genomic_sequence

        outfile, filename_fragment = tempfile.mkstemp()
        os.write(outfile, ">%s\n%s\n" % (sbjct_token, genomic_sequence))
        os.close(outfile)

        outfile, filename_peptide = tempfile.mkstemp()
        os.write(outfile, ">%s\n%s\n" % (query_token, peptide_sequence))
        os.close(outfile)

        statement = string.join(map(str, (
            self.mExecutable,
            self.mOptions,
            self.mOutputOptions,
            options,
            "--query",
            filename_peptide,
            "--target",
            filename_fragment,
        )), " ")

        if self.mLogLevel >= 2:
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
            return None

        if self.mLogLevel >= 3:
            print"# received %i lines from %s" % (len(lines), self.mExecutable)
            print lines

        if not self.mKeepTemp:
            os.remove(filename_fragment)
            os.remove(filename_peptide)

        try:
            result = self.mParser.Parse(
                lines, peptide_sequence, genomic_sequence)
        except PredictionParser.InputError, msg:
            print "# ERROR: parsing error in input:"
            print "# received %i lines from %s" % (len(lines), self.mExecutable)
            for l in lines:
                print "#", l[:-1]
            return None

        if not result:
            return None

        # shift aligned region to fragment
        if bracket_from:
            result.ShiftGenomicRegion(bracket_from)
        if peptide_from:
            result.ShiftQueryRegion(peptide_from)

        return result

    def ProcessPeptideSequence(self, sequence):
        """perform preprocessing of the peptide sequence.

        Here: mask non-standard amino acids
        """
        return re.sub("[^ACDEFGHIKLMNPQRSTVWY]", "X", sequence.upper())

    def SetLogLevel(self, loglevel):
        """set loglevel of module."""
        old_level = self.mLogLevel
        self.mLogLevel = loglevel
        return old_level


class PredictorGenewise(Predictor):

    def __init__(self):
        Predictor.__init__(self)
        self.mParser = PredictionParser.PredictionParserGenewise()
        self.mExecutable = "genewise"
        self.mOptions = "-pseudo -init endbias"
        self.mOutputOptions = "-quiet -sum -gff -trans -pep -alb"


class PredictorExonerate (Predictor):

    def __init__(self):
        Predictor.__init__(self)
        self.mParser = PredictionParser.PredictionParserExonerate()
        self.mParser.SetRestrictToStrand("+")
        self.mExecutable = "exonerate"
        self.mOptions = "-Q protein -T dna -m p2g --softmasktarget TRUE --forcegtag TRUE --forwardcoordinates FALSE"
        self.mOutputOptions = '--showalignment FALSE --showvulgar FALSE --showsugar FALSE --showcigar FALSE ' + \
                              ' --ryo "diy\t%S\t%ql\t%r\t%pi\t%ps\t%V\n" --showtargetgff FALSE --showquerygff FALSE'
