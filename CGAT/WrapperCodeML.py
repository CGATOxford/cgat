'''
WrapperCodeML.py -
======================================================

:Tags: Python

Code
----

'''
import sys
import string
import re
import tempfile
import subprocess
import shutil
import random
import traceback
from types import *
from CGAT import Experiment as Experiment
from CGAT import TreeTools as TreeTools
from CGAT import IOTools as IOTools
from CGAT import Tree as Tree
from CGAT import Mali as Mali
from CGAT import Genomics as Genomics
from CGAT import RateEstimation as RateEstimation


class Error(Exception):

    """Base class for exceptions in this module."""

    def __str__(self):
        return str(self.message)

    def _get_message(self):
        return self._message

    def _set_message(self, message):
        self._message = message
    message = property(_get_message, _set_message)


class ParsingError(Error):

    """Exception raised for errors while parsing

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message, line=None):
        if line:
            self.message = message + " at line " + line
        else:
            self.message = message


class UsageError(Error):

    """Exception raised for errors while starting

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


class CodeMLBranchInfo:

    """result with branch information."""

    def __init__(self, branch1, branch2, kaks, ka, ks, ndn, sds, n, s):

        self.mBranch1 = branch1
        self.mBranch2 = branch2
        self.mKaks = kaks
        self.mKa = ka
        self.mKs = ks
        self.mSds = sds
        self.mNdn = ndn
        self.mN = n
        self.mS = s

    def __str__(self):
        return "\t".join((self.mBranch1, self.mBranch2,
                          str(self.mKaks), str(self.mKa), str(self.mKs),
                          str(self.mNdn), str(self.mSds)))


class CodeMLResult:

    mMaxLineLength = 100

    def __init__(self):
        self.mWarnings = []
        self.mCheckConvergence = False

    def truncateLine(self, line):
        if len(line) > self.mMaxLineLength:
            return "%s..." % line[:self.mMaxLineLength]
        else:
            return line

    def mapNames(self, map_old2new):
        """map all names."""

        self.mTreeKs.relabel(map_old2new)
        self.mTreeKa.relabel(map_old2new)
        self.mTreeKaks.relabel(map_old2new)
        self.mTreeSds.relabel(map_old2new)
        self.mTreeNdn.relabel(map_old2new)
        self.mTreeS.relabel(map_old2new)
        self.mTreeN.relabel(map_old2new)

    def __str__(self):

        s = []
        members = self.__dict__
        keys = list(members.keys())
        keys.sort()

        for key in keys:

            if key[0] == 'm':

                m = members[key]
                if type(m) in (ListType, TupleType):
                    for x in range(len(m)):
                        s.append(
                            ("%-40s: %s" % (key,
                                            self.truncateLine(str(m[x])))))
                elif type(m) in (DictType,):
                    for x, y in list(m.items()):
                        s.append(
                            ("%-40s: %s: %s" % (key, str(x),
                                                self.truncateLine(str(y)))))
                else:
                    s.append(("%-40s: %s" % (key, self.truncateLine(str(m)))))

        return "\n".join(s)

    def buildTreesFromBranchInfo(self):
        """build trees from branch info.

        assumes that codeml prints out parent -> child
        """
        # t is a map of parent to children
        self.mTreeKs = TreeTools.Graph2Tree(
            [(x.mBranch1, x.mBranch2, x.mKs) for x in self.mBranchInfo])
        self.mTreeKa = TreeTools.Graph2Tree(
            [(x.mBranch1, x.mBranch2, x.mKa) for x in self.mBranchInfo])
        self.mTreeKaks = TreeTools.Graph2Tree(
            [(x.mBranch1, x.mBranch2, x.mKaks) for x in self.mBranchInfo])
        self.mTreeSds = TreeTools.Graph2Tree(
            [(x.mBranch1, x.mBranch2, x.mSds) for x in self.mBranchInfo])
        self.mTreeNdn = TreeTools.Graph2Tree(
            [(x.mBranch1, x.mBranch2, x.mNdn) for x in self.mBranchInfo])
        self.mTreeS = TreeTools.Graph2Tree(
            [(x.mBranch1, x.mBranch2, x.mS) for x in self.mBranchInfo])
        self.mTreeN = TreeTools.Graph2Tree(
            [(x.mBranch1, x.mBranch2, x.mN) for x in self.mBranchInfo])


class BaseMLResult(CodeMLResult):

    """result object for BaseML."""

    def __init__(self):
        CodeMLResult.__init__(self)
        self.mAlpha = "na"
        self.mNumGammaBins = 0
        self.mKappa = "na"


class CodeMLResultSites(CodeMLResult):

    """result with site specific information."""

    def __init__(self, num_sequences, model):
        CodeMLResult.__init__(self)
        self.mNumSequences = num_sequences
        self.mModel = model
        self.mNEB = CodeMLResultPositiveSites()
        self.mBEB = CodeMLResultPositiveSites()


class CodeMLResultPositiveSites(CodeMLResult):

    def __init__(self):
        CodeMLResult.__init__(self)
        self.mPositiveSites = []


class CodeMLResultPositiveSite:

    def __init__(self, residue, aa, p, w, stderr):

        self.mResidue = residue
        self.mAA = aa
        self.mProbability = p
        self.mOmega = w
        self.mStdErr = stderr

    def __str__(self):
        return "\t".join(("%i" % self.mResidue,
                          self.mAA,
                          "%5.2f" % self.mProbability,
                          "%5.2f" % self.mOmega,
                          "%5.2f" % self.mStdErr))


class CodeMLResultPairs(CodeMLResult):

    """results for a pairwise codeml run."""

    def __init__(self):
        CodeMLResult.__init__(self)
        self.mPairs = []

    def fromResult(self, result):
        """build pairwise results from tree."""

        self.mPairs = []

        nodes = result.mTreeKs.get_terminals()
        ids = [result.mTreeKs.node(x).get_data().taxon for x in nodes]
        for n1 in range(0, len(nodes)):
            for n2 in range(0, n1):

                name1, name2 = nodes[n1], nodes[n2]
                pair = CodeMLResultPair()
                pair.mName1, pair.mName2 = ids[n1], ids[n2]
                pair.mLogLikelihood = result.mLogLikelihood

                pair.mTau, pair.mS, pair.mN, pair.mKa, pair.mKs = \
                    (0,
                     result.mTreeSds.distance(name1, name2),
                     result.mTreeNdn.distance(name1, name2),
                     result.mTreeKa.distance(name1, name2),
                     result.mTreeKs.distance(name1, name2))

                pair.mKaks = pair.mKa / pair.mKs
                pair.mKappa = result.mKappa

                pair.mRn, pair.mRs, pair.mRn0, pair.mRs0, pair.mBranchLength = None, None, None, None, None
                self.mPairs.append(pair)


class CodeMLResultPair(CodeMLResult):

    """results for a pairwise comparison."""

    def __init__(self):
        CodeMLResult.__init__(self)
        self.mError = ""


class CodeMLAncestralSequence:

    """an ancestral sequence."""

    def __init__(self, sequence, accuracy_per_site, accuracy_per_sequence):
        self.mSequence = sequence
        self.mAccuracyPerSite = accuracy_per_site
        self.mAccuracyPerSequence = accuracy_per_sequence

    def __str__(self):
        return "\t".join(("%f" % self.mAccuracyPerSite,
                          "%f" % self.mAccuracyPerSequence,
                          self.mSequence))


class CodeML:

    """
    """

    def __init__(self):

        self.mFilenameControl = "codeml.ctl"

        self.mWarnings = []

        self.mDefaults = {"noisy": "9",
                          "verbose": "1",
                          "runmode": "0",
                          "seqtype": "2",
                          "CodonFreq": "2",
                          "clock": "0",
                          "model": "0",
                          "NSsites": "0",
                          "icode": "0",
                          "fix_kappa": "0",
                          "kappa": "2",
                          "fix_omega": "0",
                          "omega": "0.4",
                          "fix_alpha": "1",
                          "alpha": "0",
                          "fix_rho": "1",
                          "rho": "0",
                          "RateAncestor": "0",
                          "Small_Diff": ".5e-6",
                          }

        self.mOptions = {
            "noisy": {"0": "how much rubbish on screen",
                      "1": "how much rubbish on screen",
                      "2": "how much rubbish on screen",
                      "3": "how much rubbish on screen",
                      "9":  "how much rubbish on screen"
                      },
            "verbose": {"0": "concise output",
                        "1": "detailed output",
                        "2": "too much output"
                        },
            "runmode": {"0": "user tree",
                        "1": "semi automatic",
                        "2": "automatic",
                        "3": "stepwise addition",
                        "4": "PerturbationNNI",
                        "5": "PerturbationNNI",
                        "-2": "pairwise"
                        },
            "seqtype": {"1": "Sequence type is codons",
                        "2": "Sequence type is AAs",
                        "3": "Sequence type is codons translated to AAs",
                        },
            "CodonFreq": {"0": "Uniform codon frequencies: 1/61 each",
                          "1": "Codon frequencies given by F1X4",
                          "2": "Codon frequencies given by F3X4",
                          "3": "Codon frequencies given by codon table",
                          },
            "ndata": {"#": "number of data sets to analyse",
                      },
            "clock": {"0": "no clock",
                      "1": "clock",
                      "2": "local clock",
                      "3": "Combined Analysis",
                      },
            "aaDist": {"0": "equal",
                       "+": "geometric",
                       "-": "linear",
                       "1": "G1974",
                       "2": "Miyata",
                       "3": "c",
                       "4": "p",
                       "5": "v",
                       "6": "a",
                       },
            "aaRatefile": {"jones.dat": "Jones",
                           "dayhoff.data": "Dayhoff",
                           "wag.dat": "Wag",
                           "mtmam.dat": "mtmam",
                           },
            "model": {"0": "Codon: one omega per branch - AA: poisson ",
                      "1": "Codon: branch specific omega - AA: proportional",
                      "2": "Codon: 2 or more omega - AA: Empirical",
                      "3": "Empirical+F",
                      "6": "FromCodon",
                      "7": "AAClasses",
                      "8": "REVaa_0",
                      "9": "REVaa(nr=189)",
                      },
            "NSsites": {"0": "no site-specific variation of omega",
                        "1": "neutral",
                        "2": "selection",
                        "3": "discrete",
                        "4": "freqs",
                        "5": "gamma",
                        "6": "2gamma",
                        "7": "beta",
                        "8": "beta&w",
                        "9": "beta&gamma",
                        "10": "beta&gamma+1",
                        "11": "beta&normal>1",
                        "12": "0&2normal>1",
                        "13": "3normal>0",
                        },
            "icode": {"0": "Genetic code: universal",
                      "1": "mammalian mt",
                      "2": "yeast mt.",
                      "3": "mold mt.",
                      "4": "invertebrate mt.",
                      "5": "ciliate nuclear",
                      "6": "echinoderm mt.",
                      "7": "euplotid mt.",
                      "8": "alternative yeast nu.",
                      "9": "ascidian mt.",
                      "10": "blepharisma nu."},
            "Mgene": {"0": "All rates equal between genes (if G option in sequence file: different, but proportional branch lengths)",
                      "1": "separate",
                      "2": "different pi",
                      "3": "different kappa",
                      "4": "all different",
                      },
            "fix_kappa": {"0": "kappa is to be estimated",
                          "1": "kappa is fixed",
                          },
            "kappa": {"2": "default initial or fixed kappa value.",
                      "#": "fixed or initial kappa",
                      },
            "fix_omega": {"0": "omega is to be estimated",
                          "1": "omega is fixed",
                          },
            "omega": {"0.4": "default initial of fixed omega value",
                      "#": "initial or fixed omega",
                      },
            "fix_alpha": {"0": "alpha is to be estimated",
                          "1": "alpha is fixed",
                          },
            "alpha": {"0": "inifinity: fixed alpha",
                      "#": "discrete gamma model with shape parameter alpha"
                      },
            "fix_rho": {"0": "rho is to be estimated",
                        "1": "rho is fixed",
                        },

            "fix_blength": {"0": "ignore branch lengths on given tree",
                            "-1": "use random branch lengths for given tree",
                            "1": "use branch lengths on given tree as initial values",
                            "2": "assume branch lengths on given tree to be fixed"},

            "rho": {"0": "no correlation between sites",
                    "#": "initial or fixed rho."
                    },
            "Malpha": {"0": "same alpha for genes",
                       "1": "different alpha for genes",
                       },
            "ncatG": {"#": "number of categories for discrete gamma model in site-specific models",
                      },
            "getSE":  {"0": "do not output standard errors of estimates",
                       "1": "output standard errors of estimates",
                       },
            "RateAncestor": {"0": "do not calculate ancestral states",
                             "1": "calculate ancestral states (only works with runmode 0)",
                             "2": "?"},
            "Small_Diff": {".5e-6": "default value for small difference",
                           "#": "small difference",
                           },
            "cleandata": {"0": "do not remove sites with ambiguity characters (default)",
                          "1": "remove sites with ambiguity characters",
                          },
            "method": {"0": "simultaneous calculation of parameters",
                       "1": "one branch at a time calculation of parameters",
                       },
        }

        self.mExecutable = "codeml"

        self.mFilenameTree = None
        self.mFilenameSequences = None

    def GetOptions(self):
        """return options in pretty format"""
        result = ["# Options for CodeML"]
        for var in list(self.mOptions.keys()):
            result.append("# %-40s: %s" % (var, self.mOptions[var]))
        return string.join(result, "\n")

    # ------------------------------------------------------------------------
    def SetOption(self, option, value):
        self.mDefaults[option] = value

    # ------------------------------------------------------------------------
    def AddOptions(self, parser):
        """add options to an OptionParser object."""
        parser.add_option("--paml-set-alpha", dest="paml_alpha", type="float",
                          help="initial alpha value [default=%default].")

        parser.add_option("--paml-fix-alpha", dest="paml_fix_alpha", action="store_true",
                          help="do not estimate alpha [default=%default].")

        parser.add_option("--paml-set-kappa", dest="paml_kappa", type="float",
                          help="initial kappa value [default=%default].")

        parser.add_option("--paml-fix-kappa", dest="paml_fix_kappa", action="store_true",
                          help="do not estimate kappa [default=%default].")

        parser.add_option("--paml-set-clean-data", dest="paml_clean_data", type="choice",
                          choices=("0", "1"),
                          help="PAML should cleanup data:  0=only gaps within pair are removed, 1=columns in the mali with gaps are removed [default=%default].")

        parser.set_defaults(
            paml_kappa=None,
            paml_fix_kappa=False,
            paml_alpha=None,
            paml_fix_alpha=False,
            paml_clean_data=None)

    # ------------------------------------------------------------------------
    def SetOptions(self, options):
        """set options from the command line."""

        if options.paml_kappa is not None:
            self.SetOption("kappa", str(options.paml_kappa))

        if options.paml_fix_kappa:
            self.SetOption("fix_kappa", "1")

        if options.paml_alpha is not None:
            self.SetOption("alpha", str(options.paml_alpha))

        if options.paml_fix_alpha:
            self.SetOption("fix_alpha", "1")

        if options.paml_clean_data:
            self.SetOption("cleandata", options.paml_clean_data)

    # ------------------------------------------------------------------------
    def WriteAlignment(self, mali):
        """write alignment in Phylip format."""

        outfile = open(self.mTempdir + "/" + self.mFilenameSequences, "w")

        mali.writeToFile(outfile, format="phylip")

        outfile.close()

    # ------------------------------------------------------------------------
    def WriteTree(self, tree):
        """write tree to file. The root of the tree is removed.
        """

        TreeTools.Unroot(tree)

        outfile = open(self.mTempdir + "/" + self.mFilenameTree, "w")
        outfile.write("%s\n" % TreeTools.Tree2Newick(tree))
        outfile.close()

    # ------------------------------------------------------------------------
    def writeControlFile(self,
                         outfile,
                         filename_sequences="input",
                         filename_output="output",
                         filename_tree=None,
                         options={},
                         ):
        """write a codeml.ctl file into outfile.
        """

        outfile.write("seqfile = %s\n" % (filename_sequences))

        if filename_tree:
            outfile.write("treefile = %s\n" % (filename_tree))

        outfile.write("outfile = %s\n\n" % (filename_output))

        written_options = {}

        # write set options
        for key, vv in list(options.items()):

            if key not in self.mOptions:
                raise IndexError("illegal option %s" % key)

            values = vv.strip().split(" ")

            descriptions = []
            for value in values:
                if value in self.mOptions[key]:
                    descriptions.append(self.mOptions[key][value])
                elif "#" in self.mOptions[key]:
                    try:
                        v = float(value)
                    except ValueError:
                        raise ValueError("not a number: %s for option %s" %
                                         (value, key))
                    descriptions.append(self.mOptions[key]["#"])
                else:
                    raise ValueError(
                        "illegal value %s for option %s: possible values are: %s" %
                        (value, key, " ".join(list(self.mOptions[key].keys()))))

            written_options[key] = 1
            description = ";".join(descriptions)
            outfile.write("* %s\n   %s = %s\n\n" % (description, key, vv))

        # write default options
        for key, vv in list(self.mDefaults.items()):

            if key in written_options:
                continue

            if key not in self.mOptions:
                raise IndexError("illegal option %s" % key)

            values = vv.strip().split(" ")

            descriptions = []
            for value in values:
                if value in self.mOptions[key]:
                    descriptions.append(self.mOptions[key][value])
                elif "#" in self.mOptions[key]:
                    try:
                        v = float(value)
                    except ValueError:
                        raise ValueError(
                            "not a number: %s for option %s" %
                            (value, key))
                    descriptions.append(self.mOptions[key]["#"])
                else:
                    raise ValueError(
                        "illegal value %s for option %s: possible values are: %s" %
                        (value, key, " ".join(list(self.mOptions[key].keys()))))

            written_options[key] = 1
            description = ";".join(descriptions)
            outfile.write("* %s\n   %s = %s\n\n" % (description, key, vv))

    # ------------------------------------------------------------------------
    def Run(self, alignment,
            tree=None,
            dump=0,
            test=False,
            options={}):

        self.mTempdir = tempfile.mkdtemp()
        self.mFilenameSequences = "input"
        self.mWarnings = []

        self.mFilenameOutput = "output"

        if test:
            print("# temporary directory is %s" % self.mTempdir)

        if tree:

            # check what kind of tree is given.
            if type(tree) == StringType:
                t = tree.strip()
                if t[0] == "(" and t[-1] in ");":
                    nexus = TreeTools.Newick2Nexus(tree)
                else:
                    nexus = TreeTools.Newick2Nexus(open(tree, "r"))
                t = nexus.trees[0]
            else:
                t = tree

            # check if identifiers in tree and alignment are the same
            identifiers_mali = set(alignment.getIdentifiers())
            identifiers_tree = set(TreeTools.GetTaxa(t))

            not_in_tree = identifiers_mali.difference(identifiers_tree)

            for x in not_in_tree:
                self.mWarnings.append(
                    "sequence %s not in tree - removed from multiple alignment." % x)
                alignment.deleteEntry(x)

            self.mFilenameTree = "tree"
            self.WriteTree(t)
        else:
            self.mFilenameTree = None

        self.mNumSequences = len(alignment)
        self.WriteAlignment(alignment)
        self.mRunOptions = options

        self.writeControlFile(open(self.mTempdir + "/" + self.mFilenameControl, "w"),
                              filename_tree=self.mFilenameTree,
                              options=options)

        s = subprocess.Popen("%s" % (self.mExecutable),
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             cwd=self.mTempdir,
                             close_fds=True)

        (out, err) = s.communicate()

        if s.returncode != 0:
            raise UsageError("Error in running %s \n%s\n%s\nTemporary directory in %s" %
                             (self.mExecutable, err, out, self.mTempdir))

        lines = open("%s/%s" %
                     (self.mTempdir, self.mFilenameOutput), "r").readlines()

        try:
            rst_lines = open("%s/rst", "r").readlines()
        except IOError:
            rst_lines = None

        if len(lines) == 0:
            raise UsageError("Empty result from %s \n%s\n%s\nTemporary directory in %s" %
                             (self.mExecutable, err, out, self.mTempdir))

        if dump:
            print(
                "############################### CONTROL FILE ############################")
            print("".join(open(self.mTempdir + "/" +
                               self.mFilenameControl, "r").readlines()))
            print("# result output of %s:\n%s\n######################################" % (
                self.mExecutable, "".join(lines)))
            print(
                "############################### LOG OUTPUT ############################")
            print("# stdout output of %s:\n%s\n######################################" % (
                self.mExecutable, out))

        if not test:
            shutil.rmtree(self.mTempdir)

        return self.parseOutput(lines, out.split("\n"), rst_lines)

    def parseRst(self, inlines, result):
        """parse lines from rst file."""

        lines = self.getSection(
            inlines, "List of extant and reconstructed sequences")

        # for PAML 4 use :4 instead of :2
        del lines[:4]

        # for aborted runs or ancestral reconstruction not performed: return
        # empty
        if len(lines) == 0:
            return

        result.mAncestralSequences = {}
        ancestral_sequences = []

        while lines[0] != "\n":
            id, sequence = lines[0][:18], lines[0][18:-1]
            id = id.strip()
            sequence = re.sub("\s", "", sequence)

            del lines[0]
            # do not save terminal sequences
            if id not in result.mSequences:
                if id[:6] != "node #":
                    raise ParsingError(
                        "expected ancestral node names to start with 'node #', instead got %s" % id)
                else:
                    id = id[6:]
                    ancestral_sequences.append((id, sequence))
            else:
                result.mAncestralSequences[
                    id] = CodeMLAncestralSequence(sequence, 1.0, 1.0)

        while lines and lines[0] == "\n":
            del lines[0]

        if not re.match("Overall accuracy", lines[0]):
            raise ParsingError(
                "expected 'Overall accuracy', instead got %s" % lines[0])

        del lines[0]
        accuracy_per_site = list(
            map(float, re.split("\s+", lines[0][:-1].strip())))
        del lines[:3]
        accuracy_per_sequence = list(map(
            float, re.split("\s+", lines[0][:-1].strip())))

        for x in range(len(ancestral_sequences)):
            id, sequence = ancestral_sequences[x]
            result.mAncestralSequences[id] = CodeMLAncestralSequence(
                sequence,
                accuracy_per_site[x],
                accuracy_per_sequence[x])

        result.mAncestralTree = TreeTools.Graph2Tree(
            [(x.mBranch1, x.mBranch2, 1.0) for x in result.mBranchInfo],
            label_ancestral_nodes=True)

    def parseSequences(self, lines, result):

        # read sequences
        # start from the top and read until line start not with
        # "seed used"
        while len(lines) and (lines[0] == "" or lines[0] == " " or lines[0].startswith("seed used")):
            del lines[0]

        self.mNumSequences, nchar = re.match(
            "(\d+)\s+(\d+)", lines[0]).groups()
        del lines[0]
        del lines[0]
        result.mSequences = {}
        self.mMapIdentifier2Position = {}
        self.mMapPosition2Identifier = []
        while lines[0]:
            try:
                identifier, sequence = re.match(
                    "(\S+)\s+(.+)", lines[0]).groups()
            except AttributeError:
                raise ParsingError("parsing error", lines[0])

            result.mSequences[identifier] = sequence
            self.mMapIdentifier2Position[identifier] = len(
                self.mMapIdentifier2Position)
            self.mMapPosition2Identifier.append(identifier)
            del lines[0]

    def parseVersion(self, lines, result):

        lines = self.getSection(lines, "CODONML", "BASEML")

        # read version information
        if lines[0].startswith("CODONML"):
            x = re.search(
                "CODONML \(in (.+)\)\s+\S+\s+Model:\s+(\S.+)", lines[0].strip())
            if x:
                result.mVersion, result.mModel = x.groups()
                del lines[0]

            else:
                # paml version 4.4c
                x = re.search("CODONML \(in (.+)\)\s+input", lines[0].strip())
                if x:
                    result.mVersion = x.groups()[0]
                    del lines[0]
                    try:
                        result.mModel = re.search(
                            "Model: (.*)$", lines[0].strip()).groups()[0]
                    except AttributeError:
                        raise ParsingError("pattern error", lines[0])
                else:
                    raise ParsingError("pattern error", lines[0])
                del lines[0]

            # read codon information
            try:
                result.mFrequencies = re.match(
                    "Codon frequenc.*: (\S+)", lines[0]).groups()[0]
                del lines[0]
            except AttributeError:
                pass

            # read sequence info
            try:
                result.mNumSequences, result.mLength = list(map(
                    int, re.match("ns\s*=\s+(\d+)\s+ls\s*=\s*(\d+)", lines[0]).groups()))
                del lines[0]
            except AttributeError:
                raise ParsingError(
                    "parsing error: expected ns = 123 ls = 123", lines[0])

        elif lines[0].startswith("BASEML"):
            try:
                result.mVersion, result.mModel = re.search(
                    "BASEML \(in (.+)\)\s+\S+\s+(\S.+)", lines[0].strip()).groups()
                del lines[0]
            except AttributeError:
                raise ParsingError("pattern error", lines[0])
        else:
            raise ParsingError("unknown PAML program", lines[0])

        known_versions = ("paml 3.14b, May 2005",
                          "paml 3.15, November 2005",
                          "paml version 4, June 2007",
                          "paml version 4.4c, August 2010",
                          "paml version 4.8a, August 2014")

        if result.mVersion not in known_versions:
            raise ValueError("unknown paml version '%s'" % result.mVersion)

    def parseSitePatterns(self, inlines, result):
        # read site patterns
        # todo: reading of site patterns

        lines = self.getSection(inlines, "Printing out site pattern counts")
        del lines[0]
        self.nextSection(lines)

        del lines[0]
        # result.mNumSitePatterns = int( re.match("# site patterns = (\d+)",
        # lines[0]).group(1) )
        result.mNumSitePatterns = 0
        del lines[0]

    def parseSequenceDifferences(self, lines, result):
        pass

    def parseCodonUsage(self, inlines, result):

        lines = self.getSection(inlines, "Codon usage in sequences")

    def checkSection(self, lines, section_start):
        """check if section starts with string section_start."""

        if not lines:
            raise ParsingError("section error - ran out of lines")

        if not lines[0].startswith(section_start):
            raise ParsingError("section error - expected '%s'" % section_start,
                               "\n".join(lines[:3]))

    def getSection(self, lines, *args):
        """check if section starts with string section_start."""

        for x in range(len(lines)):
            for a in args:
                if lines[x].startswith(a):
                    return lines[x:]
        else:
            raise ParsingError(
                "section error - can't find section starting with '%s'" %
                str(args))

    def parseNeiGojobori(self, inlines, result):

        lines = self.getSection(inlines,
                                "Nei & Gojobori 1986. dN/dS (dN, dS)")

    def parseResultsKaks(self, inlines, result):

        lines = self.getSection(inlines, "TREE")
        del lines[0]

        if re.match("check convergence..", lines[0]):
            del lines[0]
            result.mCheckConvergence = True
        else:
            result.mCheckConvergence = False

        if re.match("This is a rooted tree.", lines[0]):
            result.mCheckRootedTree = False
            del lines[0]
            self.nextSection(lines)
        else:
            result.mCheckRootedTree = True

        try:
            a, b, c = re.match(
                "lnL\(ntime:\s*(\d+)\s*np:\s*(\d+)\):\s*(\S+)\s+\S+", lines[0]).groups()
        except AttributeError:
            raise ParsingError("parsing error", lines[0])

        result.mNTime, result.mNumParameters, result.mLogLikelihood = int(
            a), int(b), float(c)
        del lines[0]

        #################################################################
        while not lines[0] == "Detailed output identifying parameters":
            del lines[0]
        del lines[0]
        while not lines[0]:
            del lines[0]

        #################################################################
        if re.match("kappa", lines[0]):
            result.mKappa = float(
                re.match("kappa \(ts/tv\)\s*=\s*(\S+)", lines[0]).groups()[0])
            del lines[0]
            self.nextSection(lines)

        if re.match("Parameters in beta", lines[0]):
            del lines[0]
            result.mBeta = {}
            while lines[0]:

                data = re.split("[\s=]+", re.sub("[()]", "", lines[0]))
                del lines[0]
                for k, val in [(data[x], data[x + 1]) for x in range(0, len(data), 2)]:
                    result.mBeta[k] = float(val)

            self.nextSection(lines)
        #################################################################
        if re.match("dN/dS for site classes", lines[0]):
            result.mSiteClassesNumSites = re.search(
                "dN/dS for site classes \(K=(\d+)\)", lines[0]).groups()[0]
            self.skipSection(lines)

            # note: spaces might be missing, because it seems to a fixed with format.
            # This is ok for probabilities (allways less than 1) but a problem for
            # w, if it gets larger than 100. Field width seems to be 9.
            self.mSiteClassesProbabilities = list(map(
                float, re.split("\s+", lines[0])[1:]))
            del lines[0]

            ncats = len(self.mSiteClassesProbabilities)
            s = lines[0][3:]
            self.mSiteClassesOmega = [float(s[x:x + 9])
                                      for x in range(0, len(s), 9)]

            del lines[0]
            self.nextSection(lines)

        #################################################################
        if lines[0][:len("omega (dN/dS)")] == "omega (dN/dS)":
            result.mOmega = float(
                re.match("omega \(dN/dS\)\s+=\s+(\S+)", lines[0]).groups()[0])
            self.skipSection(lines)

        if lines[0] != "dN & dS for each branch":
            raise ParsingError("section error", lines[0])
        del lines[0]
        while not lines[0]:
            del lines[0]
        del lines[0]
        del lines[0]

        result.mBranchInfo = []

        while lines[0]:

            branch, t, S, N, dnds, dn, ds, sds, ndn = re.split("\s+", lines[0])

            # start counting from 0, correct below to for ancestral nodes
            # to use PAML 1-based numbering.
            a, b = [int(x) - 1 for x in branch.split("..")]

            if a < result.mNumSequences:
                branch1 = self.mMapPosition2Identifier[a]
            else:
                branch1 = str(a + 1)

            if b < result.mNumSequences:
                branch2 = self.mMapPosition2Identifier[b]
            else:
                branch2 = str(b + 1)

            result.mBranchInfo.append(
                CodeMLBranchInfo(branch1, branch2,
                                 float(dnds), float(dn), float(ds),
                                 float(ndn), float(sds),
                                 int(float(S)), int(float(N))))
            del lines[0]

        result.buildTreesFromBranchInfo()

        while not lines[0]:
            del lines[0]

        #################################################################
        # the following is optional (only when site specific model = 0)
        if re.match("tree length for dN", lines[0]):
            result.mTreeLengthKa = float(
                re.match("tree length for dN:\s+(\S+)", lines[0]).groups()[0])
            del lines[0]
            result.mTreeLengthKs = float(
                re.match("tree length for dS:\s+(\S+)", lines[0]).groups()[0])
            del lines[0]

        #################################################################
        if result.mModel in ("free dN/dS Ratios for branches",):
            while not lines[0]:
                del lines[0]
            del lines[0]
            result.mPamlTreeKs = lines[0]
            del lines[0]
            del lines[0]
            result.mPamlTreeKa = lines[0]
            del lines[0]

        self.nextSection(lines)

    def nextSection(self, lines):
        while len(lines) and not lines[0]:
            del lines[0]

    def skipSection(self, lines):

        while len(lines) and lines[0]:
            del lines[0]
        while len(lines) and not lines[0]:
            del lines[0]

    def saveSummary(self, result, lines, lines_log=None):

        # save summary
        result.mResult = "".join(lines)

        if lines_log:
            if lines_log[0] and lines_log[0][-1] == "\n":
                lines_log = [x[:-1].strip() for x in lines_log]
            result.mLog = "\n".join(lines_log)
            self.parseLog(lines_log, result)

    def parseLog(self, lines_log, result):
        """parse log output."""

        rx = re.compile(
            "\s+rN\s*=\s*(\S+)\s+rS\s*=\s*(\S+)\s+rN\*\s*=\s*(\S+)\s+rS\*\s*=\s*(\S+)\s+blength\s*=\s*(\S+)")

        result.mRn, result.mRs, result.mRn0, result.mRs0, result.mBranchLength = None, None, None, None, None

        for line in lines_log:

            x = rx.match(line)
            if x:
                result.mRn, result.mRs, result.mRn0, result.mRs0, result.mBranchLength = list(map(
                    float, x.groups()))

    def parseOutput(self, lines, lines_log=None, rst_lines=None):
        """parse CodeML output. This is rather tricky, as paml output is as
        freeformat as it can get.  Also, there is a log file and an
        output file. Proceed sequentially through file.

        """

        result = CodeMLResult()

        self.saveSummary(result, lines, lines_log)

        # chop, strip and remove comments
        lines = [x[:-1].strip() for x in lines]

        if len(lines) == 0:
            raise ParsingError("empty input")

        # parse by sections
        self.parseVersion(lines, result)

        self.parseSequences(lines, result)

        self.parseSitePatterns(lines, result)

        self.parseSequenceDifferences(lines, result)

        self.parseCodonUsage(lines, result)

        self.parseResultsKaks(lines, result)

        if rst_lines:
            self.parseRst(rst_lines, result)

        return result


class CodeMLSites (CodeML):

    def __init__(self):
        CodeML.__init__(self)

    def parseOutput(self, lines, lines_log=None, rst_lines=None):
        """parse codeml output for site-specific analysis."""

        result = CodeMLResult()
        self.saveSummary(result, lines, lines_log)

        # chop and strip
        lines = [x[:-1].strip() for x in [x for x in lines if x[0] != "#"]]

        self.parseSequences(lines, result)

        self.parseSitePatterns(lines, result)

        self.parseVersion(lines, result)

        self.parseSequenceDifferences(lines, result)

        self.parseCodonUsage(lines, result)

        self.parseNeiGojobori(lines, result)

        result.mSites = {}

        # todo: fix parser
        while 1:
            if not lines:
                break

            if not re.match("Model", lines[0]):
                raise ParsingError("section error", lines[0])

            id, name = re.match("Model (\d+): (.+)", lines[0]).groups()
            self.skipSection(lines)

            if id not in result.mSites:
                result.mSites[id] = CodeMLResultSites(
                    num_sequences=result.mNumSequences,
                    model=result.mModel)

            r = result.mSites[id]
            r.mSiteModelId, r.mSiteModelName = id, name

            # parse model info
            self.parseResultsKaks(lines, r)

            self.parseSites(lines, r)

        if rst_lines:
            self.parseRst(rst_lines, result)

        return result

    def parseSitesBayes(self, lines, result):

        if re.match("Naive Empirical Bayes \(NEB\) analysis", lines[0]):
            naive = True
        elif re.match("Bayes Empirical Bayes \(BEB\) analysis", lines[0]):
            naive = False
        else:
            raise ParsingError("section error", lines[0])

        del lines[0]

        if re.match("Time used:", lines[0]):
            return

        if re.match("Bayes Empirical Bayes \(BEB\) analysis", lines[0]):
            return

        del lines[:4]
        while 1:
            if not lines[0]:
                break

            if naive:
                residue, aa, p, w = re.match(
                    "(\d+)\s+(\S)\s+([0-9.]+)[*]{0,2}\s+([0-9.]+)", lines[0]).groups()
                stderr = 0
            else:
                residue, aa, p, w, stderr = re.match(
                    "(\d+)\s+(\S)\s+([0-9.]+)[*]{0,2}\s+([0-9.]+) \+\- ([0-9.]+)", lines[0]).groups()
            del lines[0]
            result.mPositiveSites.append(
                CodeMLResultPositiveSite(int(residue), aa, float(p),
                                         float(w), float(stderr)))

        self.nextSection(lines)

    def parseGrids(self, lines, result):
        """parse grid information."""

        if re.match("The grid", lines[0]):
            del lines[:2]

            result.mGrid = {}
            while lines[0]:
                k, vals = lines[0].split(":")
                data = re.split("\s+", vals.strip())
                result.mGrid[k.strip()] = list(map(float, data))
                del lines[0]

            self.nextSection(lines)

        if re.match("Posterior on the grid", lines[0]):
            del lines[:2]

            result.mGridPosterior = {}
            while lines[0]:
                k, vals = lines[0].split(":")
                data = re.split("\s+", vals.strip())
                result.mGridPosterior[k.strip()] = list(map(float, data))
                del lines[0]

            self.nextSection(lines)

        if re.match("Posterior for p0-p1", lines[0]):
            del lines[:2]

            result.mPosteriorP0P1 = []
            while lines[0]:
                result.mPosteriorP0P1.append(lines[0])
                del lines[0]

            del lines[0]

            result.PosteriorP0P1Sum = re.match(
                "sum of density on p0-p1\s+=\s+([0-9.]+)", lines[0]).groups()[0]
            del lines[0]

            self.nextSection(lines)

    def parseSites(self, lines, result):
        """parse site specific model results."""

        while 1:

            if not lines:
                break

            if re.match("Time used:", lines[0]):
                del lines[0]
                break

            if re.match("Naive Empirical Bayes \(NEB\) analysis", lines[0]):
                self.parseSitesBayes(lines, result.mNEB)

            elif re.match("Bayes Empirical Bayes \(BEB\) analysis", lines[0]):
                self.parseSitesBayes(lines, result.mBEB)

            elif re.match("The grid", lines[0]):
                self.parseGrids(lines, result)

            elif re.match("Note:", lines[0]):
                result.mWarnings.append(lines[0])
                self.skipSection(lines)

            else:
                raise ParsingError("section error", lines[0])

        self.nextSection(lines)


class CodeMLPairwise (CodeML):

    def __init__(self):
        CodeML.__init__(self)

        self.mDefaults["runmode"] = "-2"

    def parseLog(self, lines_log, result):
        """parse log output.

        This routine collects the rho values for each pair.
        """
        self.mRhos = []

        rx = re.compile(
            "\s+rN\s*=\s*(\S+)\s+rS\s*=\s*(\S+)\s+rN\*\s*=\s*(\S+)\s+rS\*\s*=\s*(\S+)\s+blength\s*=\s*(\S+)")

        for line in lines_log:

            x = rx.match(line)
            if x:
                self.mRhos.append(list(map(float, x.groups())))

    def parseOutput(self, lines, lines_log=None, rst_lines=None):
        """parse codeml output for pairwise rate calculation."""

        result = CodeMLResultPairs()

        self.saveSummary(result, lines, lines_log)

        # chop and strip
        lines = [x[:-1].strip() for x in [x for x in lines if x[0] != "#"]]

        self.parseVersion(lines, result)

        # skip over the next parts, we are not interested
        while lines and lines[0][:19] != "pairwise comparison":
            del lines[0]

        del lines[0]

        self.nextSection(lines)

        self.parsePairs(lines, result)

        return result

    def parsePairs(self, lines, result):
        """parse pairwise results."""

        while lines:
            pair = CodeMLResultPair()

            n1, name1, n2, name2 = re.match(
                "(\d+)\s*\((\S+)\)\s*\.\.\.\s*(\d+)\s*\((\S+)\)", lines[0]).groups()
            del lines[0]

            pair.mName1, pair.mName2 = name1, name2

            pair.mLogLikelihood = float(
                re.match("lnL =\s*(\S+)", lines[0]).groups()[0])
            del lines[0]

            # the next values are tau, kappa and omega (this might
            # change, if some of these are fixed?)
            values = list(map(float, re.split("\s+", lines[0].strip())))
            del lines[0]

            if len(values) == 3:
                pair.mTau, pair.mKappa, pair.mOmega = values
            elif "fix_omega" in self.mRunOptions and self.mRunOptions["fix_omega"] == "1":
                pair.mTau, pair.mKappa = values
                pair.mOmega = "na"
            elif "fix_kappa" in self.mRunOptions and self.mRunOptions["fix_kappa"] == "1":
                pair.mTau, pair.mOmega = values
                pair.mKappa = "na"
            else:
                raise ParsingError("can't parse omega and kappa", lines[0])

            del lines[0]

            tokens = re.split("[\s=]+", lines[0].strip())
            pair.mTau, pair.mS, pair.mN, pair.mKaks, pair.mKa, pair.mKs = list(map(
                float, [tokens[x] for x in range(1, len(tokens), 2)]))
            del lines[0]

            if self.mRhos:
                pair.mRn, pair.mRs, pair.mRn0, pair.mRs0, pair.mBranchLength = self.mRhos[
                    len(result.mPairs)]
            else:
                pair.mRn, pair.mRs, pair.mRn0, pair.mRs0, pair.mBranchLength = None, None, None, None, None

            result.mPairs.append(pair)

            self.nextSection(lines)


class BaseML(CodeML):

    def __init__(self):

        CodeML.__init__(self)

        self.mFilenameControl = "baseml.ctl"
        self.mDefaults = {
            "noisy": "9",
            "verbose": "0",
            "runmode": "0",
            "model": "4",
            "Mgene": "0",
            "clock": "0",
            "fix_kappa": "0",
            "kappa": "5",
            "fix_alpha": "1",
            "alpha": "0",
            "Malpha": "0",
            "ncatG": "5",
            "nparK": "0",
            "nhomo": "0",
            "getSE": "0",
            "RateAncestor": "0",
            "Small_Diff": "7e-6",
            "cleandata": "1",
            "method": "0",
        }

        self.mOptions = {
            "nhomo": {"0": "homogeneous",
                      "1": "homogeneous",
                      "2": "kappa for branches",
                      "3": "N1",
                      "3": "N2",
                      },

            "nparK": {"0": "no rate classes",
                      "1": "rate class model - rK",
                      "2": "rate class model - rK & fK",
                      "3": "rate class model - rK & MK(1/K)",
                      "4": "rate class model - rK & MK",
                      },

            "noisy": {"0": "how much rubbish on screen",
                      "1": "how much rubbish on screen",
                      "2": "how much rubbish on screen",
                      "3": "how much rubbish on screen",
                      "9":  "how much rubbish on screen"
                      },
            "verbose": {"0": "concise output",
                        "1": "detailed output",
                        "2": "too much output"
                        },
            "runmode": {"0": "user tree",
                        "1": "semi automatic",
                        "2": "automatic",
                        "3": "stepwise addition",
                        "4": "PerturbationNNI",
                        "5": "PerturbationNNI",
                        "-2": "pairwise"
                        },
            "seqtype": {"1": "Sequence type is codons",
                        "2": "Sequence type is AAs",
                        "3": "Sequence type is codons translated to AAs",
                        },
            "CodonFreq": {"0": "Uniform codon frequencies: 1/61 each",
                          "1": "Codon frequencies given by F1X4",
                          "2": "Codon frequencies given by F3X4",
                          "3": "Codon frequencies given by codon table",
                          },
            "ndata": {"#": "number of data sets to analyse",
                      },
            "clock": {"0": "no clock",
                      "1": "clock",
                      "2": "local clock",
                      "3": "Combined Analysis",
                      },
            "aaDist": {"0": "equal",
                       "+": "geometric",
                       "-": "linear",
                       "1": "G1974",
                       "2": "Miyata",
                       "3": "c",
                       "4": "p",
                       "5": "v",
                       "6": "a",
                       },
            "aaRatefile": {"jones.dat": "Jones",
                           "dayhoff.data": "Dayhoff",
                           "wag.dat": "Wag",
                           "mtmam.dat": "mtmam",
                           },
            "model": {"0": "JC69",
                      "1": "K80",
                      "2": "F81",
                      "3": "F84",
                      "4": "HKY85",
                      "5": "T92",
                      "6": "TN93",
                      "7": "REV",
                      "8": "UNREST",
                      "9": "REVu",
                      "10": "UNRESTu",
                      },

            "NSsites": {"0": "no site-specific variation of omega",
                        "1": "neutral",
                        "2": "selection",
                        "3": "discrete",
                        "4": "freqs",
                        "5": "gamma",
                        "6": "2gamma",
                        "7": "beta",
                        "8": "beta&w",
                        "9": "beta&gamma",
                        "10": "beta&gamma+1",
                        "11": "beta&normal>1",
                        "12": "0&2normal>1",
                        "13": "3normal>0",
                        },
            "icode": {"0": "Genetic code: universal",
                      "1": "mammalian mt",
                      "2": "yeast mt.",
                      "3": "mold mt.",
                      "4": "invertebrate mt.",
                      "5": "ciliate nuclear",
                      "6": "echinoderm mt.",
                      "7": "euplotid mt.",
                      "8": "alternative yeast nu.",
                      "9": "ascidian mt.",
                      "10": "blepharisma nu."},
            "Mgene": {"0": "All rates equal between genes (if G option in sequence file: different, but proportional branch lengths)",
                      "1": "separate",
                      "2": "different pi",
                      "3": "different kappa",
                      "4": "all different",
                      },

            "fix_kappa": {"0": "kappa is to be estimated",
                          "1": "kappa is fixed",
                          },
            "kappa": {"2": "default initial or fixed kappa value.",
                      "#": "fixed or initial kappa",
                      },
            "fix_omega": {"0": "omega is to be estimated",
                          "1": "omega is fixed",
                          },
            "omega": {"0.4": "default initial of fixed omega value",
                      "#": "initial or fixed omega",
                      },
            "fix_alpha": {"0": "alpha is to be estimated",
                          "1": "alpha is fixed",
                          },
            "alpha": {"0": "inifinity: fixed alpha",
                      "#": "discrete gamma model with shape parameter alpha"
                      },
            "fix_rho": {"0": "rho is to be estimated",
                        "1": "rho is fixed",
                        },

            "fix_blength": {"0": "ignore branch lengths on given tree",
                            "-1": "use random branch lengths for given tree",
                            "1": "use branch lengths on given tree as initial values",
                            "2": "assume branch lengths on given tree to be fixed"},

            "rho": {"0": "no correlation between sites",
                    "#": "initial or fixed rho."
                    },
            "Malpha": {"0": "same alpha for genes",
                       "1": "different alpha for genes",
                       },
            "ncatG": {"#": "number of categories for discrete gamma model in site-specific models",
                      },
            "getSE":  {"0": "do not output standard errors of estimates",
                       "1": "output standard errors of estimates",
                       },
            "RateAncestor": {"0": "do not calculate ancestral states",
                             "1": "calculate ancestral states (only works with runmode 0)",
                             "2": "?"},
            "Small_Diff": {".5e-6": "default value for small difference",
                           "#": "small difference",
                           },
            "cleandata": {"0": "do not remove sites with ambiguity characters (default)",
                          "1": "remove sites with ambiguity characters",
                          },
            "method": {"0": "simultaneous calculation of parameters",
                       "1": "one branch at a time calculation of parameters",
                       },
        }

        self.mExecutable = "baseml"

    # ------------------------------------------------------------------------
    def AddOptions(self, parser):
        """add options to an OptionParser object."""
        CodeML.AddOptions(self, parser)

        parser.add_option("--baseml-model", dest="baseml_model", type="choice",
                          choices=("T92", "JC69", "F84",
                                   "K80", "F81", "HKY85", "TN93", "REV", "UNREST", "REVU",
                                   "UNRESTU"),
                          help="baseml model [default=%default].")

        parser.set_defaults(baseml_model="REV")

    # ------------------------------------------------------------------------
    def SetOptions(self, options):
        """set options from the command line."""
        CodeML.SetOptions(self, options)

        map_model2index = {}
        for key, val in list(self.mOptions["model"].items()):
            map_model2index[val] = key

        if options.baseml_model.upper() in map_model2index:
            self.SetOption("model", map_model2index[options.baseml_model])
        else:
            raise "unknown model for baseml: %s" % options.baseml_model

    def parseOutput(self, lines, lines_log=None, rst_lines=None):
        """parse BASEML output. This is rather tricky, as paml output is as
        freeformat as it can get.  Also, there is a log file and an
        output file. Proceed sequentially through file.

        """

        result = BaseMLResult()

        self.saveSummary(result, lines, lines_log)

        # chop and strip
        lines = [x[:-1].strip() for x in lines]
        # lines = map(lambda x: x[:-1].strip(), filter( lambda x: x[0] != "#",
        # lines ))

        result = BaseMLResult()

        self.parseVersion(lines, result)

        # self.parseSitePatterns( lines, result )

        # self.parseSequenceDifferences( lines, result )

        self.parseFrequencies(lines, result)

        self.parseTree(lines, result)

        self.parseParameters(lines, result)

        return result

    def parseFrequencies(self, inlines, result):
        """parse frequency section."""

        lines = self.getSection(inlines, "Frequencies..")

        del lines[0]
        del lines[0]
        result.mBaseFrequencies = {}
        for x in range(self.mNumSequences):
            data = re.split("\s+", lines[0])
            del lines[0]
            f = {}
            for n, y in enumerate(("T", "C", "A", "G")):
                f[y] = data[n + 1]

            result.mBaseFrequencies[data[0]] = f

        self.nextSection(lines)

        if not re.match("Homogeneity statistic", lines[0]):
            raise ParsingError(
                "wrong section, expected 'Homogeneity statistic'", lines[0])

        result.mX2, result.mG = list(map(float, re.search(
            "Homogeneity statistic: X2 = (\S+) G = (\S+)", lines[0]).groups()))
        del lines[0]
        self.nextSection(lines)

        data = re.split("\s+", lines[0])
        del lines[0]
        f = {}
        for n, y in enumerate(("T", "C", "A", "G")):
            f[y] = data[n + 1]

        result.mBaseFrequencies["average"] = f

        self.nextSection(lines)

        self.nextSection(lines)

    def parseTree(self, lines, result):

        lines = self.getSection(lines, "TREE")
        del lines[0]

        try:
            a, b, c = re.match(
                "lnL\(ntime:\s*(\d+)\s*np:\s*(\d+)\):\s*(\S+)\s+\S+", lines[0]).groups()
        except AttributeError:
            raise ParsingError("parsing error", lines[0])

        result.mNTime, result.mNumParameters, result.mLogLikelihood = int(
            a), int(b), float(c)
        del lines[0]

        while lines and not re.match("tree length", lines[0]):
            del lines[0]

        result.mTreeLength = re.match(
            "tree length\s+=\s+(\S+)", lines[0]).groups()[0]
        del lines[0]

        self.nextSection(lines)
        # skip over tree topology
        self.skipSection(lines)

        # read rate tree
        result.mTree = lines[0]
        del lines[0]

        # convert tree to matrix
        tree = TreeTools.Newick2Tree(result.mTree)

        result.mDistanceMatrix = {}

        nodes = tree.get_terminals()
        for x in range(len(nodes) - 1):

            key1 = tree.node(nodes[x]).data.taxon

            if key1 not in result.mDistanceMatrix:
                result.mDistanceMatrix[key1] = {}

            for y in range(x + 1, len(nodes)):

                key2 = tree.node(nodes[y]).data.taxon
                if key2 not in result.mDistanceMatrix:
                    result.mDistanceMatrix[key2] = {}

                d = tree.distance(nodes[x], nodes[y])

                result.mDistanceMatrix[key1][key2] = d
                result.mDistanceMatrix[key2][key1] = d

        self.nextSection(lines)

    def parseParameters(self, inlines, result):

        try:
            lines = self.getSection(
                inlines, "Detailed output identifying parameters")
        except ParsingError:
            return

        # read parameters
        del lines[0]
        self.nextSection(lines)

        if re.match("Parameters \(kappa\)", lines[0]):
            del lines[0]
            result.mKappa = float(lines[0])
            del lines[0]
            self.nextSection(lines)
        else:
            result.mKappa = "na"

        if not lines:
            return

        if re.match("Parameters  in the rate matrix", lines[0]):
            while not lines[0].startswith("Rate parameters"):
                del lines[0]
            d = re.split(
                "\s+", re.search("Rate parameters: +(.*)", lines[0]).groups()[0])
            result.mRevRateBaseParameters = list(map(float, d))
            del lines[0]
            result.mRevBaseFrequencies = list(map(
                float, re.split("\s+", re.search("Base frequencies: +(.*)", lines[0]).groups()[0])))
            del lines[0]
            try:
                result.mRevK = float(
                    re.match("Rate matrix Q, Average Ts/Tv =\s*(\S+)", lines[0]).groups()[0])
            except AttributeError:
                raise ParsingError("could not parse average Ts/Tv", lines[0])
            result.mKappa = float(result.mRevK)

            del lines[0]
            self.mRevQ = []
            for x in range(0, 4):
                self.mRevQ.append(list(map(float, re.split("\s+", lines[0]))))
                del lines[0]
            self.nextSection(lines)

        if not len(lines):
            return

        if re.match("check convergence..", lines[0]):
            del lines[0]
            result.mCheckConvergence = True
            if len(lines) == 0:
                return

        try:
            result.mNumGammaBins, result.mAlpha = re.match(
                "alpha \(gamma, K=(\d+)\)\s+=\s+(\S+)", lines[0]).groups()
        except AttributeError:
            raise ParsingError("could not parse alpha and gamma", lines[0])

        result.mAlpha, result.mNumGammaBins = float(
            result.mAlpha), int(result.mNumGammaBins)
        del lines[0]
        r = re.split("\s+", lines[0])[1:]
        del lines[0]
        f = re.split("\s+", lines[0])[1:]
        result.mGammaBins = []
        for x in range(result.mNumGammaBins):
            result.mGammaBins.append({"r": float(r[x]), "f": float(f[x])})


class Evolver:

    """interface class for running evolver.
    """

    def __init__(self):

        self.mExecutable = "evolver"
        random.seed()
        x = random.randint(0, 10000001)
        while x % 2 == 0:
            x = random.randint(0, 10000001)

        self.mSeed = x
        self.mReplicates = 1
        self.mNucleotides = 1000
        self.mNSequences = 2
        self.mTree = "(1:1,2:1);"
        self.mScale = -1
        self.mOmega = 0.3
        self.mKappa = 2
        self.mCodonTable = None
        self.mDs = 0.5

        # runmode 6 = codons (5=na, 7=aa)
        self.mRunMode = 6

    # ------------------------------------------------------------------------
    def setOmega(self, omega):
        self.mOmega = omega

    def setKappa(self, kappa):
        self.mKappa = kappa

    def setLength(self, length):
        self.mNucleotides = length

    def setDs(self, ds):
        self.mDs = ds

    def setReplicates(self, replicates):
        self.mReplicates = replicates

    # ------------------------------------------------------------------------
    def writeControlFile(self, outfile):
        """write control file to outfile."""

        self.calculateScale(self.mDs)

        if self.mCodonTable is None:
            raise ValueError("please supply a codon table")

        self.mNSequences = len(
            TreeTools.Newick2Nexus(self.mTree).trees[0].get_taxa())

        outfile.write("0          * paml format output\n")
        outfile.write("%i         * random number seed\n" % self.mSeed)
        outfile.write("%i %i %i   * nsequences nnucleotides nreplicates\n\n" %
                      (self.mNSequences,
                       self.mNucleotides,
                       self.mReplicates))

        outfile.write(
            "%f          * tree length (-1 = unscaled)\n" % self.mScale)
        outfile.write(self.mTree + "\n\n")
        outfile.write("%f          * omega\n" % self.mOmega)
        outfile.write("%f          * kappa\n" % self.mKappa)
        outfile.write("\n")

        for c1 in ("T", "C", "A", "G"):
            for c2 in ("T", "C", "A", "G"):
                for c3 in ("T", "C", "A", "G"):
                    codon = c1 + c2 + c3
                    outfile.write(" %10.9f" % self.mCodonTable[codon])
                outfile.write("\n")

        outfile.write("// end of file\n")

    # ------------------------------------------------------------------------
    def fromMali(self, mali):
        """compute codon table from a multiple alignment."""
        self.mCodonTable = {}
        for c1 in ("T", "C", "A", "G"):
            for c2 in ("T", "C", "A", "G"):
                for c3 in ("T", "C", "A", "G"):
                    codon = c1 + c2 + c3
                    self.mCodonTable[codon] = 0

        for id, i in list(mali.items()):
            s = i.mString
            for x in range(0, len(s), 3):
                codon = s[x:x + 3].upper()
                try:
                    self.mCodonTable[codon] += 1
                except KeyError:
                    continue

        total = sum(self.mCodonTable.values())
        for codon in list(self.mCodonTable.keys()):
            self.mCodonTable[codon] = float(self.mCodonTable[codon]) / total

    # ------------------------------------------------------------------------
    def setUniformFrequencies(self):
        """use uniform codon frequencies."""
        self.mCodonTable = {}

        frequency = 1.0 / 61.0

        for c1 in ("T", "C", "A", "G"):
            for c2 in ("T", "C", "A", "G"):
                for c3 in ("T", "C", "A", "G"):
                    codon = c1 + c2 + c3
                    if Genomics.IsStopCodon(codon):
                        self.mCodonTable[codon] = 0.0
                    else:
                        self.mCodonTable[codon] = frequency

    # ------------------------------------------------------------------------
    def calculateScale(self, ds):
        """calculate tree scale for a given dS value.

        The branch scale is given by:

        t = 3 dS * ps + 3 omega * dS * (1-ps)
        t = 3 dS * (ps + omega (1 - ps )
        """

        if self.mCodonTable is None:
            raise ValueError("please supply a codon table")

        # number of synonymous/non-synonymous sites
        Q, t = RateEstimation.getQMatrix(
            self.mCodonTable, self.mKappa, 1.0, self.mKappa, 1.0)

        rI, rV, rS, rN = RateEstimation.countSubstitutions(self.mCodonTable, Q)

        ps = rS / (rS + rN)

        self.mScale = 3.0 * self.mDs * (ps + self.mOmega * (1 - ps))

    # ------------------------------------------------------------------------
    def setTree(self, tree):
        """set tree."""
        self.mTree = tree

    # ------------------------------------------------------------------------
    def run(self, ds=None, tree=None, test=False, dump=False):
        """run evolver."""

        self.mTempdir = tempfile.mkdtemp()

        self.mWarnings = []

        self.mFilenameOutput = "mc.paml"

        if test:
            print("# temporary directory is %s" % self.mTempdir)

        if tree:
            # check what kind of tree is given.
            if type(tree) == StringType:
                t = tree.strip()
                if t[0] == "(" and t[-1] in ");":
                    self.mTree = t
                else:
                    nexus = TreeTools.Newick2Nexus(open(tree, "r"))
                    self.mTree = nexus.trees[0].to_string()

        if ds:
            self.mDs = ds

        self.mFilenameControl = "evolver.ctl"
        self.writeControlFile(
            open(self.mTempdir + "/" + self.mFilenameControl, "w"))

        if dump:
            print("################### control file input  ########")
            infile = open(self.mTempdir + "/" + self.mFilenameControl, "r")
            print("".join(infile.readlines()))
            infile.close()

        s = subprocess.Popen("%s %i %s" %
                             (self.mExecutable, self.mRunMode,
                              self.mFilenameControl),
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             cwd=self.mTempdir,
                             close_fds=True)

        (out, err) = s.communicate()

        if s.returncode != 0:
            raise UsageError(
                "Error in running %s \n%s\n%s\nTemporary directory in %s" %
                (self.mExecutable, err, out, self.mTempdir))

        lines = open("%s/%s" %
                     (self.mTempdir, self.mFilenameOutput), "r").readlines()

        if len(lines) == 0:
            raise UsageError(
                "Empty result from %s \n%s\n%s\nTemporary directory in %s" %
                (self.mExecutable, err, out, self.mTempdir))

        if dump:
            print("# result output of %s:\n%s\n######################################" % (
                self.mExecutable, "".join(lines)))
            print(
                "############################### LOG OUTPUT ############################")
            print("# stdout output of %s:\n%s\n######################################" % (
                self.mExecutable, out))

        if not test:
            shutil.rmtree(self.mTempdir)

        return self.parseOutput(lines, out.split("\n"))

    def parseOutput(self, lines, lines_log=None, rst_lines=None):

        result = []
        chunk = []
        for line in lines:
            line = line.strip()
            if not line:
                continue
            if re.match("\d+", line):
                if chunk:
                    mali = Mali.Mali()
                    mali.readFromFile(chunk, format="phylip")
                    result.append(mali)
                chunk = []
            chunk.append(line + "\n")

        if chunk:
            mali = Mali.Mali()
            mali.readFromFile(chunk, format="phylip")
            result.append(mali)

        return result

    def printFrequencies(self, outfile):

        if not self.mCodonTable:
            outfile.write("# no frequency table defined\n")
        else:
            codons = list(self.mCodonTable.keys())
            outfile.write("# codon table used in evolver:\n")
            for codon in codons:
                outfile.write("# %s\t%6.4f\n" %
                              (codon, self.mCodonTable[codon]))


class EvolverBaseml(Evolver):

    """interface class for running evolver for nucleotides.
    """

    def __init__(self, *args, **kwargs):

        Evolver.__init__(self, *args, **kwargs)

        self.mFrequencies = None

        # runmode 5 = codons (5=na, 7=aa)
        self.mRunMode = 5

        self.mKappa = 2.0
        self.mModel = "K80"
        self.mModels = {"JC69": 0,
                        "K80": 1,
                        "F81": 2,
                        "F84": 3,
                        "HKY85": 4,
                        "T92": 5,
                        "TN93": 6,
                        "REV": 7}

        # 0 is same rate over all sites
        self.mAlpha = 0

        # 0 is continuous gamma
        self.mGammaCategories = 0

    # ------------------------------------------------------------------------
    def setModel(self, model):
        self.mModel = model.upper()

    # ------------------------------------------------------------------------
    def setUniformFrequencies(self):
        self.mFrequencies = {}
        for c in ("T", "C", "A", "G"):
            self.mFrequencies[c] = 0.25

    # ------------------------------------------------------------------------
    def fromMali(self, mali):
        """compute frequencies from a multiple alignment."""
        self.mFrequencies = {}
        for c in ("T", "C", "A", "G"):
            self.mFrequencies[c] = 0

        for id, i in list(mali.items()):
            for c in i.mString:
                try:
                    self.mFrequencies[c] += 1
                except KeyError:
                    continue

        total = sum(self.mFrequencies.values())
        for c in list(self.mFrequecies.keys()):
            self.mFrequencies[c] = float(self.mFrequencies[c]) / total

    # ------------------------------------------------------------------------
    def getParameters(self):
        """get parameters for a model.

        From the MCbase.dat:
        Parameter kappa or rate parameters in the substituton model:
        For TN93, two kappa values are required, while for REV, 5 values
        (a,b,c,d,e) are required (see Yang 1994 for the definition of these
        parameters).
        The kappa parameter is defined differently under HKY85 (when k=1 means
        no transition bias) and under F84 (when k=0 means no bias).
        JC69 and F81 are considered species cases of HKY85, so use 1 for kappa
        for those two models.  Notation is from my two papers in JME in 1994.
        """

        if self.mModel in ("JC69", "F81"):
            self.mKappa = 1.0

        if self.mModel == "REV":
            return (10, 5, 1, 2, 3)

        elif self.mModel == "TN93":
            raise NotImplementedError("not implemented")
            return (k1, k2)

        else:
            return (self.mKappa,)

    # ------------------------------------------------------------------------
    def writeControlFile(self, outfile):
        """write control file to outfile."""

        if self.mModel in ("JC69", "K80"):
            self.setUniformFrequencies()

        if self.mFrequencies is None:
            raise ValueError("please supply nucleotide frequencies")

        self.mNSequences = len(
            TreeTools.Newick2Nexus(self.mTree).trees[0].get_taxa())

        outfile.write("0          * paml format output\n")
        outfile.write("%i         * random number seed\n" % self.mSeed)
        outfile.write("%i %i %i   * nsequences nnucleotides nreplicates\n\n" %
                      (self.mNSequences,
                       self.mNucleotides,
                       self.mReplicates))

        outfile.write("%f          * tree length (-1 = unscaled)\n" % self.mDs)
        outfile.write(self.mTree + "\n\n")
        outfile.write("%i          * model: 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85, 5:T92, 6:TN93, 7:REV\n" %
                      self.mModels[self.mModel])
        outfile.write("%s          * kappa or rate parameters in model\n" %
                      " ".join(map(str, self.getParameters())))
        outfile.write("%f %i       * <alpha>  <#categories for discrete gamma>\n" %
                      (self.mAlpha, self.mGammaCategories))
        outfile.write("\n")

        for c in ("T", "C", "A", "G"):
            outfile.write(" %10.9f" % self.mFrequencies[c])
        outfile.write("\n")
        for c in ("T", "C", "A", "G"):
            outfile.write(" %s" % c)
        outfile.write("\n")

        outfile.write("// end of file\n")

    def printFrequencies(self, outfile):

        if not self.mFrequencies:
            outfile.write("# no frequency table defined\n")
        else:
            outfile.write("# frequency table used in evolver:\n")
            for c in list(evolver.mFrequencies.keys()):
                outfile.write("# %s\t%6.4f\n" % (c, evolver.mFrequencies[c]))


def getOptions(options):
    """translate command line options to PAML options."""

    codeml_options = {}

    if options.analysis == "branch-specific-kaks":
        codeml_options["seqtype"] = "1"
        codeml_options["model"] = "1"

    elif options.analysis == "branch-fixed-kaks":
        codeml_options["seqtype"] = "1"
        codeml_options["model"] = "0"

    elif options.analysis == "branch-all-but-one-fixed-kaks":
        codeml_options["seqtype"] = "1"
        codeml_options["model"] = "2"
        if not tree:
            raise ValueError("please supply a tree for this mode.")
        if not options.filename_output_tree:
            raise ValueError(
                "please speficy filename-output-tree as location "
                "(relative to this script) for trees.")

    elif options.analysis == "site-specific-kaks":
        codeml_options["ncatG"] = "10"
        codeml_options["getSE"] = "1"
        codeml_options["seqtype"] = "1"
        codeml_options["NSsites"] = "0 3 1 2 7 8"
        codeml_options["model"] = "0"
        codeml_options["CodonFreq"] = "2"

    elif options.analysis == "pairwise":
        codeml_options["seqtype"] = "1"
        codeml_options["model"] = "0"
        codeml_options["runmode"] = "-2"

    if options.multiple_genes:
        codeml_options["Malpha"] = "0"
        codeml_options["Mgene"] = "0"

    if options.omega is not None:
        codeml_options["omega"] = str(options.omega)

    if options.estimate_ancestors:
        codeml_options["RateAncestor"] = "1"

    if options.codon_frequencies is not None:
        c = options.codon_frequencies.upper()
        if c in ("UNIFORM", "FEQUAL"):
            a = "0"
        elif c == "F1X4":
            a = "1"
        elif c == "F3X4":
            a = "2"
        elif c == "F61":
            a = "3"
        else:
            a = options.codon_frequencies
        codeml_options["CodonFreq"] = a
    if options.method is not None:
        codeml_options["method"] = str(options.method)

    if options.optimization_threshold is not None:
        codeml_options["Small_Diff"] = str(options.optimization_threshold)

    if options.clean_data:
        codeml_options["cleandata"] = options.clean_data

    return codeml_options


def runEvolver(options):
    """run evolver."""

    if options.evolver_model == "codon":
        evolver = Evolver()
        if options.omega:
            evolver.setOmega(options.omega)
    else:
        evolver = EvolverBaseml()
        evolver.setModel(options.evolver_model)

    if options.kappa:
        evolver.setKappa(options.kappa)

    if options.evolver_tree:
        evolver.setTree(options.evolver_tree)

    if options.evolver_length:
        evolver.setLength(options.evolver_length)

    if options.evolver_ds:
        evolver.setDs(options.evolver_ds)

    if options.evolver_replicates:
        evolver.setReplicates(options.evolver_replicates)

    if options.filename_sequences != "input":
        mali = Mali.Mali()
        if options.loglevel >= 1:
            options.stdlog.write(
                "# reading multiple alignment from %s\n" % (options.filename_sequences))
            options.stdlog.flush()
        mali.readFromFile(open(options.filename_sequences, "r"),
                          format="fasta")
        if options.loglevel >= 1:
            options.stdlog.write("# calculating frequency table\n")
            options.stdlog.flush()

        evolver.fromMali(mali)

        if options.loglevel >= 2:
            evolver.printFrequencies(options.stdlog)
    else:
        if options.loglevel >= 1:
            options.stdlog.write("# using uniform codon frequencies.\n")

        evolver.setUniformFrequencies()

    if options.write_control_file:
        evolver.writeControlFile(sys.stdout)
        return

    result = evolver.run(
        test=options.test,
        dump=options.dump)

    for r in result:
        if options.evolver_pairwise:
            ids = r.getIdentifiers()
            for x in range(0, len(ids)):
                for y in range(0, x):
                    options.stdout.write(">%s\n%s\n>%s\n%s\n" %
                                         (ids[x],
                                          r.getSequence(ids[x]),
                                          ids[y],
                                          r.getSequence(ids[y])))
        else:
            r.writeToFile(options.stdout, format="fasta")

    if options.loglevel >= 2:
        options.stdlog.write("# Parameters used for evolver\n")
        options.stdlog.write("# dS=%f t=%f kappa=%f omega=%f\n" %
                             (evolver.mDs,
                              evolver.mScale,
                              evolver.mKappa,
                              evolver.mOmega))
        options.stdlog.write("# tree: %s\n" % evolver.mTree)

if __name__ == "__main__":

    parser = E.OptionParser(
        version="%prog version: $Id$")

    parser.add_option(
        "--write-control-file", dest="write_control_file", action="store_true",
        help="write a control file.")

    parser.add_option(
        "--flavour", dest="flavour", type="choice",
        choices=("codeml", "baseml", "evolver"),
        help="codeml or baseml or evolver, that is the question.")

    parser.add_option(
        "--analysis", dest="analysis", type="choice",
        choices=("branch-specific-kaks", "branch-fixed-kaks",
                 "branch-all-but-one-fixed-kaks", "site-specific-kaks",
                 "pairwise"),
        help="choose an analysis scenario.")

    parser.add_option(
        "--multiple-genes", dest="multiple_genes", action="store_true",
        help="analyse multiple genes.")
    parser.add_option(
        "-t", "--tree-nh-file", dest="filename_tree", type="string",
        help="filename with tree information. Evolver: tree and number of sequences.")
    parser.add_option(
        "-i", "--filename-sequences", dest="filename_sequences", type="string",
        help="filename with sequences. Evolver: determines codon frequencies.")
    parser.add_option(
        "-o", "--filename-output", dest="filename_output", type="string",
        help="filename for output information.")
    parser.add_option(
        "--filename-clusters", dest="filename_clusters", type="string",
        help="filename for cluster information.")
    parser.add_option(
        "--filename-output-tree", dest="filename_output_tree", type="string",
        help="filename pattern for trees to be output.")
    parser.add_option(
        "-l", "--filename-log", dest="filename_log", type="string",
        help="filename for logging information.")
    parser.add_option(
        "--filename-pattern-control", dest="filename_pattern_control",
        type="string",
        help="filename with pattern for control files to create.")
    parser.add_option(
        "--output-pattern-id", dest="output_pattern_id", type="string",
        help="output pattern for id.")
    parser.add_option(
        "-p", "--parse-output", dest="parse_output", type="choice",
        choices=("none", "all", "terminal-nodes", "likelihood",
                 "ancestral-sequence",
                 "ks-tree", "ka-tree", "kaks-tree", "sds-tree",
                 "ndn-tree", "s-tree", "n-tree",
                 "mali", "sequences"),
        help="parse output")
    parser.add_option(
        "-b", "--parse-batch", dest="parse_batch", type="string",
        help="""supply a batch file for output parsing. The batch "
        "file contains the following (tab-separated): "
        "1. name, 2. output-filename, 3. logfile-filename, "
        "4. rst-filename, 5. mali filename(fasta) (all options "
        "after the second are optional and can be empty.)""")

    parser.add_option(
        "--parse-tabular", dest="parse_tabular", action="store_true",
        help="""output parsing results in tabular format (default = fasta-like format).""")
    parser.add_option(
        "--filename-rst", dest="filename_rst", type="string",
        help="filename with ancestral data. Needed for parsing ancestral sequences.")
    parser.add_option("--set-omega", dest="omega", type="float",
                      help="initial omega value. Evolver: omega values used for simulation.")
    parser.add_option("--set-kappa", dest="kappa", type="float",
                      help="initial kappa value. Evolver: kappa values used for simulation.")
    parser.add_option("--set-codon-frequencies", dest="codon_frequencies", type="string",
                      help="set codon frequencies.")
    parser.add_option("--set-method", dest="method", type="int",
                      help="set paml optimization method [0|1].")
    parser.add_option("--set-optimization-threshold", dest="optimization_threshold", type="string",
                      help="set paml optimization threshold.")
    parser.add_option("--estimate-ancestors", dest="estimate_ancestors", action="store_true",
                      help="estimat ancestral sequences.")
    parser.add_option("--filename-read-tree", dest="input_filename_tree", type="string",
                      help="filename with tree information to read (synonym for input-filename-tree).")
    parser.add_option("--input-filename-tree", dest="input_filename_tree", type="string",
                      help="filename with tree information.")
    parser.add_option("--test", dest="test", action="store_true",
                      help="run a test.")
    parser.add_option("--dump", dest="dump", action="store_true",
                      help="dump output.")
    parser.add_option("--set-clean-data", dest="clean_data", type="choice",
                      choices=("0", "1"),
                      help="PAML should cleanup data:  0=only gaps within pair are removed, 1=columns in the mali with gaps are removed.")

    parser.add_option("--filter-output",  dest="filter_output", type="choice", action="append",
                      choices=("ndn", "sds"),
                      help="filter output. Give values as options to --filter_parameters.")

    parser.add_option("--filter-parameters", dest="filter_parameters", type="string",
                      help="comma-separated list of filter values.")

    parser.add_option("--evolver-length", dest="evolver_length", type="int",
                      help="evolver: number of nucleotides/amino acids/codons to simulate.")

    parser.add_option("--evolver-replicates", dest="evolver_replicates", type="int",
                      help="evolver: number of replicates to simulate.")

    parser.add_option("--evolver-model", dest="evolver_model", type="choice",
                      choices=(
                          "T92", "JC69", "F84", "K80", "F81", "HKY85", "TN93", "REV", "codon"),
                      help="model for evolver[default=%default].")

    parser.add_option("--evolver-ds", dest="evolver_ds", type="float",
                      help="evolver: ds to use to scale tree.")

    parser.add_option("--evolver-tree", dest="evolver_tree", type="string",
                      help="Evolver: tree to use. Determines the number of sequences. Default are two sequences.")

    parser.add_option("--filename-mali", dest="filename_mali", type="string",
                      help="Filename with mali in fasta format. Used for special purpose: map species names back to genes for 1:1 ortholog cluster.")

    parser.add_option("--map-tsv-file", dest="filename_map_old2new", type="string",
                      help="Filename with map of identifiers mapping from old to new identifiers.")

    parser.add_option("--invert-map", dest="invert_map", action="store_true",
                      help="invert the mapping, so that old2new is read from new2old [%default].")

    parser.add_option("--evolver-pairwise", dest="evolver_pairwise", action="store_true",
                      help="output results from evolver simulation as all-on-all sequence pairs [%default].")

    parser.set_defaults(
        write_control_file=False,
        analysis=None,
        filename_tree=None,
        filename_output="output",
        filename_output_tree=None,
        filename_sequences="input",
        filename_pattern_control="-",
        filename_clusters=None,
        filename_log=None,
        filename_rst=None,
        filename_mali=None,
        filename_map_old2new=None,
        multiple_genes=False,
        parse_output=None,
        omega=None,
        kappa=None,
        codon_frequencies=None,
        method=None,
        optimization_threshold=None,
        output_pattern_id="%s",
        input_filename_tree=None,
        test=False,
        dump=False,
        flavour="codeml",
        clean_data=None,
        filter_output=None,
        filter_parameters=None,
        evolver_tree=None,
        evolver_length=None,
        evolver_ds=None,
        evolver_model="codon",
        estimate_ancestors=False,
        parse_batch=None,
        parse_tabular=False,
        separator="|",
        invert_map=False,
        evolver_pairwise=False,
    )

    (options, args) = Experiment.Start(parser)

    if options.filter_parameters is not None:
        options.filter_parameters = options.filter_parameters.split(",")
        options.filter_parameters.reverse()

    if options.flavour == "codeml":
        if options.analysis in ("branch-specific-kaks", "branch-fixed-kaks",
                                "branch-all-but-one-fixed-kaks",
                                "pairwise"):
            codeml = CodeML()

        elif options.analysis == "site-specific-kaks":
            codeml = CodeMLSites()

        else:
            raise ValueError("unknown analysis scenario %s" % options.analysis)
    elif options.flavour == "baseml":
        if 1:
            codeml = BaseML()
        else:
            raise ValueError("unknown analysis scenario %s" % options.analysis)
    elif options.flavour == "evolver":

        runEvolver(options)
        Experiment.Stop()
        sys.exit(0)

    else:
        raise ValueError("unknown flavour %s" % options.flavour)

    if options.input_filename_tree:
        nexus = TreeTools.Newick2Nexus(open(options.input_filename_tree, "r"))
        Tree.updateNexus(nexus)
        tree = nexus.trees[0]
    else:
        tree = None

    ##########################################################################
    ##########################################################################
    ##########################################################################
    # build options
    ##########################################################################
    codeml_options = getOptions(options)

    ##########################################################################
    ##########################################################################
    ##########################################################################
    # write control file(s)
    ##########################################################################
    if options.write_control_file:

        if options.filename_clusters:
            clusters, nerrors = IOTools.ReadList(
                open(options.filename_clusters, "r"))
        else:
            clusters = ("", )

        nwritten = 0

        for cluster in clusters:

            id = options.output_pattern_id % cluster
            fn_sequences = re.sub("%s", id, options.filename_sequences)
            fn_output = re.sub("%s", id, options.filename_output)
            if options.filename_tree:
                fn_tree = re.sub("%s", id, options.filename_tree)
            else:
                fn_tree = None

            fn_ctl = re.sub("%s", id, options.filename_pattern_control)

            if options.analysis == "branch-all-but-one-fixed-kaks":

                branch = 0
                # setup branch specific models.
                # You need to create a new tree for each one
                for node_id, node in list(tree.chain.items()):
                    if node.prev is None:
                        continue

                    branch += 1

                    fn_branch_ctl = fn_ctl % branch
                    fn_branch_tree = fn_tree % branch
                    fn_branch_out = fn_output % branch

                    old_name = node.data.taxon

                    if node.succ == []:
                        node.data.taxon += " # 1 "
                    else:
                        node.data.taxon = " $ 1 "

                    outfile = open(options.filename_output_tree % branch, "w")
                    outfile.write(TreeTools.Tree2Newick(
                        tree, with_branch_lengths=False, write_all_taxa=True) + "\n")
                    outfile.close()

                    outfile = open(fn_branch_ctl, "w")

                    codeml.writeControlFile(outfile,
                                            filename_sequences=fn_sequences,
                                            filename_output=fn_branch_out,
                                            filename_tree=fn_branch_tree,
                                            options=codeml_options)
                    outfile.close()

                    node.data.taxon = old_name
                    nwritten += 1
            else:

                if fn_ctl == "-":
                    outfile = options.stdout
                else:
                    outfile = open(fn_ctl, "w")

                codeml.writeControlFile(outfile,
                                        filename_sequences=fn_sequences,
                                        filename_output=fn_output,
                                        filename_tree=fn_tree,
                                        options=codeml_options)

                if outfile != options.stdout:
                    outfile.close()

                nwritten += 1

        if options.loglevel >= 1:
            options.stdout.write("# written %i control files.\n" % nwritten)

    ##########################################################################
    ##########################################################################
    ##########################################################################
    # parse output
    ##########################################################################
    elif options.parse_output:

        filters = []
        if options.filter_output:
            for f in options.filter_output:
                if f == "sds":
                    p = options.filter_parameters.pop()
                    filters.append(lambda x: x.mSds < p)
                elif f == "ndn":
                    p = options.filter_parameters.pop()
                    filters.append(lambda x: x.mNdn < p)

        batch = []
        ninput, noutput, nerrors = 0, 0, 0

        if options.parse_batch:
            batch = []
            for line in open(options.parse_batch, "r"):
                if line[0] == "#":
                    continue
                data = line[:-1].split("\t")
                id, filename_log, filename_rst, filename_mali, filename_map_old2new = None, None, None, None, None
                if len(data) == 2:
                    filename_out, id = data
                elif len(data) == 3:
                    filename_out, id, filename_log = data
                elif len(data) == 4:
                    filename_out, id, filename_log, filename_rst = data
                elif len(data) == 5:
                    filename_out, id, filename_log, filename_rst, filename_mali = data
                elif len(data) == 6:
                    filename_out, id, filename_log, filename_rst, filename_mali, filename_map_old2new = data
                else:
                    raise ValueError(
                        "please supply at least two arguments (filename + id): %s" % str(data))

                batch.append(
                    (id, filename_out, filename_log, filename_rst, filename_mali, filename_map_old2new))

        else:
            batch.append((None, None, options.filename_log, options.filename_rst,
                          options.filename_mali, options.filename_map_old2new))

        first = True

        for cluster_id, filename_out, filename_log, filename_rst, filename_mali, filename_map_old2new in batch:

            if options.loglevel >= 2:
                options.stdlog.write("# processing: id=%s, out=%s, log=%s, rst=%s, mali=%s\n" %
                                     (cluster_id, filename_out, filename_log, filename_rst, filename_mali))

            try:
                if filename_log:
                    log_lines = open(filename_log, "r").readlines()
                else:
                    log_lines = None

                if filename_rst:
                    rst_lines = open(filename_rst, "r").readlines()
                else:
                    rst_lines = None

                if filename_out:
                    out_lines = open(filename_out, "r").readlines()
                else:
                    out_lines = sys.stdin.readlines()

            except IOError as msg:
                options.stdlog.write(
                    "file not found error for cluster id %s: %s\n" % (str(cluster_id), msg))
                traceback.print_exc()
                nerrors += 1
                continue

            try:
                result = codeml.parseOutput(out_lines, log_lines, rst_lines)
            except ParsingError as msg:
                options.stdlog.write(
                    "parsing error for cluster id %s: %s\n" % (str(cluster_id), msg))
                traceback.print_exc()
                nerrors += 1
                continue

            if cluster_id and first and options.parse_tabular:
                options.stdout.write("id\t")

            if filename_mali:
                mali = Mali.Mali()
                infile = open(filename_mali, "r")
                mali.readFromFile(infile)
                map_old2new = {}
                for id in mali.getIdentifiers():
                    species = id.split(options.separator)[0]
                    map_old2new[species] = id

                result.mapNames(map_old2new)
                infile.close()

            if filename_map_old2new:
                map_old2new = IOTools.readMap(open(filename_map_old2new))
                if options.invert_map:
                    map_old2new = IOTools.invert_dictionary(
                        map_old2new, make_unique=True)
                result.mapNames(map_old2new)

            if options.parse_output == "all":
                output = "%s" % str(result)

            elif options.parse_output == "likelihood":
                if first and options.parse_tabular:
                    options.stdout.write("lnL\tnp\n")
                output = "%f\t%i" % (
                    result.mLogLikelihood, result.mNumParameters)

            elif options.parse_output == "terminal-nodes":
                if first and options.parse_tabular:
                    options.stdout.write("branch\tka\tks\tkaks\n")
                for branch in result.mBranchInfo:
                    if branch.mBranch1 in result.mSequences:
                        output = "%s\t%s\t%s\t%s" % (branch.mBranch1,
                                                     branch.mKa, branch.mKs, branch.mKaks)
                    elif branch.mBranch2 in result.mSequences:
                        output = "%s\t%s\t%s\t%s" % (branch.mBranch2,
                                                     branch.mKa, branch.mKs, branch.mKaks)

            elif options.parse_output == "ks-tree":
                output = TreeTools.Tree2Newick(
                    result.mTreeKs, with_branch_lengths=True, write_all_taxa=True)

            elif options.parse_output == "ka-tree":
                output = TreeTools.Tree2Newick(
                    result.mTreeKa, with_branch_lengths=True, write_all_taxa=True)

            elif options.parse_output == "kaks-tree":
                output = TreeTools.Tree2Newick(
                    result.mTreeKaks, with_branch_lengths=True, write_all_taxa=True)

            elif options.parse_output == "sds-tree":
                output = TreeTools.Tree2Newick(
                    result.mTreeSds, with_branch_lengths=True, write_all_taxa=True)

            elif options.parse_output == "ndn-tree":
                output = TreeTools.Tree2Newick(
                    result.mTreeNdn, with_branch_lengths=True, write_all_taxa=True)

            elif options.parse_output == "s-tree":
                output = TreeTools.Tree2Newick(
                    result.mTreeS, with_branch_lengths=True, write_all_taxa=True)

            elif options.parse_output == "n-tree":
                output = TreeTools.Tree2Newick(
                    result.mTreeN, with_branch_lengths=True, write_all_taxa=True)

            elif options.parse_output == "sequences":
                output = []
                for id, s in list(result.mSequences.items()):
                    output.append(">%s\n%s" % (id, s))
                output = "\n".join(output)

            elif options.parse_output == "ancestral-sequence":

                if not hasattr(result, "mAncestralSequences"):
                    raise ValueError(
                        "no ancestral sequences defined in result")

                # return most ancestral sequence
                # mid-point root the tree with ks and then return the sequence closest to the root.
                # this is not the true ancestral sequence, as PAML works on
                # unrooted trees.
                result.mTreeKs.root_midpoint()
                succ = result.mTreeKs.node(result.mTreeKs.root).succ
                assert(len(succ) == 2)
                if result.mTreeKs.node(succ[0]).data.branchlength < result.mTreeKs.node(succ[1]).data.branchlength:
                    s = succ[0]
                else:
                    s = succ[1]

                t = result.mAncestralTree.node(s).data.taxon
                a = result.mAncestralSequences[t]

                if first and options.parse_tabular:
                    options.stdout.write(
                        "taxon\tacc_site\tacc_seq\tsequence\n")

                output = "%s\t%f\t%f\t%s" % (t,
                                             a.mAccuracyPerSite,
                                             a.mAccuracyPerSequence,
                                             a.mSequence)

            elif options.parse_output == "mali":

                # output multiple alignment information
                if first and options.parse_tabular:
                    options.stdout.write("nseq\tlength\tnsite_patterns\n")

                output = "%i\t%i\t%i" % (
                    result.mNumSequences, result.mLength, result.mNumSitePatterns)

            if options.parse_tabular:
                if cluster_id:
                    options.stdout.write("%s\t" % cluster_id)
            else:
                if cluster_id:
                    options.stdout.write(">%s\n" % (cluster_id))

            options.stdout.write("%s\n" % (output))

            noutput += 1
            first = False

        if options.loglevel >= 1:
            options.stdlog.write(
                "# ninput=%i, noutput=%i, nerrors=%i\n" % (ninput, noutput, nerrors))

    ##########################################################################
    ##########################################################################
    ##########################################################################
    # do a small test
    ##########################################################################
    elif options.test:

        alignment = (
            ('Hsa_Human',
             "AAGGTCTTTGAAAGGTGTGAGTTGGCCAGAACTCTGAAAAGATTGGGAATGGATGGCTACAGGGGAATCAGCCTAGCAAACTGGATGTGTTTGGCCAAATGGGAGAGTGGTTACAACACACGAGCTACAAACTACAATGCTGGAGACAGAAGCACTGATTATGGGATATTTCAGATCAATAGCCGCTACTGGTGTAATGATGGCAAAACCCCAGGAGCAGTTAATGCCTGTCATTTATCCTGCAGTGCTTTGCTGCAAGATAACATCGCTGATGCTGTAGCTTGTGCAAAGAGGGTTGTCCGTGATCCACAAGGCATTAGAGCATGGGTGGCATGGAGAAATCGTTGTCAAAACAGAGATGTCCGTCAGTATGTTCAAGGTTGTGGAGTG"),
            ("Hla_gibbon",
             "AAGGTCTTTGAAAGGTGTGAGTTGGCCAGAACTCTGAAAAGATTGGGAATGGATGGCTACAGGGGAATCAGCCTAGCAAACTGGATGTGTTTGGCCAAATGGGAGAGTGGTTATAACACACGAGCTACAAACTACAATCCTGGAGACAGAAGCACTGATTATGGGATATTTCAGATCAATAGCCGCTACTGGTGTAATGATGGCAAAACCCCAGGAGCAGTTAATGCCTGTCATTTATCCTGCAATGCTTTGCTGCAAGATAACATCGCCGATGCTGTAGCTTGTGCAAAGAGGGTTGTCCGCGATCCACAAGGCATTAGAGCATGGGTGGCATGGAGAAATCGTTGTCAAAACAGAGATCTCCGTCAGTATATTCAAGGTTGTGGAGTA"),
            ("Cgu/Can_colobus",
             "AAGATCTTTGAAAGGTGTGAGTTGGCCAGAACTCTGAAAAAATTGGGACTGGATGGCTACAAGGGAGTCAGCCTAGCAAACTGGGTGTGTTTGGCCAAATGGGAGAGTGGTTATAACACAGACGCTACAAACTACAATCCTGGAGATGAAAGCACTGATTATGGGATATTTCAGATCAATAGCCGCTACTGGTGTAATAATGGCAAAACCCCAGGAGCAGTTAATGCCTGTCATATATCCTGCAATGCTTTGCTGCAAAATAACATCGCTGATGCTGTAGCTTGTGCAAAGAGGGTTGTCAGTGATCCACAAGGCATTCGAGCATGGGTGGCATGGAAAAAGCACTGTCAAAACAGAGATGTCAGTCAGTATGTTGAAGGTTGTGGAGTA"),
            ("Pne_langur",
             "AAGATCTTTGAAAGGTGTGAGTTGGCCAGAACTCTGAAAAAATTGGGACTGGATGGCTACAAGGGAGTCAGCCTAGCAAACTGGGTGTGTTTGGCCAAATGGGAGAGTGGTTATAACACAGAAGCTACAAACTACAATCCTGGAGACGAAAGCACTGATTATGGGATATTTCAGATCAATAGCCGCTACTGGTGTAATAATGGCAAAACCCCAGGAGCAGTTGATGCCTGTCATATATCCTGCAGTGCTTTGCTGCAAAACAACATCGCTGATGCTGTAGCTTGTGCAAAGAGGGTTGTCAGTGATCCACAAGGCGTTCGAGCATGGGTGGCATGGAGAAATCACTGTCAAAACAAAGATGTCAGTCAGTACGTTAAAGGTTGTGGAGTG"),
            ("Mmu_rhesus",
             "AAGATCTTTGAAAGGTGTGAGTTGGCCAGAACTCTGAAAAGATTGGGACTGGATGGCTACAGGGGAATCAGCCTAGCAAACTGGGTGTGTTTGGCCAAATGGGAGAGTAATTATAACACACAAGCTACAAACTACAATCCTGGAGACCAAAGCACTGATTATGGGATATTTCAGATCAATAGCCACTACTGGTGTAATAATGGCAAAACCCCAGGAGCAGTTAATGCCTGTCATATATCCTGCAATGCTTTGCTGCAAGATAACATCGCTGATGCTGTAACTTGTGCAAAGAGGGTTGTCAGTGATCCACAAGGCATTAGAGCATGGGTGGCATGGAGAAATCACTGTCAAAACAGAGATGTCAGTCAGTATGTTCAAGGTTGTGGAGTG"),
            ("Ssc_squirrelM",
             "AAGGTCTTCGAAAGGTGTGAGTTGGCCAGAACTCTGAAAAGGCTTGGAATGGATGGCTACAGGGGAATCAGCCTAGCAAACTGGATGTGTTTGGCCAAATGGGAGAGTGACTATAACACACGTGCTACAAACTACAATCCTGGAGACCAAAGCACTGATTATGGGATATTTCAGATCAATAGCCACTATTGGTGTAATAATGGCAGAACCCCAGGAGCAGTTAATGCCTGTCATATATCCTGCAATGCTTTGCTGCAAGATGACATCACTCAAGCTGTGGCCTGTGCAAAGAGGGTTGTCCGTGATCCACAAGGCATTAGAGCATGGGTGGCATGGAAAGCTCATTGTCAAAACAGAGATGTCAGTCAGTATGTTCAAGGTTGTGGAGTA"),
            ("Cja_marmoset",
             "AAGGTCTTTGAAAGGTGTGAGTTGGCCAGAACTCTGAAAAGGTTTGGACTGGATGGCTACAGGGGAATCAGCCTAGCAAACTGGATGTGTTTGGCCAAATGGGAGAGTGATTATAACACACGTGCTACAAACTACAATCCTGGAGACCAAAGCACTGATTATGGGATATTTCAGATCAATAGCCACTATTGGTGTAACAATGGCAGAACCCCAGGAGCAGTTAATGCCTGTCATATATCCTGCAATGCTTTGCTGCAAGATGACATCACTGAAGCTGTGGCCTGTGCAAAGAGGGTTGTCCGCGATCCACAAGGCATTAGGGCATGGGTGGCATGGAAAGCTCATTGTCAAAACAGAGATGTCAGTCAGTATGTTCAAGGTTGTGGAGTA"),
        )

        mali = Mali.Mali()
        for key, s in alignment:
            mali.addSequence(key, 0, len(s), s)

        tree = "((Hsa_Human, Hla_gibbon),((Cgu/Can_colobus, Pne_langur),Mmu_rhesus), (Ssc_squirrelM, Cja_marmoset));"

        result = codeml.Run(mali, tree,
                            options=codeml_options,
                            dump=options.dump)

        print(result)
    else:

        mali = Mali.Mali()

        mali.readFromFile(
            open(options.filename_sequences, "r"), format="fasta")

        # run it on input and dump out the output
        result = codeml.Run(mali, options.filename_tree,
                            options=codeml_options,
                            dump=options.dump)

        for x in codeml.mWarnings:
            options.stdlog.write("# PAML warning: %s\n" % x)

        options.stdout.write(result.mResult + "\n")

        if options.filename_log:
            outfile = open(options.filename_log, "w")
        else:
            outfile = options.stdlog
            outfile.write(
                "############################### LOG OUTPUT ############################\n")

        outfile.write(result.mLog + "\n")

        if outfile != options.stdlog:
            outfile.close()

    Experiment.Stop()
