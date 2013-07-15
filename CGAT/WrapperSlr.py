################################################################################
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
#################################################################################
'''
WrapperSlr.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import os, sys, string, re, tempfile, subprocess, optparse, shutil

from types import *

"""Wrapper for CodeSlr
"""

import Experiment
import Mali
import TreeTools

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

    def __str__(self):
        return str(self.message)

class ParsingError(Error):
    """Exception raised for errors while parsing

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message, line):
        self.message = message + " at line " + line

class UsageError(Error):
    """Exception raised for errors while starting

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class Result:
    
    mMaxLineLength = 100
    
    def __init__(self):
        self.mWarnings = []

    def truncateLine( self, line ):
        if len(line) > self.mMaxLineLength:
            return "%s..." %  line[:self.mMaxLineLength]
        else:
            return line
        
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
                        s.append( ("%-40s: %s" % (key, self.truncateLine(str(m[x]))) ))
                elif type(m) in (DictType,):
                    for x, y in m.items():
                        s.append( ("%-40s: %s: %s" % (key, str(x), self.truncateLine(str(y)))))
                else:
                    s.append( ("%-40s: %s" % (key, self.truncateLine(str(m)))))
                    
        return "\n".join(s)

class SlrResult(Result):
    """result of an SLR run."""
    def __init__(self):
        Result.__init__(self)
        self.mSites = []
        self.mNSitesSynonymous = 0
        self.mNSitesGaps = 0
        self.mNSitesSingleChar = 0                
        self.mLog = ""
        self.mResult = ""
    
class SlrResultSite(Result):
    """result with site specific information.

    mResult can be either
    -,--        for negatively selected sites
    0           for neutral sites
    +,++        for positively selected sites
    None        for sites, which can't be estimated: gapped sites/synonymous sites/etc.
    """
    def __init__(self,
                 residue, lneutral, loptimal, w,
                 statistic, pvalue, adjusted_pvalue,
                 result, warning, note, is_synonymous):
        
        Result.__init__(self)        
        self.mResidue = residue
        self.mLikelihoodNeutral = lneutral
        self.mLikelihoodOptimal = loptimal
        self.mOmega = w
        self.mStatistic = statistic
        self.mPValue = pvalue
        self.mAdjustedPValue = adjusted_pvalue

        self.mResult = result
        self.mIsSynonymous = is_synonymous
        self.mHasWarning = warning 
        self.mNote = note

        self.mIsSynonymous = note == "Synonymous"
        
        if self.mResult:
            self.mIsPositive = self.mResult[0] == "+"
            self.mIsNegative = self.mResult[0] == "-"
        else:
            self.mIsPositive = False
            self.mIsNegative = False            
        
    def __str__(self):
        return "\t".join( ( "%i" % self.mResidue,
                            "%s" % self.mResult,
                            "%s" % str(self.mHasWarning),
                            "%5.2f" % self.mOmega,
                            "%5.2f" % self.mPValue,
                            "%5.2f" % self.mAdjustedPValue,                            
                            "%5.2f" % self.mLikelihoodNeutral,
                            "%5.2f" % self.mLikelihoodOptimal,
                            "%5.2f" % self.mStatistic,
                            "%s" % self.mNote) )

    def isSynonymous(self ):
        return self.mIsSynonymous
    
    def isNegative(self, threshold = 0.05, use_adjusted = False):
        """returns true, if site is under negative selection."""
        if not self.mIsNegative:
            return False
            
        if use_adjusted:
            p = self.mAdjustedPValue
        else:
            p = self.mPValue

        return p < threshold 

    def isPositive(self, threshold = 0.05, use_adjusted = False):
        """returns true, if site is under positive selection."""
        
        if not self.mIsPositive:
            return False
            
        if use_adjusted:
            p = self.mAdjustedPValue
        else:
            p = self.mPValue

        return p < threshold 
            
class Slr:
    """
    """
    mExecutable = "Slr"

    def __init__( self ):

        self.mFilenameTree = None
        self.mFilenameSequences = None

        self.mDefaults = { "positive_only": "0",
                           "initial_kappa": "2.0",
                           "initial_omega": "0.1" }

        self.mOptions = {
            "positive_only" : {  "0": "output positve and neutral sites",
                                 "1": "only output positive sites"
                                 },
            "initial_kappa"   : { "#" : "initial value for kappa",
                                  },
            "initial_omega"  : { "#" : "initial value for omega",
                                 },
            }

    def GetOptions( self ):
        """return options in pretty format"""
        result = ["# Options for Slr" ]
        for var in self.mOptions.keys():
            result.append("# %-40s: %s" % (var, self.mOptions[var]))
        return string.join(result, "\n")

    ##------------------------------------------------------------------------
    def WriteAlignment( self, mali ):
        """write alignment in Phylip format."""

        mali.mapIdentifiers( self.mMapOld2New )
        outfile = open( self.mTempdir + "/" + self.mFilenameSequences, "w" )
        mali.writeToFile( outfile, format="phylip" )
        outfile.close()

    ##------------------------------------------------------------------------
    def WriteTree( self, tree ):
        """write tree to file.
        """

        nexus = TreeTools.Newick2Nexus( tree )
        t = nexus.trees[0]
        TreeTools.MapTaxa( t, self.mMapOld2New )
        
        outfile = open( self.mTempdir + "/" + self.mFilenameTree, "w" )
        outfile.write("%i 1\n" % self.mNumSequences )
        outfile.write("%s\n" % TreeTools.Tree2Newick(t))
        outfile.close()

    ##------------------------------------------------------------------------
    def BuildMap( self, mali ):
        """build map of identifiers.
        SLR requires the taxa in the tree to be identical to the
        position (not the name) in the multiple alignment file (and
        not the ids). Thus choose sequential numbering for ids.
        """
        
        self.mMapOld2New = {}
        self.mMapNew2Old = {}
        x = 0
        for id in mali.getIdentifiers():
            x += 1            
            new_id = "%i" % x
            self.mMapOld2New[id] = new_id
            self.mMapNew2Old[new_id] = id

    ##------------------------------------------------------------------------        
    def Run( self, alignment, tree = None,
             dump = 0,
             test = False,
             options = {} ):

        self.mTempdir = tempfile.mkdtemp()
        self.mFilenameSequences = "input"

        self.mFilenameOutput = "output"        
        self.mNumSequences = len(alignment)

        self.BuildMap( alignment )
        
        self.WriteAlignment( alignment )

        self.mNumSequenes = alignment.getLength()
        
        if test:
            print "# temporary directory is %s" % self.mTempdir
            
        if tree:
            self.mFilenameTree = "tree"
            
            ## check what kind of tree is given.
            if type(tree) == StringType:
                t = tree.strip()
                if not (t[0] == "(" and t[-1] in ");"):
                    tree = "".join(open( tree, "r").readlines())
                self.WriteTree( tree )                    
            
        self.mFilenameControl = "%s.ctl" % self.mExecutable

        ## build options
        ## write set options
        set_options = {}
        for key, vv in options.items():

            if key not in self.mOptions:
                raise IndexError, "illegal option %s" % key

            values = vv.strip().split(" ")
            
            descriptions = []
            for value in values:
                if value in self.mOptions[key]:
                    descriptions.append( self.mOptions[key][value] )
                elif "#" in self.mOptions[key]:
                    try:
                        v = float(value)
                    except ValueError:
                        raise ValueError, "not a number: %s for option %s" % (value, key)
                    descriptions.append( self.mOptions[key]["#"] )
                else:
                    raise ValueError, "illegal value %s for option %s: possible values are: %s" % (value, key,
                                                                                                   " ".join(self.mOptions[key].keys()) )

            set_options[key] = vv


        ## write default options
        for key, vv in self.mDefaults.items():
            
            if key in set_options: continue

            if key not in self.mOptions:
                raise IndexError, "illegal option %s" % key

            values = vv.strip().split(" ")
            
            descriptions = []
            for value in values:
                if value in self.mOptions[key]:
                    descriptions.append( self.mOptions[key][value] )
                elif "#" in self.mOptions[key]:
                    try:
                        v = float(value)
                    except ValueError:
                        raise ValueError, "not a number: %s for option %s" % (value, key)
                    descriptions.append( self.mOptions[key]["#"] )
                else:
                    raise ValueError, "illegal value %s for option %s: possible values are: %s" % (value, key,
                                                                                                   " ".join(self.mOptions[key].keys()) )

            set_options[key] = vv

        # write a dummy control file
        outfile = open( self.mTempdir + "/slr.ctl", "w" )
        outfile.write("\n")
        outfile.close()
            
        statement = " ".join( ( self.mExecutable,
                                "-seqfile %s" % self.mFilenameSequences,
                                "-treefile %s" % self.mFilenameTree,
                                "-outfile %s" % self.mFilenameOutput,
                                "-reoptimize 1", 
                                "-positive_only %s" % set_options["positive_only"],
                                "-kappa %s" % set_options["initial_kappa"],
                                "-omega %s" % set_options["initial_omega"] ) )        

        s = subprocess.Popen( statement,
                              shell = True,
                              stdout = subprocess.PIPE,
                              stderr = subprocess.PIPE,
                              cwd = self.mTempdir,
                              close_fds = True)                              

        (out, err) = s.communicate()
        
        if s.returncode != 0:
            raise UsageError, "Error in running %s \n%s\n%s\nTemporary directory in %s" % (statement, err, out, self.mTempdir)

        if dump:
            print "# stdout output of %s:\n%s\n######################################" % (self.mExecutable, out)

        lines = open( "%s/%s" % (self.mTempdir, self.mFilenameOutput), "r").readlines()

        if len(lines) == 0:
            raise UsageError, "Empty result from %s \n%s\n%s\nTemporary directory in %s" % (statement, err, out, self.mTempdir)            

        if dump: 
            print "# result output of %s:\n%s\n######################################" % (self.mExecutable, "".join(lines))        

        if not test:
            shutil.rmtree( self.mTempdir )

        return self.parseOutput( lines, out.split("\n") )

    def nextSection( self, lines):
        if lines:
            while not lines[0]: del lines[0]
        
    def skipSection( self, lines ):

        while lines[0]: del lines[0]
        while not lines[0]: del lines[0]

    def parseOutput( self, lines, lines_log = None ):
        """parse Slr output.

        if lines_log is given, the output from the logfile is parsed
        as well.
        """
        
        result = SlrResult()
        
        # chop
        lines = map(lambda x: x[:-1], lines )
        result.mResult = "\n".join(lines)

        self.parseSites( lines, result )
        
        if lines_log:
            if lines_log[0][-1] == "\n":
                lines_log = map(lambda x: x[:-1].strip(), lines_log )
            result.mLog = "\n".join(lines_log)
            self.parseLog( lines_log, result )

        return result        

    def parseLog( self, lines, result ):
        """parse log file.
        """
        ## ignore everything before the results section
        ## (that includes the supplied parameters)

        while not re.match("\*\*\*", lines[0]):
            del lines[0]
        del lines[0]

        result.mLogLikelihood = float(re.match("# lnL = (\S+)", lines[0]).groups()[0])
        del lines[0]            

        result.mTree = lines[0]
        del lines[0]

        del lines[0]

        if re.search("nan", lines[0]):
            result.mKappa, result.mOmega = 0, 0
        else:
            try:
                result.mKappa, result.mOmega = map(float, re.match("# Kappa\s+=\s+(\S+)\s*Omega\s+=\s+(\S+)", lines[0]).groups())
            except AttributeError:
                raise ParsingError( "parsing error for kappa and omega", lines[0] )
                
        del lines[0]

        if re.search("nan", lines[0]):
            result.mTreeLength, result.mAverageBranchLength, \
            result.mMinBranchLength, result.mMaxBranchLength = 0, 0, 0, 0
        else:
            try:
                result.mTreeLength, result.mAverageBranchLength, \
                                    result.mMinBranchLength, result.mMaxBranchLength = map(float, \
                re.match("# Tree length =\s*(\S+),\s*average branch length\s*=\s*(\S+)\s*\(min=(\S+),\s*max=(\S+)\)", lines[0]).groups())
            except AttributeError:
                raise ParsingError( "parsing error for tree", lines[0] )
        
        del lines[0]

        result.mScaleFactor = float(re.match("# Scaling tree to neutral evolution. Factor =\s*(\S+)", lines[0]).groups()[0])

        while not re.match("# Significance", lines[0]) : del lines[0]
        del lines[0]

        result.mNPositiveSites = {}

        for x in range(4):
            try:
                cat, n1, n2 = re.match("# (.+)\s+(\d+)\s+(\d+)", lines[0]).groups()
            except AttributeError:
                raise ParsingError( "parsing for categories", lines[0] )                
            del lines[0]
            result.mNPositiveSites[cat.strip()] = (int(n1), int(n2))

        result.mNNegativeSites = {}
        
        while not re.match("# Significance", lines[0]) : del lines[0]
        del lines[0]

        for x in range(4):
            cat, n1, = re.match("# (.+)\s+(\d+)", lines[0]).groups()
            del lines[0]
            result.mNNegativeSites[cat.strip()] = int(n1)
        
    def parseSites( self, lines, result ):
        """parse sites file.

        seems to be fixed with format. Note, there seems
        to be a warning field between Result and Note without
        a heading.
        """

        fields = (0, 5, 14, 23, 32, 41, 52, 63, 69, 70)
        
        for line in lines[1:]:
            data = []

            for x in range( len(fields) - 1):
                data.append( line[fields[x]:fields[x+1]] )
            if len(line) >= fields[-1]:
                data.append( line[fields[-1]:] )
            else:
                data.append( "")

            data = map( lambda x: x.strip(), data)
            ## blank out synonymous sites and those with all gaps
            is_synonymous = False
            if data[9] in ("Synonymous", "All gaps", "Single char"):
                if data[9] == "Synonymous":
                    result.mNSitesSynonymous += 1
                    is_synonymous = True
                elif data[9] == "All gaps":
                    result.mNSitesGaps += 1
                    data[7] = None                    
                elif data[9] == "Single char":
                    result.mNSitesSingleChar += 1
                    data[7] = None
            elif data[7] == "":
                data[7] = "0"
                
            site = SlrResultSite(
                residue = int(data[0]),
                lneutral = float(data[1]),
                loptimal = float(data[2]),
                w = float(data[3]),
                statistic = float(data[4]),
                pvalue = float(data[5]),
                adjusted_pvalue = float(data[6]),
                result = data[7],
                warning = data[8] == "!",
                note = data[9],
                is_synonymous = is_synonymous)
            
            result.mSites.append( site)

def getOptions( options ):
    """translate command line options to PAML options."""

    slr_options = {}
    
    if options.omega != None:
        slr_options[ "initial_omega" ] = str(options.omega)
    if options.kappa != None:
        slr_options[ "initial_kappa" ] = str(options.kappa)
    if options.positive_only != None:
        slr_options[ "positive_only" ] = str(options.positive_only)

    return slr_options

if __name__ == "__main__":
    
    parser = E.OptionParser( version = "%prog version: $Id: WrapperSlr.py 2781 2009-09-10 11:33:14Z andreas $" )

    parser.add_option("-t", "--filename-tree", dest="filename_tree", type="string",
                      help="filename with tree information."  )
    parser.add_option("-i", "--filename-sequences", dest="filename_sequences", type="string",
                      help="filename with sequences."  )
    parser.add_option("-p", "--parse-output", dest="parse_output", action="store_true",
                      help="parse output")
    parser.add_option("-l", "--filename-log", dest="filename_log", type="string",
                      help="filename for logging information."  )
    parser.add_option( "--set-omega", dest="omega", type="float",
                      help="initial omega value.")
    parser.add_option( "--set-kappa", dest="kappa", type="float",
                      help="initial kappa value.")
    parser.add_option( "--set-positive_only", dest="positive_only", type="int",
                      help="only output positive sites.")
    parser.add_option( "--dump", dest="dump", action="store_true",
                      help="dump output."  )

    parser.set_defaults(
        analysis = None,
        filename_tree = None,
        filename_sequences = "input",
        multiple_genes = False,
        parse_output = False,
        filename_log = None,
        omega = None,
        kappa = None,
        positive_only = None,
        dump = False,
        )

    (options, args) = Experiment.Start( parser )

    slr = Slr()

    slr_options = getOptions(options)

    mali = Mali.Mali()
    mali.readFromFile( open(options.filename_sequences, "r"),
                       format = "fasta")

    if not options.filename_tree:
        raise "please supply a tree."
    
    result = slr.Run( mali, options.filename_tree,
                      dump = options.dump )

    options.stdout.write( result.mResult + "\n" )

    if options.filename_log:
        outfile = open(options.filename_log, "w" )
    else:
        outfile = options.stdlog
        
    outfile.write( result.mLog + "\n" )

    if outfile != options.stdlog:
        outfile.close()
    
    Experiment.Stop()
