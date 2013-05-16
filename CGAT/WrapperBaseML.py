################################################################################
#   Gene prediction pipeline 
#
#   $Id: WrapperBaseML.py 2765 2009-09-04 16:55:18Z andreas $
#
#   Copyright (C) 2004 Andreas Heger
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
"""
Wrapper for BaseML
==================

"""
import os, sys, string, re, tempfile, subprocess

import Genomics

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

    def __init__(self, message):
        self.message = message

class UsageError(Error):
    """Exception raised for errors while starting

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class BaseMLResult:
    def __init__(self):
        pass

class BaseML:
    
    """
    Parser and outputs configuration files like this::

                seqfile = brown.nuc
                treefile = brown.trees

                outfile = mlb       * main result file
                noisy = 9     * 0,1,2,3: how much rubbish on the screen
                verbose = 0   * 1: detailed output, 0: concise output
                runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                              * 3: StepwiseAddition; (4,5):PerturbationNNI -2: pairwise

                model = 4   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
                            * 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu

                Mgene = 0   * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff

        *        ndata = 5
                clock = 0   * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
            fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below
                kappa = 5  * initial or fixed kappa

            fix_alpha = 0   * 0: estimate alpha; 1: fix alpha at value below
                alpha = 0.5   * initial or fixed alpha, 0:infinity (constant rate)
               Malpha = 0   * 1: different alpha's for genes, 0: one alpha
                ncatG = 5   * # of categories in the dG, AdG, or nparK models of rates
                nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK

                nhomo = 0   * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2
                getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
         RateAncestor = 1   * (0,1,2): rates (alpha>0) or ancestral states

           Small_Diff = 7e-6
            cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
        *        icode = 0  * (with RateAncestor=1. try "GC" in data,model=4,Mgene=4)
        *  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed
                method = 0  * 0: simultaneous; 1: one branch at a time
    """

    mOptions = {
        "noisy"         : 9,
        "verbose"       : 0,
        "runmode"       : 0, 
        "model"         : 4,     
        "Mgene"         : 0,
        "clock"         : 0,
        "fix_kappa"     : 0,
        "kappa"         : 5,
        "fix_alpha"     : 0, 
        "alpha"         : 0.5, 
        "Malpha"        : 0,
        "ncatG"         : 5,
        "nparK"         : 0,
        "nhomo"         : 0,
        "getSE"         : 0,
        "RateAncestor"  : 1,
        "Small_Diff"    : "7e-6",
        "cleandata"     : 1,
        "method"        : 0,
        }
    mExecutable = "baseml"

    def __init__( self ):

        self.mFilenameTree = None
        self.mFilenameSequences = None

    def GetOptions( self ):
        """return options in pretty format"""
        result = ["# Options for BaseML" ]
        for var in self.mOptions.keys():
            result.append("# %-40s: %s" % (var, self.mOptions[var]))
        return string.join(result, "\n")
        
    def SetOption( self, option, value ):
        self.mOptions[option] = value

    def WriteAlignment( self, alignment ):
        """write alignment in Phylip format."""

        outfile = open( self.mTempdir + "/" + self.mFilenameSequences, "w" )

        outfile.write( " %4i %i\n" % (len(alignment), len(alignment[0][1])) )
        for i, s in alignment:
            outfile.write("%s\n%s\n" % (i, s))

        outfile.close()

    def WriteTree( self, alignment ):
        """write tree to file."""

        outfile = open( self.mTempdir + "/" + self.mFilenameTree, "w" )

        if len(alignment) != 2:
            raise UsageError( "only 2 sequences right now!")

        outfile.write("(%s,%s);\n" % (alignment[0][0], alignment[1][0]))
        outfile.close()
        
    def WriteControlFile( self ):
        """write baseml.ctl"""
        
        outfile = open( self.mTempdir + "/baseml.ctl", "w" )

        if self.mFilenameSequences:
            outfile.write( "seqfile = %s\n" % (self.mFilenameSequences) )
        if self.mFilenameTree:
            outfile.write( "treefile = %s\n" % (self.mFilenameTree) )

        outfile.write( "outfile = %s\n" % (self.mFilenameOutput) )

        for o,a in self.mOptions.items():
            outfile.write( "%s = %s\n" % (o, str(a)))
            
        outfile.close()
    
    def Run( self, alignment, dump_result = 0 ):

        self.mTempdir = tempfile.mkdtemp()
        self.mFilenameSequences = "input"
        self.mFilenameTree = "tree"
        self.mFilenameOutput = "output"        
        self.mNumSequences = len(alignment)
        
        self.WriteAlignment( alignment )
        self.WriteTree( alignment )
        self.WriteControlFile()
        
        statement = "cd %s; %s;" % (self.mTempdir, self.mExecutable )

        p = subprocess.Popen( statement , 
                              shell=True, 
                              stdin=subprocess.PIPE, 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, 
                              close_fds=True)

        (file_stdout, file_stdin, file_stderr) = (p.stdin, p.stdout, p.stderr)
        
        file_stdin.close()
        lines = file_stdout.readlines()
        lines_stderr = file_stderr.readlines()
        exit_code = file_stdout.close()
        file_stderr.close()
        file_stdout.close()

        lines = open( "%s/%s" % (self.mTempdir, self.mFilenameOutput), "r").readlines()

        if dump_result: print string.join(lines, "")
        os.system( "rm -rf %s" % self.mTempdir )
        
        return self.ParseResult( lines )

    def ParseResult( self, lines ):
        """parse BaseML output.
        """

        """
BASEML (in paml 3.14, January 2004)  input  HKY85 dGamma (ncatG=5)
ns = 2          ls = 60
# site patterns = 7
   15    7   18   16    1    2    1

sequence 1                           AGCTGCT
sequence 2                           ....ATC




Frequencies..
                                    T      C      A      G
sequence 1                     0.2833 0.3333 0.2500 0.1333
sequence 2                     0.3000 0.3167 0.2667 0.1167

Average                        0.2917 0.3250 0.2583 0.1250

# constant sites:     56 (93.33%)
ln Lmax (unconstrained) = -93.644144

Distances:HKY85 (kappa)  (alpha set at 0.50)
This matrix is not used in later m.l. analysis.

sequence 1
sequence 2         0.0823(999.0000)

====================================================================
BASEML (in paml 3.14, January 2004)  input  JC69 dGamma (ncatG=5)
ns = 2          ls = 60
# site patterns = 7
   15    7   18   16    1    2    1

sequence 1                           AGCTGCT
sequence 2                           ....ATC




Frequencies..
                                    T      C      A      G
sequence 1                     0.2833 0.3333 0.2500 0.1333
sequence 2                     0.3000 0.3167 0.2667 0.1167

Average                        0.2917 0.3250 0.2583 0.1250

# constant sites:     56 (93.33%)
ln Lmax (unconstrained) = -93.644144

Distances: JC69 (alpha set at 0.50)
This matrix is not used in later m.l. analysis.

sequence 1
sequence 2         0.0767
BASEML (in paml 3.14b, May 2005)  input  HKY85 dGamma (ncatG=5)
ns = 2          ls = 314
# site patterns = 16
   77   19   61   66    4    3   20   25    7    6    6    8    5    3    2
   2
seq1                                 AGGTTAACTG CCATGC
seq2                                 .A..GTG.AC ATCCTG
Frequencies..
T      C      A      G
seq1                           0.2548 0.1306 0.3344 0.2803
seq2                           0.2516 0.1242 0.3471 0.2771
Homogeneity statistic: X2 = 0.00044 G = 0.00044
Average                        0.2532 0.1274 0.3408 0.2787
# constant sites:    229 (72.93%)
ln Lmax (unconstrained) = -672.510333
Distances:HKY85 (kappa)  (alpha set at 0.50)
This matrix is not used in later m.l. analysis.
seq1
seq2               0.5990( 4.7311)
TREE #  1:  (1, 2);  MP score: 85.00
lnL(ntime:  2  np:  4):   -679.913467   +0.000000
  3..1     3..2
    0.18769  0.16256  3.25628 999.00000
tree length =   0.35024
(seq1, seq2);
(seq1: 0.187685, seq2: 0.162558);
Detailed output identifying parameters
Parameters (kappa) in the rate matrix (HKY85) (Yang 1994 J Mol Evol 39:105-111):
  3.25628
alpha (gamma, K=5) = 999.00000
r:   0.95611  0.98294  0.99967  1.01660  1.04468
f:   0.20000  0.20000  0.20000  0.20000  0.20000

################
Explanation
TREE #  1:  (1, 2);  MP score: 85.00
lnL(ntime:  2  np:  4):   -679.913467   +0.000000
    3..1     3..2   # branches
    0.18769  0.16256  3.25628 999.00000 # branch1 branch2 kappa alpha
tree length =   0.35024                 # distance
(seq1, seq2);
(seq1: 0.187685, seq2: 0.162558);       # tree format
Detailed output identifying parameters  # 
Parameters (kappa) in the rate matrix (HKY85) (Yang 1994 J Mol Evol 39:105-111):
  3.25628
alpha (gamma, K=5) = 999.00000         # gamma distribution, here: five bins, uniform
r:   0.95611  0.98294  0.99967  1.01660  1.04468
f:   0.20000  0.20000  0.20000  0.20000  0.20000

If alpha is fixed, there are only three values for alpha.

"""
        result = BaseMLResult()        
        # chop
        lines = map(lambda x: x[:-1], lines )
        lines = filter( lambda x: x != "", lines)

        result.mVersion, result.mModel = re.search("BASEML \(in (.+)\)\s+input\s+(\S.+)", lines[0]).groups()

        try:
            if result.mVersion == "paml 3.14, January 2004":
                result.mSitePatterns = map(int, re.split("\s+", lines[0].strip()))
                while not re.match("^Frequencies", lines[0] ): del lines[0]
                result.mFrequencies = {}
                result.mMatrix = {}
                for x in range(self.mNumSequences):
                    id = lines[0][:30].strip()
                    result.mMatrix[id] = {}
                    result.mFrequencies[id] = map(float, re.split( "\s+", lines[0][30:].strip()))
                    del lines[0]
                result.mFrequencies['average'] = map(float, re.split( "\s+", lines[0][30:].strip()))
                result.mConstantSites, result.mConstantSitesPercent = \
                                       re.search( "# constant sites:\s+(\S+)\s\((.+)%\)", lines[0] ).groups()
                result.mLmax  = re.search( "ln Lmax \(unconstrained\) = (\S+)", lines[0] ).groups()
                previous_ids = []
                for x in range(self.mNumSequences):
                    id = lines[0][:17].strip()
                    rest = lines[0][17:].strip()
                    if rest:
                        values = re.findall("[\d.+-]+", rest)                    
                        if self.mOptions['model'] == 4:
                            values = map(float, values[::2])
                        else:
                            values = map(float, values )
                        for i in range(len(previous_ids)):
                            result.mMatrix[previous_ids[i]][id] = values[i]
                            result.mMatrix[id][previous_ids[i]] = values[i]                
                    previous_ids.append(id)
            elif result.mVersion == "paml 3.14b, May 2005":
                # print string.join(lines,"\n")
                while not re.match("^# site patterns", lines[0] ): del lines[0]
                del lines[0]
                result.mSitePatterns = map(int, re.split("\s+", lines[0].strip()))

                result.mFrequencies = {}
                result.mMatrix = {}
                while not re.match("^Frequencies", lines[0] ): del lines[0]
                del lines[:2]
                for x in range(self.mNumSequences):
                    id = lines[0][:30].strip()
                    result.mMatrix[id] = {}
                    result.mFrequencies[id] = map(float, re.split( "\s+", lines[0][30:].strip()))
                    del lines[0]
                result.mX2, result.mG = map(float, re.search( "Homogeneity statistic: X2 = (\S+) G = (\S+)", lines[0] ).groups())
                del lines[0]
                result.mFrequencies['average'] = map(float, re.split( "\s+", lines[0][30:].strip()))
                del lines[0]
                result.mConstantSites, result.mConstantSitesPercent = \
                                       map(float, re.search( "# constant sites:\s+(\S+)\s\((.+)%\)", lines[0] ).groups())

                while lines and not re.match("^TREE", lines[0] ):
                    del lines[0]

                result.mBranch1, result.mBranch2, result.mKappa, result.mAlpha = 0.0, 0.0, 0.0, 0.0                                                        
                ## sometimes, nothing after tree follows,
                ## for example, if there are no variable sites, no estimation possible
                ## take only branch lengths from the next line, because alpha and kappa might or might
                ## not be present
                del lines[:3]
                if lines:
                    data = map(float, re.split("\s+", lines[0].strip()))

                    result.mBranch1, result.mBranch2 = data[:2]

                    if self.mOptions["fix_alpha"] == 0:
                        self.mKappa, self.mAlpha = data[2:4]
                    else:
                        self.mKappa = data[2]
                    
                result.mDistance = result.mBranch1 + result.mBranch2
                keys = result.mMatrix.keys()
                result.mMatrix[keys[0]][keys[1]] = result.mDistance
                result.mMatrix[keys[1]][keys[0]] = result.mDistance            
            else:
                raise ParsingError( "unknown BaseML version" )

        except Exception, inst:
            raise ParsingError( str(inst) + ":" + lines[0] )
            
        return result

if __name__ == "__main__":
    
    baseml = BaseML()

    alignment = ( ("sequence 1", "AAGCTTCACCGGCGCAGTCATTCTCATAATCGCCCACGGACTTACATCCTCATTACTATT"),
                  ("sequence 2", "AAGCTTCACCGGCGCAATTATCCTCATAATCGCCCACGGACTTACATCCTCATTATTATT") )
                  
    baseml.Run( alignment )
    
        
