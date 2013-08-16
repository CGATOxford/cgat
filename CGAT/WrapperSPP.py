################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
WrapperSPP.py - wrap spp peak caller
========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Runs the spp peak caller.

The workflow follows the tutorial at:

http://compbio.med.harvard.edu/Supplements/ChIP-seq/tutorial.html

Usage
-----

Documentation
-------------

Code
----

'''

import os, sys, re, optparse, tempfile, shutil, subprocess
import collections

import Experiment as E
import IOTools

## for zinba
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError

def bamToBed( infile, outfile ):
    '''convert bam to bed with bedtools.'''

    statement = "bamToBed -i %(infile)s > %(outfile)s" % locals()

    E.debug( "executing statement '%s'" % statement )

    retcode = subprocess.call(  statement,
                                cwd = os.getcwd(), 
                                shell = True )
    if retcode < 0:
        raise OSError( "Child was terminated by signal %i: \n%s\n" % (-retcode, statement ))

    return outfile

ZinbaPeak = collections.namedtuple( "ZinbaPeak", "contig unrefined_start unrefined_end strand posterior summit height refined_start refined_end median fdr" )

def iteratePeaks( infile ):
    '''iterate of zinba peaks in infile.'''
    
    for line in infile:

        if line.startswith("#"): continue
        if line.startswith("PEAKID\tChrom"): continue
        # skip empty lines
        if line.startswith("\n"): continue

        data = line[:-1].split("\t")

        if len(data) != 12:
            raise ValueError( "could not parse line %s" % line )

        # I assume these are 1-based coordinates
        data[2] = max(int(data[2]) - 1, 0)
        # end
        data[3] = int(data[3])
        # posterior
        data[5] = float(data[5])
        # summit
        data[6] = max(int(data[6]) - 1, 0)
        # height
        data[7] = int(data[7])
        # refined_start
        data[8] = max(int(data[8]) - 1, 0)
        # end
        data[9] = int(data[9])
        # median
        data[10] = int(data[10])
        # qvalue
        data[11] = float(data[11])

        yield ZinbaPeak._make( data[1:] )

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-f", "--input-format", dest="input_format", type="choice",
                      choices = ("bam",),
                      help="input file format [default=%default]."  )
    
    parser.add_option("-w", "--window-size", dest="window_size", type="int",
                      help="window size [default=%default]."  )

    parser.add_option("-c", "--control-filename", dest="control_filename", type="string",
                      help="filename of input/control data in bed format [default=%default]."  )

    parser.add_option("-t", "--threads", dest="threads", type="int",
                      help="number of threads to use [default=%default]."  )

    parser.add_option("-q", "--fdr-threshold", dest="fdr_threshold", type="float",
                      help="fdr threshold [default=%default]."  )
 
    parser.add_option("-z", "--z-threshold", dest="z_threshold", type="float",
                      help="z threshold [default=%default]."  )
 
    parser.add_option( "--bin", dest="bin", type="int",
                       help="bin tags within the specified number of basepairs to speed up calculation;"
                       " increasing bin size decreases the accuracy of the determined parameters [default=%default]")

    parser.add_option( "--srange-min", dest="srange_min", type="float",
                       help = "srange gives the possible range for the size of the protected region;"
                               " srange should be higher than tag length; making the upper boundary too high"
                               " will increase calculation time [%default]" )

    parser.add_option( "--srange-max", dest="srange_max", type="float",
                       help = "srange gives the possible range for the size of the protected region;"
                               " srange should be higher than tag length; making the upper boundary too high"
                               " will increase calculation time [%default]" )
    
    parser.set_defaults(
        input_format = "bam",
        threads = 1,
        fdr_threshold = 0.05,
        window_size = 1000,
        offset = 125,
        srange_min=50,
        srange_max=500,
        bin=5,
        z_threshold = 3,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if len(args) != 2:
        raise ValueError("please specify a filename with sample data and an output file")

    filename_sample, filename_output = args[0], args[1]
    filename_control = options.control_filename
    
    # load Zinba
    R.library( 'spp' )
    R.library( 'snow' )
    
    # read data
    E.info( "reading data" )
    R('''chip.data <- read.bam.tags('%s')''' % filename_sample)
    R('''input.data <- read.bam.tags('%s')''' % filename_control)
    R('''cluster = makeCluster( %i )''' % (options.threads) )

    E.info( "computing binding characteristics" )
    # get binding info from cross-correlation profile

    # srange gives the possible range for the size of the protected region;
    # srange should be higher than tag length; making the upper boundary too high will increase calculation time

    # bin - bin tags within the specified number of basepairs to speed up calculation;
    # increasing bin size decreases the accuracy of the determined parameters
    srange_min, srange_max = options.srange_min, options.srange_max
    bin = options.bin
    R('''binding.characteristics <- get.binding.characteristics(chip.data,
                                          srange=c(%(srange_min)i,%(srange_max)i),
                                          bin=%(bin)s,
                                          cluster=cluster);''' % locals())
    # print out binding peak separation distance
    options.stdout.write( "shift\t%i\n" % R('''binding.characteristics$peak$x''')[0])

    ##################################################
    ##################################################
    ##################################################
    E.info( "plot cross correlation profile" )
    # plot cross-correlation profile
    R('''pdf(file="%s.crosscorrelation.pdf",width=5,height=5)''' % filename_output )
    R('''par(mar = c(3.5,3.5,1.0,0.5), mgp = c(2,0.65,0), cex = 0.8);''')
    R('''plot(binding.characteristics$cross.correlation,type='l',xlab="strand shift",ylab="cross-correlation");''')
    R('''abline(v=binding.characteristics$peak$x,lty=2,col=2)''')
    R('''dev.off();''')

    E.info( "selecting informative tags based on the binding characteristics" )
    # select informative tags based on the binding characteristics
    R('''chip.data <- select.informative.tags(chip.data,binding.characteristics);''')
    R('''input.data <- select.informative.tags(input.data,binding.characteristics);''')

    E.info( "outputting broad peaks" )
    window_size, z_threshold = options.window_size, options.z_threshold
    R('''broad.clusters <- get.broad.enrichment.clusters(chip.data,input.data,
                                           window.size=%(window_size)i,
                                           z.thr=%(z_threshold)f,
                                           tag.shift=round(binding.characteristics$peak$x/2))''' % locals())
    # write out in broadPeak format
    R('''write.broadpeak.info(broad.clusters,"%s.broadpeak.txt")''' % filename_output )

    # binding detection parameters
    # desired FDR (1%). Alternatively, an E-value can be supplied to the method calls below instead of the fdr parameter
    # the binding.characteristics contains the optimized half-size for binding detection window
    R('''detection.window.halfsize <- binding.characteristics$whs;''')

    # determine binding positions using wtd method
    E.info( "determining binding positions using wtd method" )
    fdr = options.fdr_threshold
    R('''bp <- find.binding.positions(signal.data=chip.data,control.data=input.data,
                            fdr=%(fdr)f,whs=detection.window.halfsize,cluster=cluster)''' % locals())
    options.stdout.write( "detected_peaks\t%i\n" % R('''sum(unlist(lapply(bp$npl,function(d) length(d$x))))''')[0])

    # output detected binding positions
    R('''output.binding.results(bp,"%s.summit.txt");''' % filename_output)

    R('''bp <- add.broad.peak.regions(chip.data,input.data,bp,
                        window.size=%(window_size)i,z.thr=%(z_threshold)f)''' % locals())
    # output using narrowPeak format
    R('''write.narrowpeak.binding(bp,"%s.narrowpeak.txt")''' % filename_output )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

