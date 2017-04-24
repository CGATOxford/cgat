
'''
WrapperIDR.py - wrap idr analysis
=================================

:Tags: Python

Purpose
-------

Wrap the IDR analysis workflow.

For more information, see

https://sites.google.com/site/anshulkundaje/projects/idr


Usage
-----

Documentation
-------------

Code
----

This script requires some R functions within the package
called 

functions-all-clayton-12-13.r

These functions have been copied into a script called
WrapperIDR.r in the scripts directory.

'''

import os
import sys

from CGAT import Experiment as E

from rpy2.robjects import r as R


def runIDR(options, peakfile1, peakfile2):
    '''run IDR analysis.

    This code is taken from the R script

    batch-consistency-analysis.r
    '''

    if options.half_width is not None:
        R.assign("half.width", options.half_width)
    else:
        R('''half.width = NULL''')
    R.assign("overlap.ratio", options.overlap_ratio)
    R.assign("is.broadpeak", options.is_broadpeak)
    R.assign("sig.value", options.signal_value)

    dirname = os.path.dirname(__file__)
    R.source(os.path.join(dirname, "WrapperIDR.r"))

    # read the length of the chromosomes, which will be used to concatenate
    # chr's
    R('''chr.size <- read.table('%s', sep='\t')''' %
      options.filename_chromosome_table)

    output_prefix = options.output_prefix
    output_uri = output_prefix + "-uri.sav"
    output_em = output_prefix + "-em.sav"
    output_overlapped_peaks = output_prefix + "-overlapped-peaks.txt"
    output_peaks_above_idr = output_prefix + "-npeaks-aboveIDR.txt"

    # process data, summit: the representation of the location of summit
    E.info("loading data")
    R('''rep1 <- process.narrowpeak('%(peakfile1)s', chr.size, 
                     half.width=half.width, summit="offset", broadpeak=is.broadpeak)''' % locals())
    R('''rep2 <- process.narrowpeak('%(peakfile2)s', chr.size, 
                     half.width=half.width, summit="offset", broadpeak=is.broadpeak)''' % locals())

    E.info("replicate 1: read %s: %i peaks, %i after filtering" %
           (peakfile1,
            R('''nrow(rep1$data.ori)''')[0],
            R('''nrow(rep1$data.cleaned)''')[0]))
    E.info("replicate 2: read %s: %i peaks, %i after filtering" %
           (peakfile2,
            R('''nrow(rep2$data.ori)''')[0],
            R('''nrow(rep2$data.cleaned)''')[0]))

    E.info("computing correspondence profile (URI)")

    R('''uri.output <- compute.pair.uri(rep1$data.cleaned, rep2$data.cleaned, 
                                        sig.value1=sig.value, sig.value2=sig.value, 
                                        overlap.ratio=overlap.ratio)''')
    E.info("saving correspondence profile to %s" % output_uri)
    R('''save(uri.output, file='%(output_uri)s') ''' % locals())

    E.info("computing EM procedure for inference")
    R('''em.output <- fit.em(uri.output$data12.enrich, fix.rho2=T)''')
    E.info("saving EM to %s" % output_em)
    R('''save(em.output, file='%(output_em)s') ''' % locals())

    # write em output into a file
    # cat(paste("EM estimation for the following files\n", peakfile1, "\n", peakfile2, "\n", sep=""))

    options.stdout.write("em_estimation\n%s\n" %
                         str(R('''em.output$em.fit$para''')))

    # add on 3-29-10
    # output both local idr and IDR
    E.info("writing overlapped peaks to %s" % output_overlapped_peaks)
    R('''idr.local <- 1-em.output$em.fit$e.z''')
    R('''IDR <- c()''')
    R('''o <- order(idr.local)''')
    R('''IDR[o] <- cumsum(idr.local[o])/c(1:length(o))''')
    R('''
    write.out.data <- data.frame(chr1=em.output$data.pruned$sample1[, "chr"],
                                 start1=em.output$data.pruned$sample1[, "start.ori"],
                                 stop1=em.output$data.pruned$sample1[, "stop.ori"],
                                 sig.value1=em.output$data.pruned$sample1[, "sig.value"],
                                 chr2=em.output$data.pruned$sample2[, "chr"],
                                 start2=em.output$data.pruned$sample2[, "start.ori"],
                                 stop2=em.output$data.pruned$sample2[, "stop.ori"],
                                 sig.value2=em.output$data.pruned$sample2[, "sig.value"],
                                 idr.local=1-em.output$em.fit$e.z, IDR=IDR)
    ''')
    R('''write.table(write.out.data, file='%(output_overlapped_peaks)s')''' %
      locals())

    # number of peaks passing IDR range (0.01-0.25)
    E.info("computing number of peaks at various thresholds")
    R('''IDR.cutoff <- seq(0.01, 0.25, by=0.01)''')
    R('''idr.o <- order(write.out.data$idr.local)''')
    R('''idr.ordered <- write.out.data$idr.local[idr.o]''')
    R('''IDR.sum <- cumsum(idr.ordered)/c(1:length(idr.ordered))''')
    R('''
    IDR.count <- c()
    n.cutoff <- length(IDR.cutoff)
    for(i in 1:n.cutoff){
        IDR.count[i] <- sum(IDR.sum <= IDR.cutoff[i])
        }
    ''')

    # write the number of peaks passing various IDR ranges into a file
    E.info(
        "writing number of peaks above IDR cutoffs in range [0.01, 0.25] to %s" % output_peaks_above_idr)
    R('''idr.cut <- data.frame( cutoff=IDR.cutoff, count=IDR.count)''')
    R('''write.table(idr.cut, file='%(output_peaks_above_idr)s', quote=F, 
                     row.names=F, col.names=T, sep='\t')''' % locals())

    R('''mar.mean <- get.mar.mean(em.output$em.fit)''')
    options.stdout.write("marginal mean of two components\n%s\n)" %
                         R('''print(mar.mean)'''))


def plotIDR(output_file, input_prefixes):
    '''create IDR plots.

    This code is taken from the R script

    batch-consistency-plot.r

    within the IDR package.
    '''

    dirname = os.path.dirname(__file__)
    R.source(os.path.join(dirname, "WrapperIDR.r"))

    R('''df.txt = 10''')

    R('''uri.list <- list()
         uri.list.match <- list()
         ez.list <- list()
         legend.txt <- c()
         em.output.list <- list()
         uri.output.list <- list()''')

    npair = len(input_prefixes)
    for x, input_prefix in enumerate(input_prefixes):

        R.load(input_prefix + "-uri.sav")
        R.load(input_prefix + "-em.sav")
        i = x + 1

        R( '''uri.output.list[[%(i)i]] <- uri.output;
              em.output.list[[%(i)i]] <- em.output;
              # reverse =T for error rate;''' % locals())
        R('''
              ez.list[[%(i)i]] <- get.ez.tt.all(em.output, uri.output.list[[%(i)i]]$data12.enrich$merge1,
                                        uri.output.list[[%(i)i]]$data12.enrich$merge2);''' % locals())
        R('''
              # URI for all peaks
              uri.list[[%(i)i]] <- uri.output$uri.n;

              # URI for matched peaks
              uri.match <- get.uri.matched(em.output$data.pruned, df=df.txt);
              uri.list.match[[%(i)i]] <- uri.match$uri.n;
         ''' % locals() )

        legend = "%(i)i = %(input_prefix)s" % locals()
        R('''
              legend.txt[%(i)i] <- '%(legend)s';
        ''' % locals())

    R.pdf(output_file)
    R('''par(mfcol=c(2,3), mar=c(5,6,4,2)+0.1)''')
    R('''plot.uri.group(uri.list, NULL, file.name=NULL, c(1:%(npair)i), title.txt="all peaks");
         plot.uri.group(uri.list.match, NULL, file.name=NULL, c(1:%(npair)i), title.txt="matched peaks");
         plot.ez.group(ez.list, plot.dir=NULL, file.name=NULL, legend.txt=c(1:%(npair)i), y.lim=c(0, 0.6));
         plot(0, 1, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n"); 
         legend(0, 1, legend.txt, cex=0.6);''' % locals())
    R["dev.off"]()


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-o", "--output-prefix", dest="output_prefix", type="string",
                      help="output filename prefix [default=%default].")

    parser.add_option("-c", "--chromosome-table", dest="filename_chromosome_table", type="string",
                      help="filename with tab separated list of chromosome names [default=%default].")

    parser.add_option("--action", dest="action", type="choice",
                      choices=("plot", "run"),
                      help="action to perform [default=%default]")

    parser.add_option("-s", "--signal-value", dest="signal_value", type="string",
                      help="use either p.value or sig.value as ranking measure [default=%default]")

    parser.add_option("-r", "--overlap-ratio", dest="overlap_ratio", type="int",
                      help="a value between 0 and 1 that controls how much two peaks have to overlap to be called as the same [default=%default]")

    parser.set_defaults(
        action="plot",
        output_prefix="output",
        half_width=None,
        overlap_ratio=0,
        is_broadpeak=False,
        signal_value="signal.value",
        filename_chromosome_table="genome_table.txt",
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.action == "plot":
        plotIDR(options.output_prefix + ".pdf", args)
    elif options.action == "run":
        if len(args) != 2:
            raise ValueError("require exactly two replicates")
        runIDR(options, args[0], args[1])

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
