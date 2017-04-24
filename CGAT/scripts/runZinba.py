'''
runZinba.py - wrap zinba peak caller
========================================

:Tags: Python

Purpose
-------

Runs the zinba peak caller.

Note that for zinba to work, a mappability file needs to have
been created.

The output bed files contain the P-value as their score field.

Usage
-----

Documentation
-------------

Requirements:

* zinba >= 2.01

Code
----

'''

import os
import sys
import tempfile
import shutil
import subprocess

from rpy2.robjects import r as R

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.BamTools as BamTools


def bamToBed(infile, outfile,
             min_insert_size=0,
             max_insert_size=1000):
    '''convert bam to bed with bedtools.'''
    scriptsdir = "/ifs/devel/andreas/cgat/scripts"

    if BamTools.isPaired(infile):
        # output strand as well
        statement = ['cat %(infile)s '
                     '| python %(scriptsdir)s/bam2bed.py '
                     '--merge-pairs '
                     '--min-insert-size=%(min_insert_size)i '
                     '--max-insert-size=%(max_insert_size)i '
                     '--log=%(outfile)s.log '
                     '--bed-format=6 '
                     '> %(outfile)s' % locals()]
    else:
        statement = "bamToBed -i %(infile)s > %(outfile)s" % locals()

    E.debug("executing statement '%s'" % statement)

    retcode = subprocess.call(statement,
                              cwd=os.getcwd(),
                              shell=True)
    if retcode < 0:
        raise OSError("Child was terminated by signal %i: \n%s\n" %
                      (-retcode, statement))

    return outfile


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-f", "--input-format", dest="input_format",
                      type="choice",
                      choices=("bed", "bam"),
                      help="input file format [default=%default].")

    parser.add_option("-s", "--fragment-size", dest="fragment_size",
                      type="int",
                      help="fragment size, used for the extension parameter "
                      "in Zinba [default=%default].")

    parser.add_option("-m", "--zinba-mappability-dir", dest="mappability_dir",
                      type="string",
                      help="mappability_dir [default=%default].")

    parser.add_option("-b", "--bit-file", dest="bit_filename",
                      type="string",
                      help="2bit genome filename [default=%default].")

    parser.add_option("-c", "--control-filename", dest="control_filename",
                      type="string",
                      help="filename of input/control data in bed format "
                      "[default=%default].")

    parser.add_option("-i", "--zinba-index-dir", dest="index_dir", type="string",
                      help="index directory [default=%default].")

    parser.add_option("-t", "--threads", dest="threads", type="int",
                      help="number of threads to use [default=%default].")

    parser.add_option("-q", "--fdr-threshold", dest="fdr_threshold",
                      type="float",
                      help="fdr threshold [default=%default].")

    parser.add_option("-a", "--zinba-alignability-threshold",
                      dest="alignability_threshold", type="int",
                      help="alignability threshold [default=%default].")

    parser.add_option("-p", "--aggregate-by-contig", dest="per_contig",
                      action="store_true",
                      help="run analysis per chromosome [default=%default]")

    parser.add_option("-w", "--temp-dir", dest="tempdir", type="string",
                      help="use existing directory as temporary directory "
                      "[default=%default].")

    parser.add_option("--keep-temp", dest="keep_temp", action="store_true",
                      help="keep temporary directory [default=%default]")

    parser.add_option("--action", dest="action", type="choice",
                      choices=("full", "count", "predict", "model"),
                      help="action to perform [default=%default]")

    parser.add_option("--zinba-improvement", dest="improvement", type="float",
                      help="relative improvement of likelihood until "
                      "convergence [default=%default]")

    parser.add_option("--min-insert-size", dest="min_insert_size", type="int",
                      help="minimum insert size for paired end data "
                      "[default=%default]")

    parser.add_option("--max-insert-size", dest="max_insert_size", type="int",
                      help="maximum insert size for paired end data "
                      "[default=%default]")

    parser.set_defaults(
        input_format="bed",
        fragment_size=200,
        mappability_dir=None,
        threads=1,
        alignability_threshold=1,
        bit_filename=None,
        fdr_threshold=0.05,
        tempdir=None,
        winsize=250,
        offset=125,
        cnvWinSize=1e+05,
        cnvOffset=2500,
        per_contig=False,
        keep_temp=False,
        min_insert_size=0,
        max_insert_size=1000,
        filelist="files.list",
        selectchr="chr19",
        action="full",
        improvement=0.00001,
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if len(args) != 2:
        raise ValueError(
            "please specify a filename with sample data and an output file")

    filename_sample, filename_output = args[0], args[1]

    filename_sample = os.path.abspath(filename_sample)
    filename_output = os.path.abspath(filename_output)

    if options.control_filename:
        filename_control = os.path.abspath(options.control_filename)
    else:
        filename_control = None

    # load Zinba
    R.library('zinba')

    if not options.tempdir:
        tmpdir = tempfile.mkdtemp()
    else:
        tmpdir = options.tempdir

    E.info("changing to temporary directory %s" % tmpdir)
    os.chdir(tmpdir)

    if options.input_format == "bam":
        E.info("converting bam files to bed")
        if not os.path.exists(os.path.join(tmpdir, "sample.bed")):
            filename_sample = bamToBed(
                filename_sample,
                os.path.join(tmpdir, "sample.bed"))
        else:
            E.info("using existing file %(tmpdir)s/sample.bed" %
                   locals())
            filename_sample = os.path.join(
                tmpdir, "sample.bed")
        if filename_control:
            if not os.path.exists(os.path.join(tmpdir, "control.bed")):
                filename_control = bamToBed(
                    filename_control,
                    os.path.join(tmpdir, "control.bed"))
            else:
                E.info("using existing file %(tmpdir)s/control.bed" %
                       locals())
                filename_control = os.path.join(
                    tmpdir, "control.bed")

    fragment_size = options.fragment_size
    threads = options.threads
    bit_filename = options.bit_filename
    mappability_dir = options.mappability_dir
    fdr_threshold = options.fdr_threshold
    tol = options.improvement

    contigs = E.run(
        "twoBitInfo %(bit_filename)s %(tmpdir)s/contig_sizes" % locals())
    contig2size = dict(
        [x.split() for x in IOTools.openFile(
            os.path.join(tmpdir, "contig_sizes"))])

    outdir = filename_output + "_files"
    E.info('saving intermediate results in %s' % outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    filelist = os.path.join(outdir, filename_output + ".list")
    modelfile = os.path.join(outdir, filename_output + ".model")
    winfile = os.path.join(outdir, filename_output + ".wins")
    winSize = 250
    offset = 125
    cnvWinSize = 100000
    cnvOffset = 0
    winGap = 0
    peakconfidence = 1.0 - fdr_threshold
    selectchr = options.selectchr

    if not os.path.exists(os.path.join(tmpdir, "basecount")):
        E.info("computing counts")

        R('''basealigncount(inputfile='%(filename_sample)s',
        outputfile='%(tmpdir)s/basecount',
        extension=%(fragment_size)i,
        filetype='bed',
        twoBitFile='%(bit_filename)s' )
        ''' % locals())
    else:
        E.info("using existing counts")

    # tried incremental updates
    # for contig, size in contig2size.iteritems():
    #     for size in
    #     fn = os.path.join( tmpdir, "sample_%(contig)s_win%(size)ibp_offset(offset)ibp.txt" % locals() )
    if options.action == "count":

        E.info("computing window counts only - saving results in %s" % outdir)
        R('''buildwindowdata(
                     seq='%(filename_sample)s',
                     align='%(mappability_dir)s',
                     input='%(filename_control)s',
                     twoBit='%(bit_filename)s',
                     winSize=%(winSize)i,
                     offset=%(offset)i,
                     cnvWinSize=%(cnvWinSize)i,
                     cnvOffset=%(cnvOffset)i,
                     filelist='%(filelist)s',
                     filetype='bed',
                     extension=%(fragment_size)s,
                     outdir='%(outdir)s/') ''' % locals())

    elif options.action == "model":

        # The important option is buildwin = 0
        # parameterized for broad == FALSE and input present
        # see zinba.R
        # model selection only on chr19.
        R('''run.zinba( 
                filelist='%(filelist)s',
                formula=NULL,formulaE=NULL,formulaZ=NULL,
                outfile='%(filename_output)s',
                seq='%(filename_sample)s',
                input='%(filename_control)s',
                filetype='bed',  
                align='%(mappability_dir)s',
                twoBit='%(bit_filename)s',
                extension=%(fragment_size)s,
                winSize=%(winSize)i,
                offset=%(offset)i,
                cnvWinSize=%(cnvWinSize)i,
                cnvOffset=%(cnvOffset)i,
                basecountfile='%(tmpdir)s/basecount',
                buildwin=0,
                threshold=%(fdr_threshold)f,
                pquant=1,
                peakconfidence=%(peakconfidence)f,
                winGap=%(winGap)i,
                tol=%(tol)f,
                initmethod="count",
                method="mixture",
                numProc=%(threads)i,
                printFullOut=1,
                interaction=FALSE,
                selectmodel=TRUE,
                selectchr='%(selectchr)s',
                selectcovs=c("input_count"),
                selecttype="complete",
                FDR=TRUE)''' % locals())

    elif options.action == "predict":

        # The important option is buildwin = 0 and selectmodel = FALSE
        # parameterized for broad == FALSE and input present
        # see zinba.R
        # model selection only on chr19.
        if not os.path.exists(modelfile):
            raise OSError("model file %s does not exist" % modelfile)

        E.info("reading model from %s" % modelfile)

        R('''
        final=read.table('%(modelfile)s', header=T, sep="\t")
        final=final[final$fail==0,]
        bestBIC=which.min(final$BIC)
        formula=as.formula(paste("exp_count~",final$formula[bestBIC]))
        formulaE=as.formula(paste("exp_count~",final$formulaE[bestBIC]))
        formulaZ=as.formula(paste("exp_count~",final$formulaZ[bestBIC]))
        cat("Background formula is:\n\t")
        print(formula)
        cat("Enrichment formula is:\n\t")
        print(formulaE)
        cat("Zero-inflated formula is:\n\t")
        print(formulaE)
        ''' % locals())

        E.info("predicting peaks")

        R('''run.zinba(
                filelist='%(filelist)s',
                outfile='%(filename_output)s',
                seq='%(filename_sample)s',
                input='%(filename_control)s',
                filetype='bed',
                align='%(mappability_dir)s',
                twoBit='%(bit_filename)s',
                extension=%(fragment_size)s,
                winSize=%(winSize)i,
                offset=%(offset)i,
                cnvWinSize=%(cnvWinSize)i,
                cnvOffset=%(cnvOffset)i,
                basecountfile='%(tmpdir)s/basecount',
                buildwin=0,
                threshold=%(fdr_threshold)f,
                pquant=1,
                winGap=%(winGap)i,
                initmethod="count",
                selectchr='%(selectchr)s',
                tol=%(tol)f,
                method="mixture",
                numProc=%(threads)i,
                printFullOut=1,
                interaction=FALSE,
                selectmodel=FALSE,
                formula=formula,
                formulaE=formulaE,
                formulaZ=formulaZ,
                peakconfidence=%(peakconfidence)f,
                FDR=TRUE)''' % locals())

    elif options.action == "per_contig":

        E.info("processing per chromosome")
        for contig, size in contig2size.items():
            if contig not in ("chr16",):
                continue

            E.info("processing contig %s" % contig)
            filename_sample_contig = filename_sample + "_%s" % contig
            filename_control_contig = filename_control + "_%s" % contig
            if not os.path.exists(filename_output + "_files"):
                os.mkdir(filename_output + "_files")
            filename_output_contig = os.path.join(
                filename_output + "_files", contig)
            filename_basecounts_contig = os.path.join(
                tmpdir, "basecount_%s" % contig)

            E.run(
                "grep %(contig)s < %(filename_sample)s > %(filename_sample_contig)s" % locals())
            E.run(
                "grep %(contig)s < %(filename_control)s > %(filename_control_contig)s" % locals())

            if not os.path.exists(filename_basecounts_contig):
                E.info("computing counts")

                R('''basealigncount( inputfile='%(filename_sample_contig)s',
                                  outputfile='%(filename_basecounts_contig)s',
                                  extension=%(fragment_size)i,
                                  filetype='bed',
                                  twoBitFile='%(bit_filename)s' )
                                  ''' % locals())
            else:
                E.info("using existing counts")

            # run zinba, do not build window data
            R('''zinba( refinepeaks=1,
            seq='%(filename_sample_contig)s',
            input='%(filename_control_contig)s',
            filetype='bed',
            align='%(mappability_dir)s',
            twoBit='%(bit_filename)s',
            outfile='%(filename_output_contig)s',
            extension=%(fragment_size)s,
            basecountfile='%(filename_basecounts_contig)s',
            numProc=%(threads)i,
            threshold=%(fdr_threshold)f,
            broad=FALSE,
            printFullOut=0,
            interaction=FALSE,
            mode='peaks',
            FDR=TRUE) ''' % locals())
    elif options.action == "full":

        # run zinba, build window data and refine peaks

        # Note that zinba() uses 'chr22' to select model
        # which is not present in mouse. So call run.zinba
        # directly.
        R('''run.zinba(
        refinepeaks=1,
        buildwin=1,
        seq='%(filename_sample)s',
        input='%(filename_control)s',
        filetype='bed',
        align='%(mappability_dir)s',
        twoBit='%(bit_filename)s',
        outfile='%(filename_output)s',
        extension=%(fragment_size)s,
        winSize=%(winSize)i,
        offset=%(offset)i,
        basecountfile='%(tmpdir)s/basecount',
        numProc=%(threads)i,
        threshold=%(fdr_threshold)f,
        pquant=1,
        winGap=%(winGap)i,
        selectchr='%(selectchr)s',
        interaction=FALSE,
        method="mixture",
        cnvWinSize=%(cnvWinSize)i,
        cnvOffset=%(cnvOffset)i,
        selectmodel=TRUE,
        selectcovs=c("input_count"),
        selecttype="complete",
        initmethod="count",
        printFullOut=1,
        diff=0,
        pWinSize=200,
        peakconfidence=%(peakconfidence)f,
        FDR=TRUE) ''' % locals())

    if not (options.tempdir or options.keep_temp):
        shutil.rmtree(tmpdir)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
