##########################################################################
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
##########################################################################
'''WrapperMEDIPS.py - wrap MEDIPS methylation analysis
===================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

perform Medip-Seq data analysis using the MEDIPS R package.

Note that the appropriate UCSC genome file has to have been downloaded
previously in Bioconductor::

   source("http://bioconductor.org/biocLite.R")
   biocLite("BSgenome.Hsapiens.UCSC.hg19")
   biocLite("BSgenome.Rnorvegicus.UCSC.rn5")
   ...

Usage
-----

Documentation
-------------

Code
----

'''

import os
import sys
import tempfile
import shutil

from CGAT import Experiment as E
from CGAT import IOTools as IOTools
from CGAT import IndexedFasta as IndexedFasta

from rpy2.robjects import r as R


def compress(infile):
    '''gzip infile'''

    statement = "gzip -f %(infile)s" % locals()

    E.debug("executing statement '%s'" % statement)

    return E.run(statement)


def bigwig(infile, contig_sizes):
    '''convert infile to bigwig file'''

    if infile.endswith(".wig"):
        outfile = infile[:-4] + ".bigwig"
    else:
        outfile = infile + ".bigwig"

    tmp, filename_sizes = tempfile.mkstemp()

    os.write(tmp, "\n".join(["\t".join(map(str, x))
             for x in contig_sizes.iteritems()]))
    os.close(tmp)

    statement = "wigToBigWig -clip %(infile)s %(filename_sizes)s %(outfile)s " % locals()

    E.debug("executing statement '%s'" % statement)

    if E.run(statement):
        os.unlink(infile)

    os.unlink(filename_sizes)


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(
        version="%prog version: $Id",
        usage=globals()["__doc__"])

    parser.add_option("-u", "--ucsc-genome", dest="ucsc_genome", type="string",
                      help="UCSC genome identifier [default=%default].")

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default].")

    parser.add_option("-e", "--extend", dest="extension", type="int",
                      help="extend tags by this number of bases "
                      "[default=%default].")

    parser.add_option("-s", "--shift", dest="shift", type="int",
                      help="shift tags by this number of bases "
                      "[default=%default].")

    parser.add_option("-b", "--bin-size", dest="bin_size", type="int",
                      help="bin size of genome vector [default=%default].")

    parser.add_option("-l", "--fragment-length", dest="fragment_length",
                      type="int",
                      help="bin size of genome vector [default=%default].")

    parser.add_option("--saturation-iterations",
                      dest="saturation_iterations", type="int",
                      help="iterations for saturation analysis "
                      "[default=%default].")

    parser.add_option("-t", "--toolset", dest="toolset", type="choice",
                      action="append",
                      choices=("saturation", "coverage", "enrichment",
                               "dmr", "rms", "rpm", "all"),
                      help = "actions to perform [default=%default].")

    parser.add_option("-w", "--bigwig", dest="bigwig", action="store_true",
                      help="store wig files as bigwig files - requires a "
                      "genome file [default=%default]")

    parser.add_option("--treatment", dest="treatment_files", type="string",
                      action="append",
                      help="BAM files for treatment. At least one is required "
                      "[%default]")

    parser.add_option("--control", dest="control_files", type="string",
                      action="append",
                      help="BAM files for control for differential "
                      "methylation analysis. Optional [%default].")

    parser.add_option("--input", dest="input_files", type="string",
                      action="append",
                      help="BAM files for input correction. "
                      "Optional [%default].")

    parser.add_option("--is-not-medip",
                      dest="is_medip", action="store_false",
                      help="data is not MeDIP data and is not expected "
                      "to fit the calibration model. No CpG "
                      "density normalized rms data is computed"
                      "[default=%default].")

    parser.set_defaults(
        input_format="bam",
        ucsc_genome="Hsapiens.UCSC.hg19",
        genome_file=None,
        extend=0,
        shift=0,
        bin_size=50,
        window_size=300,
        saturation_iterations=10,
        fragment_length=700,
        toolset=[],
        bigwig=False,
        treatment_files=[],
        control_files=[],
        input_files=[],
        is_medip=True,
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv, add_output_options=True)

    if len(options.treatment_files) < 1:
        raise ValueError("please specify a filename with sample data")

    if options.bigwig and not options.genome_file:
        raise ValueError("please provide a genome file when outputting bigwig")

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta(options.genome_file)
        contig_sizes = fasta.getContigSizes()

    if len(options.toolset) == 0:
        options.toolset = ["all"]

    do_all = "all" in options.toolset

    # load MEDIPS
    R.library('MEDIPS')
    genome_file = 'BSgenome.%s' % options.ucsc_genome
    R.library(genome_file)

    bin_size = options.bin_size
    window_size = options.window_size
    extend = options.extend
    shift = options.shift
    fragment_length = options.fragment_length
    saturation_iterations = options.saturation_iterations
    uniq = "TRUE"

    if "saturation" in options.toolset or do_all:
        E.info("saturation analysis")
        for fn in options.treatment_files + options.control_files:
            R('''sr = MEDIPS.saturation(
            file='%(fn)s',
            BSgenome='%(genome_file)s',
            shift=%(shift)i,
            extend=%(extend)i,
            window_size=%(window_size)i,
            uniq=%(uniq)s,
            nit = %(saturation_iterations)i,
            nrit = 1)''' % locals())

            R.png(E.getOutputFile("%s_saturation.png" % fn))
            R('''MEDIPS.plotSaturation(sr)''')
            R('''dev.off()''')
            R('''write.table(sr$estimation, file ='%s', sep='\t')''' %
              E.getOutputFile("%s_saturation_estimation.tsv" % fn))

            outfile = IOTools.openFile("%s_saturation.tsv" % fn, "w")
            outfile.write("category\tvalues\n")
            outfile.write(
                "estimated_correlation\t%s\n" %
                ",".join(["%f" % x for x in R('''sr$maxEstCor''')]))
            outfile.write(
                "true_correlation\t%s\n" %
                ",".join(["%f" % x for x in R('''sr$maxTruCor''')]))
            outfile.write(
                "nreads\t%s\n" %
                ",".join(["%i" % x for x in R('''sr$numberReads''')]))
            outfile.close()

    if "coverage" in options.toolset or do_all:
        E.info("CpG coverage analysis")
        for fn in options.treatment_files + options.control_files:
            R('''cr = MEDIPS.seqCoverage(
            file='%(fn)s',
            BSgenome='%(genome_file)s',
            pattern='CG',
            shift=%(shift)i,
            extend=%(extend)i,
            uniq=%(uniq)s)''' % locals())

            R.png(E.getOutputFile("%s_cpg_coverage_pie.png" % fn))
            R('''MEDIPS.plotSeqCoverage(seqCoverageObj=cr,
            type = "pie", cov.level = c(0, 1, 2, 3, 4, 5))''')
            R('''dev.off()''')

            R.png(E.getOutputFile("%s_cpg_coverage_hist.png" % fn))
            R('''MEDIPS.plotSeqCoverage(seqCoverageObj=cr,
            type = "hist", t=15)''')
            R('''dev.off()''')

            # note: this file is large
            R('''write.table(cr$cov.res, file=gzfile('%s','w'),
            sep='\t')''' %
              E.getOutputFile("%s_saturation_coveredpos.tsv.gz" % fn))

    if 'enrichment' in options.toolset or do_all:
        E.info("CpG enrichment analysis")
        outfile = IOTools.openFile(E.getOutputFile("enrichment.tsv.gz"), "w")
        slotnames = (("regions.CG", "regions_CG", "%i"),
                     ("regions.C", "regions_C", "%s"),
                     ("regions.G", "regions_G", "%f"),
                     ("regions.relH", "regions_relH", "%i"),
                     ("regions.GoGe", "regions_GoGe", "%i"),
                     ("genome.CG", "genome_CG", "%s"),
                     ("genome.C", "genome_C", "%s"),
                     ("genome.G", "genome_G", "%i"),
                     ("genome.relH", "genome_relH", "%i"),
                     ("enrichment.score.relH", "enrichment_relH", "%s"),
                     ("enrichment.score.GoGe", "enrichment_GoGe", "%s"))

        outfile.write("\t".join(['sample'] +
                                [x[1] for x in slotnames]) + "\n")
        for fn in options.treatment_files + options.control_files:
            R('''ce = MEDIPS.CpGenrich(
            file='%(fn)s',
            BSgenome='%(genome_file)s',
            shift=%(shift)i,
            extend=%(extend)i,
            uniq=%(uniq)s)''' % locals())

            outfile.write("%s" % fn)
            for slotname, label, pattern in slotnames:
                value = tuple(R('''ce$%s''' % slotname))
                if len(value) == 0:
                    value = ""
                outfile.write("\t%s" % pattern % value[0])
            outfile.write("\n")
            outfile.close()

    if "dmr" in options.toolset or "correlation" in options.toolset \
       or do_all:
        # build four sets
        for x, fn in enumerate(options.treatment_files):
            R('''treatment_R%(x)i = MEDIPS.createSet(
            file='%(fn)s',
            BSgenome='%(genome_file)s',
            shift=%(shift)i,
            extend=%(extend)i,
            window_size=%(window_size)i,
            uniq=%(uniq)s)''' % locals())
        R('''treatment_set = c(%s)''' %
          ",".join(["treatment_R%i" % x
                    for x in range(len(options.treatment_files))]))

        if options.control_files:
            for x, fn in enumerate(options.control_files):
                R('''control_R%(x)i = MEDIPS.createSet(
                file='%(fn)s',
                BSgenome='%(genome_file)s',
                shift=%(shift)i,
                extend=%(extend)i,
                window_size=%(window_size)i,
                uniq=%(uniq)s)''' % locals())
            R('''control_set = c(%s)''' %
              ",".join(["control_R%i" % x
                        for x in range(len(options.control_files))]))

        # build coupling vector
        R('''CS = MEDIPS.couplingVector(pattern="CG",
        refObj = treatment_set[[1]])''')

        if "correlation" in options.toolset or do_all:
            R('''cor.matrix = MEDIPS.correlation(
            c(treatment_set, control_set))''')

            R('''write.table(cor.matrix,
            file='%s',
            sep="\t")''' % E.getOutputFile("correlation"))

        if "dmr" in options.toolset or do_all:
            # Data that does not fit the model causes
            # "Error in 1:max_signal_index : argument of length 0"
            # The advice is to set MeDIP=FALSE
            # See: http://comments.gmane.org/gmane.science.biology.informatics.conductor/52319

            if options.is_medip:
                medip = "TRUE"
            else:
                medip = "FALSE"

            R('''meth = MEDIPS.meth(
            MSet1 = treatment_set,
            MSet2 = control_set,
            CSet = CS,
            ISet1 = NULL,
            ISet2 = NULL,
            p.adj = "bonferroni",
            diff.method = "edgeR",
            prob.method = "poisson",
            MeDIP = %(medip)s,
            CNV = F,
            type = "rpkm",
            minRowSum = 1)''' % locals())

            # test windows for differential methylation
            R('''tested = MEDIPS.selectSig(meth,
            adj=T,
            ratio=NULL,
            p.value=0.1,
            bg.counts=NULL,
            CNV=F)''')

            R('''write.table(tested,
            file=gzFile('%s', 'w')
            sep="\t",
            quote=F)''' % E.getOutputFile("windows"))

            # select gain and merge adjacent windows
            R('''gain = tested[which(tested[, grep("logFC", colnames(tested))] > 0),];
                 gain_merged = MEDIPS.mergeFrames(frames=gain, distance=1)''')

            R('''write.table(gain_merged,
            file=gzFile('%s', 'w')
            sep="\t",
            row.names=FALSE,
            col.names=FALSE)''' % E.getOutputFile("gain.bed.gz"))

            # select loss and merge adjacent windows
            R('''loss = tested[which(tested[, grep("logFC", colnames(tested))] < 0),];
                 loss_merged = MEDIPS.mergeFrames(frames=gain, distance=1)''')

            R('''write.table(loss_merged,
            file=gzFile('%s', 'w')
            sep="\t",
            row.names=FALSE,
            col.names=FALSE)''' % E.getOutputFile("loss.bed.gz"))

    # if "rpm" in options.toolset or do_all:
    #     outputfile = E.getOutputFile("rpm.wig")
    #     R('''MEDIPS.exportWIG(file = '%(outputfile)s',
    #     data = CONTROL.SET, raw = T, descr = "rpm")''' %
    #       locals())
    #     if options.bigwig:
    #         bigwig(outputfile, contig_sizes)
    #     else:
    #         compress(outputfile)

    # if "rms" in options.toolset or do_all:
    #     outputfile = E.getOutputFile("rms.wig")
    #     R('''MEDIPS.exportWIG(file = '%(outputfile)s',
    #     data = CONTROL.SET, raw = F, descr = "rms")''' %
    #       locals())
    #     if options.bigwig:
    #         bigwig(outputfile, contig_sizes)
    #     else:
    #         compress(outputfile)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
