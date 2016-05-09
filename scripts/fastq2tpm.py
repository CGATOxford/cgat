'''
fastq2tpm.py - use rapid/lightweight alignment RNA seq quantification methods
=============================================================================

:Author: MikeMorgan
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Wrapper for kallisto & sailfish, need to add in Salmon
NB: I wrote this before Tom implemented his mapper-class based
approach, so that should supersede this.  I just need to get round to
actually doing it.

Usage
-----

.. This might change, for now usage is limited to creating an index with either
Kallisto or Sailfish and quantifying from fastq files.

Example::

   python fastq2tpm.py

Type::

   python fastq2tpm.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGATPipelines.PipelineScRnaseqQc as scQC
import os
import subprocess


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help")

    parser.add_option("--program", dest="program", type="choice",
                      choices=["kallisto", "sailfish"],
                      help="use either kallisto or sailfish, "
                      "for alignment-free quantification")

    parser.add_option("--method", dest="method", type="choice",
                      choices=["make_index", "quant"],
                      help="method of kallisto to run")

    parser.add_option("--index-fasta", dest="fa_index", type="string",
                      help="multi-fasta to use to make index for kallisto")

    parser.add_option("--index-file", dest="index_file", type="string",
                      help="kallisto index file to use for quantificaiton")

    parser.add_option("--use-bias", dest="bias", action="store_true",
                      help="use kallisto's bias correction")

    parser.add_option("--bootstraps", dest="bootstrap", type="int",
                      help="number of bootstraps to apply to quantification")

    parser.add_option("--seed", dest="seed", type="int",
                      help="seed number for random number genration "
                      "and bootstrapping")

    parser.add_option("--just-text", dest="text_only", action="store_true",
                      help="only output files in plain text, not HDF5")

    parser.add_option("--library-type", dest="library", type="choice",
                      choices=["ISF", "ISR", "IU", "MSF", "MSR", "MU",
                               "OSF", "OSR", "OU", "SR", "SF", "U"],
                      help="sailfish fragment library type code")

    parser.add_option("--paired-end", dest="paired", action="store_true",
                      help="data are paired end")

    parser.add_option("--kmer-size", dest="kmer", type="int",
                      help="kmer size to use for index generation")

    parser.add_option("--gene-gtf", dest="gene_gtf", type="string",
                      help="GTF file containing transcripts and gene "
                      "identifiers to calculate gene-level estimates")

    parser.add_option("--threads", dest="threads", type="int",
                      help="number of threads to use for kallisto "
                      "quantificaion")

    parser.add_option("--output-directory", dest="outdir", type="string",
                      help="directory to output transcript abundance "
                      "estimates to")

    parser.add_option("--output-file", dest="outfile", type="string",
                      help="output filename")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.method == "make_index":
        if options.program == "kallisto":
            scQC.runKallistoIndex(fasta_file=options.fa_index,
                                  outfile=options.outfile,
                                  kmer=options.kmer)
        elif options.program == "sailfish":
            scQC.runSailfishIndex(fasta_file=options.fa_index,
                                  outdir=options.outdir,
                                  threads=options.threads,
                                  kmer=options.kmer)
        else:
            E.warn("program not recognised, exiting.")

    elif options.method == "quant":
        infiles = argv[-1]
        qfiles = infiles.split(",")
        # make the output directory if it doesn't exist
        if os.path.exists(options.outdir):
            pass
        else:
            os.system("mkdir %s" % options.outdir)

        if options.program == "kallisto":
            scQC.runKallistoQuant(fasta_index=options.index_file,
                                  fastq_files=qfiles,
                                  output_dir=options.outdir,
                                  bias=options.bias,
                                  bootstrap=options.bootstrap,
                                  seed=options.seed,
                                  threads=options.threads,
                                  plaintext=options.text_only)
        elif options.program == "sailfish":
            infiles = argv[-1]
            qfiles = infiles.split(",")
            scQC.runSailfishQuant(fasta_index=options.index_file,
                                  fastq_files=qfiles,
                                  output_dir=options.outdir,
                                  paired=options.paired,
                                  library=options.library,
                                  threads=options.threads,
                                  gene_gtf=options.gene_gtf)

        else:
            E.warn("program not recognised, exiting.")
    else:
        pass

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
