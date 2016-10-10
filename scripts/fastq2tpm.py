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
import os
import subprocess
import re
import numpy as np


# ------------------------------------------------------- #
# Functions for expression quantification
# ------------------------------------------------------- #

def runSailfishIndex(fasta_file, outdir, threads,
                     kmer):
    '''
    Wrapper for sailfish index
    '''

    if fasta_file.endswith(".fa"):
        pass
    elif fasta_file.endswith(".fasta"):
        pass
    else:
        E.warn("are you sure this is a fasta file?")

    command = '''
    sailfish index --transcripts %s --out %s --threads %i --kmerSize %i
    ''' % (fasta_file, outdir, threads, kmer)

    os.system(command)


def runSailfishQuant(fasta_index, fastq_files, output_dir,
                     paired=False, library="ISF", threads=4,
                     gene_gtf=None):
    '''
    Wrapper for sailfish quant command
    '''

    decompress = False
    if len(fastq_files) > 1:
        if fastq_files[0].endswith(".gz"):
            decompress = True
        else:
            pass
    else:
        if fastq_files[0].endswith(".gz"):
            decompress = True
        else:
            pass

    # check output directory is an absolute path
    if os.path.isabs(output_dir):
        pass
    else:
        out_dir = os.path.abspath(output_dir)

    states = []
    command = " sailfish quant --index %s -l %s  -o %s " % (fasta_index,
                                                            library,
                                                            output_dir)

    states.append(command)

    if threads:
        states.append(" --threads %i " % threads)
    else:
        pass

    if gene_gtf:
        states.append(" --geneMap %s " % gene_gtf)
    else:
        pass

    # sailfish does not handle compress files natively,
    # need to decompress on the fly with advanced
    # bash syntax
    if decompress and paired:
        first_mates = tuple([fq for fq in fastq_files if re.search("fastq.1.gz",
                                                                   fq)])
        fstr_format = " ".join(["%s" for hq in first_mates])
        fdecomp_format = fstr_format % first_mates
        decomp_first = " -1 <(zcat %s)" % fdecomp_format

        states.append(decomp_first)

        second_mates = tuple([sq for sq in fastq_files if re.search("fastq.2.gz",
                                                                    sq)])
        sstr_format = " ".join(["%s" for aq in second_mates])
        sdecomp_format = sstr_format % second_mates
        decomp_second = " -2 <(zcat %s)" % sdecomp_format

        states.append(decomp_second)

    elif decompress and not paired:
        first_mates = tuple([fq for fq in fastq_files if re.search("fastq.gz",
                                                                   fq)])
        fstr_format = " ".join(["%s" for sq in first_mates])
        fdecomp_format = fstr_format % first_mates
        decomp_first = " -r <(zcat %s)" % fdecomp_format

        states.append(decomp_first)

    elif paired and not decompress:
        first_mates = tuple([fq for fq in fastq_files if re.search("fastq.1",
                                                                   fq)])
        fstr_format = " ".join(["%s" for sq in first_mates])
        fdecomp_format = fstr_format % first_mates
        decomp_first = " -1 %s " % fdecomp_format

        states.append(decomp_first)

        second_mates = tuple([sq for sq in fastq_files if re.search("fastq.2",
                                                                    sq)])
        sstr_format = " ".join(["%s" for aq in second_mates])
        sdecomp_format = sstr_format % second_mates
        decomp_second = " -2 %s " % sdecomp_format

        states.append(decomp_second)

    statement = " ".join(states)

    # subprocess cannot handle process substitution
    # therefore needs to be wrapped in /bin/bash -c '...'
    # for bash to interpret the substitution correctly

    process = subprocess.Popen(statement, shell=True,
                               executable="/bin/bash")

    stdout, stderr = process.communicate()

    if process.returncode != 0:
        raise OSError(
            "-------------------------------------------\n"
            "Child was terminated by signal %i: \n"
            "The stderr was \n%s\n%s\n"
            "-------------------------------------------" %
            (-process.returncode, stderr, statement))


def runKallistoIndex(fasta_file, outfile, kmer=31):
    '''
    Wrapper for kallisto index
    '''

    if fasta_file.endswith(".fa"):
        pass
    elif fast_file.endswith(".fasta"):
        pass
    else:
        E.warn("are you sure this is a fasta file?")

    command = "kallisto index --index=%s  %s" % (outfile,
                                                 fasta_file)

    os.system(command)


def runKallistoQuant(fasta_index, fastq_files, output_dir,
                     bias=False, bootstrap=None,
                     seed=1245, threads=None, plaintext=False):
    '''
    Wrapper for kallisto quant command
    '''

    if len(fastq_files) > 1:
        fastqs = " ".join(fastq_files)
    else:
        fastqs = fastq_files

    # check output directory is an absolute path
    if os.path.isabs(output_dir):
        pass
    else:
        out_dir = os.path.abspath(output_dir)

    states = []
    command = " kallisto quant --index=%s --output-dir=%s" % (fasta_index,
                                                              output_dir)
    states.append(command)

    if bias:
        states.append(" --use-bias ")
    else:
        pass

    if bootstrap:
        states.append(" --bootstrap=%i --seed=%i " % (bootstrap,
                                                      seed))
    else:
        pass

    if plaintext:
        states.append(" --plaintext ")
    else:
        pass

    if threads:
        states.append(" --threads=%i " % threads)
    else:
        pass

    states.append(" %s " % fastqs)

    statement = " ".join(states)

    # need to rename output files to conform to input/output
    # pattern as required.  Default name is abundance*.txt
    # when using plaintext output
    # kaliisto requires an output directory - create many small
    # directories, one for each file.
    # then extract the abundance.txt file and rename using the
    # input/output pattern

    os.system(statement)


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

    parser.set_defaults(paired=False)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.method == "make_index":
        if options.program == "kallisto":
            runKallistoIndex(fasta_file=options.fa_index,
                             outfile=options.outfile,
                             kmer=options.kmer)
        elif options.program == "sailfish":
            runSailfishIndex(fasta_file=options.fa_index,
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
            runKallistoQuant(fasta_index=options.index_file,
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
            runSailfishQuant(fasta_index=options.index_file,
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
