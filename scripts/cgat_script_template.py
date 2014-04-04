'''
cgat_script_template.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

<Overall purpose and function of the script>

Options
-------

<Options for the script, with detail of how they are combined to provide \
intended functionality>

Usage
-----
<Example use case>

Example::

   python cgat_script_template.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import pysam
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.pipeline as P


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--input-file", dest="infile", type="string",
                      help="input bam file")

    parser.add_option("--genome-dir", dest="gDir", type="string",
                      help = "genome directory containing fasta format files")
    
    parser.add_option("--genome", dest"genome", type="string",
                      help = "the genome to which reads were aligned")

    parser.add_option("--contig-file", dest="contig", type="string",
                      help = "list of contigs and their sizes in tsv format")

    parset.set_defaults(infile = None, 
                        gDir = None,
                        genome = None, 
                        contig = None)
    
    if len(args) == 0:
        args.append("-")
        
    # currently assumes that the input file is bam format
    # will need to change to either accept sam format
    # or automatically detect format

    samfile = pysam.Samfile(options.infile, "rb")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    def contigFileMaker(genome_dir, genome):
        '''Make a contig file unless a file is passed through options.contig.
        This code is ripped from pipeline_annotations function buildContigSizes'''

        genome_file = open("%s/%s.fasta" % (genome_dir, genome), "r")
        prefix = "%s/%s" % (genome_dir, genome)
        fasta  = IndexedFasta.InfexedFasta(genome_dir)
        
        with open("contig.tsv", "w") as outs:
            
            for contig, size in fasta.getContigSizes(with_synonyms = False).iteritems():
                outs.write("%s\t%i\n" % (contig, size))
        
        genome_file.close

    if options.contig == None:
        contigFileMaker(options.gDir, options.genome)
    else: continue

    contig_dict = {}


    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
