'''
assoc2plot.py - script for post-GWAS plotting
====================================================

:Author:
:Tags: Python

Purpose
-------

.. Post-gwas plotting, for visualization and QC

Usage
-----

.. Example use case

Example::

   python assoc2plot.py

Type::

   python assoc2plot.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.PipelineGWAS as gwas


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

    parser.add_option("--plot-type", dest="plot_type", type="choice",
                      choices=["manhattan", "qqplot", "epistasis"],
                      help="plot type to generate")

    parser.add_option("--resolution", dest="resolution", type="choice",
                      choices=["genome_wide", "chromosome",
                               "fine_map"],
                      help="the resolution of plotting, wether the plot "
                      "depicts the whole genome, a single chromosome or "
                      "a specific locus")

    parser.add_option("--save-path", dest="save_path", type="string",
                      help="path and filename to save image to")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    parser.set_defaults(resolution="genome_wide",
                        plot_type="manhattan")

    # if the input is a list of files, split them
    infile = argv[-1]
    infiles = infile.split(",")
    
    # need to parse epistasis output slightly differently
    if options.plot_type == "epistasis":
        epi = True
    else:
        epi = False

    if len(infiles) > 1:
        results = gwas.GWASResults(assoc_file=infiles,
                                   epistasis=epi)
    elif len(infiles) == 1:
        results = gwas.GWASResults(assoc_file=infile,
                                   epistasis=epi)
    else:
        raise IOError("no input files detected, please specifiy association "
                      "results files as the last command line argument")

    if options.plot_type == "manhattan":
        df = results.plotManhattan(resolution=options.resolution,
                                   save_path=options.save_path)
    elif options.plot_type == "qqplot":
        results.plotQQ(save_path=options.save_path,
                       resolution=options.resolution)
    elif options.plot_type == "epistasis":
        results.plotEpistasis(save_path=options.save_path,
                              resolution=options.resolution)
    else:
        pass

    # only output appended results for Manhattan plot, not qqplot
    try:
        df.to_csv(options.stdout, sep="\t", index=None)
    except UnboundLocalError:
        pass

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
