'''
pheno2plot.py - format and manipulate phenotype files
====================================================

:Author:
:Tags: Python

Purpose
-------

.. Generate plots of phenotype data

Usage
-----

.. Example use case

Example::

   python pheno2plot.py

Type::

   python pheno2plot.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
from rpy2.robjects import r as R
from rpy2.robjects import pandas2ri
import pandas as pd
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

    parser.add_option("-p", "--plot-type", dest="plot_type", type="choice",
                      choices=["histogram", "barplot", "density",
                               "boxplot", "scatter", "map",
                               "pca"],
                      help="the plot type to generate")

    parser.add_option("--plot-n-pc", dest="n_pcs", type="int",
                      help="The number of principal components to "
                      "plot")

    parser.add_option("-g", "--group-by", dest="group_by", type="string",
                      help="column header to group observations by")

    parser.add_option("-x", "--x-column", dest="x_col", type="string",
                      help="column to plot on X axis")

    parser.add_option("-y", "--y-column", dest="y_col", type="string",
                      help="column to plot on y axis")

    parser.add_option("-i", "--index_column", dest="indx", type="string",
                      help="column number that refers to the row index")

    parser.add_option("--output-file", dest="outfile", type="string",
                      help="path and filename to save plot to")

    parser.add_option("--labels", dest="labels", type="string",
                      help="a comma-separated list of axis labels. "
                      "The first 2 correspond to the X and Y-axis, "
                      "respectively, and the third is the plot title")

    parser.add_option("--metadata-file", dest="meta_file", type="string",
                      help="file containing metadata for annotating "
                      "plots with. Use `--group-labels` to define table "
                      "columns to use")

    parser.add_option("--fam-file", dest="fam_file", type="string",
                      help="Plink .fam file containing file IDs")

    parser.add_option("--xvar-labels", dest="xvar_labs", type="string",
                      help="a comma-separated list of variable X labels"
                      "only applies when X is a discrete or categorical "
                      "variable. The labels must be in the correct order")

    parser.add_option("--group-labels", dest="group_labs", type="string",
                      help="a comma-separated list of grouping variable "
                      "labels.  Can only apply when the grouping variable "
                      "is discrete or categorical.  The labels must be "
                      "input in the order of the data")

    parser.add_option("--yvar-labels", dest="yvar_labs", type="string",
                      help="a comma-separated list of variable Y labels"
                      "only applies when Y is a discrete or categorical "
                      "variable")

    parser.add_option("--var-type", dest="var_type", type="choice",
                      choices=["continuous", "categorical", "integer"],
                      help="The data type of the variables to be plotted."
                      "The default is continuous")

    parser.add_option("--coordinate-file", dest="coordinates", type="string",
                      help="file containing co-ordinates data")

    parser.add_option("--coords-id-col", dest="coord_ids", type="string",
                      help="column header containing individual IDs")

    parser.add_option("--lattitude-column", dest="lat_col", type="string",
                      help="column header containing lattitude co-ordinates")

    parser.add_option("--longitude-column", dest="long_col", type="string",
                      help="column header containing longitude co-ordinates")

    parser.add_option("--reference-value", dest="ref_val", type="string",
                      help="categorical variable level to dichotomise on")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    parser.set_defaults(y_col=None,
                        group_by=None,
                        indx=None,
                        labels="X,Y,title",
                        xvar_labs=None,
                        yvar_labs=None,
                        var_type="continuous")
    infile = argv[-1]

    df = pd.read_table(infile, sep="\t", index_col=options.indx,
                       header=0)

    if options.plot_type == "map":
        df = pd.read_table(infile, sep="\t", index_col=options.indx,
                           header=0)

        coords_df = pd.read_table(options.coordinates, sep="\t",
                                  header=0, index_col=options.indx)
        gwas.plotMapPhenotype(data=df,
                              coords=coords_df,
                              coord_id_col=options.coord_ids,
                              lat_col=options.lat_col,
                              long_col=options.long_col,
                              save_path=options.outfile,
                              xvar=options.x_col,
                              var_type=options.var_type,
                              xlabels=options.xvar_labs,
                              level=options.ref_val)

    elif options.plot_type == "pca":
        data = gwas.parseFlashPCA(pcs_file=infile,
                                  fam_file=options.fam_file)

        gwas.plotPCA(data=data,
                     nPCs=options.n_pcs,
                     point_labels=options.group_labs,
                     save_path=options.outfile,
                     headers=False,
                     metadata=options.meta_file,
                     multiplot=True)
    else:
        df = pd.read_table(infile, sep="\t", index_col=options.indx,
                           header=0)

        gwas.plotPhenotype(data=df,
                           plot_type=options.plot_type,
                           x=options.x_col,
                           y=options.y_col,
                           group=options.group_by,
                           save_path=options.outfile,
                           labels=options.labels,
                           xlabels=options.xvar_labs,
                           ylabels=options.yvar_labs,
                           glabels=options.group_labs,
                           var_type=options.var_type)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
