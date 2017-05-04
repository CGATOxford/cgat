'''data2resamples.py - template for CGAT scripts
====================================================

:Tags: Python

Purpose
-------

Generate ``n`` resampled data sets.  Randomly resampling a single
replicate with replacement at each time point.

Usage
-----

To generate n number of resampled pseudo-data sets.  For each time
point a random replicate is selected and the expression profiles for
all genes generated as a single file.  File names are output with an
iteration number for identification in later stages.  This is akin to
bootstrapping the samples rather than the genes.  Data can later be
combined by either averaging across data-sets or using consensus
clustering as implemented in pipeline_timeseries.py.

File names are based on the condition, repicates, resample number and
input gtf file.

Options
-------

All options are mandatory.

  --timepoints - a comma separated list of time points that have been measured

  --replicates - labels used to identify replicates in the input file

  --condition - if appicable a label to apply to the output file.  Especially
                useful if there are multiple conditions in your data.

  --resamples - the number of pseudo-data sets to generate.  The recommended
                number of times to resample the data is dependent on the
                application.  For clustering 100 resamples would be recommended
                if consensus clustering is used.

  --input-gtf - the gtf file used to generate counts/expression dat

  --seed - seed for pseudo-random number generation

Example::

  head IFNgamma-refcoding-expresion.tsv

  gene_id             1      2      3      4      5       6       7
  ENSMUSG00000001305  10.27  10.29  13.19  13.46  11.93   10.62   11.08
  ENSMUSG00000001674  11.9   11.76  14.38  14.32  12.88   12.20   12.27
  ENSMUSG00000003134  6.67   5.93   4.66   5.7    7.47    9.08    8.36
  ENSMUSG00000004110  13.41  12.98  13.10  8.64   7.77    9.12    9.13
  ENSMUSG00000004451  9.24   9.88   9.75   9.55   11.32   11.94   11.64
  ...
  times               0      0      0      1      1       1       3
  replicates          R1     R2     R3     R1     R2      R3      R1

   python data2resamples.py
          --time=0,1,3,6,12,24,48
          --replicates=R1,R2,R3
          --condition=IFNgamma
          --resamples=100
          --input-gtf=ensembl-refcoding74.gtf.gz
          < IFNgamma-refcoding-expression.tsv

   head IFN-refcoding-expression-resample_23-expression.tsv

  gene_id             0      1      3      6      12      24      48
  ENSMUSG00000001305  10.27  13.46  12.84  11.46  10.93   8.78    7.77
  ENSMUSG00000001674  11.76  12.20  13.55  14.32  14.88   14.28   13.84
  ENSMUSG00000003134  6.67   5.70   8.36   7.14   7.11    7.83    8.06
  ENSMUSG00000004110  12.98  8.64   9.13   7.40   7.09    6.01    6.18
  ENSMUSG00000004451  9.88   11.32  11.77  11.89  12.44   14.26   14.87

Type::

   python data2resamples.py --help

for command line help.

Command line options
--------------------

'''

import sys
import pandas as pd
import CGAT.Experiment as E
import CGAT.Timeseries as TS
import CGAT.IOTools as IOTools


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

    parser.add_option("--time", dest="timepoints", type="string",
                      help="a comma-separated list of time points measured")

    parser.add_option("--replicates", dest="reps", type="string",
                      help="a comma-separated list of replicate IDs")

    parser.add_option("--condition", dest="condition", type="string",
                      help="experimental condition")

    parser.add_option("--resamples", dest="resamples", type="string",
                      help="number of times to resample replicates to"
                      " generate pseudo datasets")

    parser.add_option("--input-gtf", dest="gtf_file", type="string",
                      help="reference gtf file")

    parser.add_option("--output-file-directory", dest="output_dir",
                      type="string", help="directory to output"
                      " resampled files to")

    # add common options (-h/--help, ...) and parse command line

    (options, args) = E.Start(parser, argv=argv)

    try:
        infile = IOTools.openFile(argv[-1], "r")
    except IOError:
        infile = options.stdin

    data_frame = pd.read_table(infile,
                               sep="\t",
                               index_col=0,
                               header=0)
    time_str = options.timepoints.split(",")
    time_points = [int(x) for x in time_str]
    replicates = options.reps.split(",")
    reps = int(options.resamples)

    its = [time_str, replicates]
    midx = pd.MultiIndex.from_product(its,
                                      names=['times', 'replicates'])

    TS.genResampleData(data_frame=data_frame,
                       multiple_index=midx,
                       replicates=reps,
                       sample_reps=replicates,
                       times=time_points,
                       condition=options.condition,
                       ref_gtf=options.gtf_file,
                       out_dir=options.output_dir,
                       seed=int(options.random_seed))

    # Write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
