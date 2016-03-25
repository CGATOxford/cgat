'''
extract_stats.py - template for CGAT scripts
====================================================

:Author: Mike Morgan
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Extract single cell metrics from sqlite databases

Usage
-----

.. Example use case

Example::

   python cgat_script_template.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

tasks
+++++

`extract_table` - extract a table from a database
                  containing relevant information

`get_coverage` - calculate gene/transcript model
                 coverage stats

`aggregate` - aggregate together multiple stats tables,
              and select relevant measures


'''

import sys
import CGAT.Experiment as E
import CGATPipelines.PipelineScRnaseqQc as scQC


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--task", dest="task", type="choice",
                      choices=["extract_table", "get_coverage",
                               "clean_table"],
                      help="task to perform")

    parser.add_option("-d", "--database", dest="database", type="string",
                      help="SQLite database containing relevant tables")

    parser.add_option("-t", "--table-name", dest="table", type="string",
                      help="table in SQLite DB to extract")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.task == "extract_table":
        out_df = scQC.getTableFromDb(db=options.database,
                                     table=options.table)

    elif options.task == "get_coverage":
        out_df = scQC.getModelCoverage(db=options.database,
                                       table_regex="(\S+)_transcript_counts")

    elif options.task == "clean_table":
        infile = argv[-1]
        out_df = scQC.cleanStatsTable(infile)

    out_df.to_csv(options.stdout,
                  sep="\t", index_label="track")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
