'''
extract_stats.py - extract and process tables from CSVDB
========================================================

:Author: Mike Morgan
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Extract tables from sqlite databases and process

Usage
-----

.. Example use case

Example::

   python extract_stats.py

Type::

   python extract_stats.py --help

for command line help.

Command line options
--------------------

tasks
+++++

`extract_table` - extract a table from a database
                  containing relevant information

`get_coverage` - calculate gene/transcript model
                 coverage stats - only works on single cell data
                 with filename format <seqRun>_<plate>_<well>_<mappper>

`aggregate` - aggregate together multiple stats tables,
              and select relevant measures

'''

import sys
import CGAT.Experiment as E
import pandas as pd
import pandas.io.sql as pdsql
import numpy as np
import re
import sqlite3 as sql


def getTableFromDb(db, table):
    '''
    Get a table from a database with pandas
    '''

    state = ''' SELECT * FROM %(table)s;''' % locals()
    dbh = sql.connect(db)
    df = pdsql.read_sql(state, dbh)
    df.index = df["track"]
    df.drop(labels="track", inplace=True,
            axis=1)

    return df


def cleanStatsTable(stats_file):
    '''
    Take in a table containing aggregated stats
    and clean by removing duplicate columns
    '''

    _df = pd.read_table(stats_file, sep="\t", header=0,
                        index_col=None, mangle_dupe_cols=False)
    # drop duplicates is case sensitive, convert all to
    # same case - SQL is not case sensitive so will throw
    # a hissy fit for same column names in different cases
    _df.columns = [cx.lower() for cx in _df.columns]
    _df = _df.T.drop_duplicates().T
    _df.index = _df["track"]
    return _df


def extractTranscriptCounts(con, table):
    '''
    Extract transcript model counts for a
    given sample

    Arguments
    ---------
    con: sqlite.connection
      An SQLite connection

    table: string
      the table to extract the transcript counts
      from.

    Returns
    -------
    coverages: pandas.Core.Series
    '''

    statement = '''
    SELECT coverage_sense_pcovered
    FROM %(table)s
    WHERE coverage_sense_nval > 0;
    ''' % locals()

    coverages = pdsql.read_sql(statement, con)
    coverages = coverages.loc[:, "coverage_sense_pcovered"]
    return coverages


def summariseOverBins(coverages, bins):
    '''
    Summarise model coverages over a set of bins

    Argumnets
    ---------
    coverages: pandas.Core.Series
      coverages over gene/transcripts

    bins: list
      values corresponding to percentage bins

    Returns
    -------
    freqs: numpy.array
      frequency array of coverages over percentiles
    '''

    freqs = np.zeros(shape=len(bins),
                     dtype=np.float64)
    for i in range(len(bins)):
        if i == 0:
            hits = coverages <= bins[i]
        else:
            hits = (coverages <= bins[i]) & (coverages > bins[i-1])

        freqs[i] = len(coverages[hits])

    return freqs


def getModelCoverage(db, table_regex, model_type="transcript"):
    '''
    Compute transcript model coverage stats

    Arguments
    ---------
    db: string
      database containing transcript counts

    table_regex: string
      regular expression for transcript count table

    model_type: string
      calculate coverages over either transcripts or
      genes.  Default is gene models

    Returns
    -------
    coverage_df: Pandas.Core.DataFrame
      model coverage stats summarised for each cell
    '''

    # need to regex for all the tables, one for each sample
    # fetch_all returns a list of tuples
    dbh = sql.connect(db)
    cursor = dbh.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")

    tab_reg = re.compile(table_regex)
    table_list = [tx[0] for tx in cursor.fetchall() if re.search(tab_reg,
                                                                 tx[0])]

    # pull out counts for each cell and compute coverages

    bins = range(0, 101)
    cov_dict = {}
    for tab in table_list:
        covs = extractTranscriptCounts(dbh, tab)
        freq_array = summariseOverBins(covs, bins)
        cov_dict[tab] = freq_array

    coverage_df = pd.DataFrame(cov_dict).T
    # create a regex group to remove superfluous characters
    # from the track names
    ix_re = re.compile("_(?P<run>\d+)_(?P<plate>\d+)_(?P<well>\d+)_(?P<mapper>\S+)_transcript_counts")
    re_matches = [re.match(ix_re, ix) for ix in coverage_df.index]
    indx = ["%s_%s-%s.%s" % rm.group(1, 2, 3, 4) for rm in re_matches]
    coverage_df.index = indx
    coverage_df.columns = ["Bin%i" % bx for bx in coverage_df.columns]
    return coverage_df


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
        out_df = getTableFromDb(db=options.database,
                                table=options.table)

    elif options.task == "get_coverage":
        out_df = getModelCoverage(db=options.database,
                                  table_regex="(\S+)_transcript_counts")

    elif options.task == "clean_table":
        infile = argv[-1]
        out_df = cleanStatsTable(infile)

    out_df.to_csv(options.stdout,
                  sep="\t", index_label="track")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
