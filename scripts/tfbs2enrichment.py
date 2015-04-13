'''
tfbs2enrichment.py - template for CGAT scripts
====================================================

:Author: Mike Morgan
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Runs transcription factor motif enrichment using the
Transfac(R) database.  Signficance testing is currently
carried out by Fisher's Exact test based on a background
gene set matched for %CG content.

Currently only for use in pipeline_transfacmatch.py, this is
not a stand-alone script.

Usage
-----

.. Example use case

Example::

   python tfbs2enrichment.py

Type::

   python tfbs2enrichment.py --help

for command line help.

Command line options
--------------------

'''

import sys
import sqlite3
import collections
import pandas.io.sql as pdsql
import numpy as np
from rpy2.robjects import r as R
from rpy2.robjects.numpy2ri import numpy2ri
import CGAT.Experiment as E
import CGAT.IOTools as IOTools


def testSignificanceOfMatrices(background,
                               foreground,
                               database,
                               match_table,
                               outfile,
                               header_line=True,
                               direction="two.sided"):
    '''
    Uses the fishers exact test to estimate the significance of
    transcription factor motifs in a foregound set of intervals
    compared to a background set
    '''
    # jj: This is crashing because of global namespace issues
    R('''rm(list = ls())''')

    # get background and foreground sets
    interval_sets = {"foreground": set(), "background": set()}
    for name, fn in (("background", background),
                     ("foreground", foreground)):
        inf = IOTools.openFile(fn)
        if header_line:
            inf.readline()
        for line in inf:
            interval_sets[name].add(line[:-1])

    # connect to database
    dbh = sqlite3.connect(database)

    # ids_tfs dict - key <sequence_id>: val <set of bound tf_ids>
    # all_tfs - set containing all the tf ids
    # NB. each tf id is only counted once for each sequence id,
    # regardless of no. times motif occurs in sequence

    ids_tfs = collections.defaultdict(set)
    all_tfs = set()

    # MM: this causes a memory error > 24G used per job!!
    # read in table in chunks, iterate with a generator
    # reduces memory usage < 1G
    # may increase processing time as millions of lines

    # for data in cc.execute("SELECT seq_id, matrix_id FROM %s" %
    #                        match_table).fetchall():
    #     seq_id, tf = data[0], data[1]
    #     ids_tfs[seq_id].add(tf)
    #     all_tfs.add(tf)

    state = "SELECT seq_id, matrix_id FROM %s" % match_table
    sql_df = pdsql.read_sql(sql=state, con=dbh, index_col=None,
                            chunksize=1000000)

    # make chunksize bigger or user defined?
    # with 1000000, takes < 2s per iteration
    for chunk in sql_df:
        _df = chunk
        sql_group = _df.groupby(by='seq_id')

        for names, groups in sql_group:
            tfs = groups[groups["seq_id"] == names]['matrix_id'].tolist()
            ids_tfs[names].update(tfs)
            all_tfs.update(tfs)

    outf = IOTools.openFile(outfile, "w")
    outf.write("matrix_id\toddsRatio\t95CI_low\t95_CI_hi\tpvalue\t"
               "nforeground\tnbackground\ttforeground\ttbackground\tqvalue\n")

    # iterate over factors and calculate significance using fishers exact test
    fisherpy = R["fisher.test"]
    padjpy = R["p.adjust"]

    results = []
    pvalues = []
    # perform hypergeometric test for each transcription factor in turn
    for tf in all_tfs:
        # initialise contingency table...
        contingency = np.zeros((2, 2))
        # populate contingency table...
        for seq_id in ids_tfs.keys():
            if seq_id in interval_sets["foreground"]:
                if tf in ids_tfs[seq_id]:
                    contingency[0, 0] += 1
                elif tf not in ids_tfs[seq_id]:
                    contingency[1, 0] += 1
            elif seq_id in interval_sets["background"]:
                if tf in ids_tfs[seq_id]:
                    contingency[0, 1] += 1
                elif tf not in ids_tfs[seq_id]:
                    contingency[1, 1] += 1
        # # Debug:
        # print tf
        # print contingency

        # run fishers exact test in R
        f = fisherpy(numpy2ri(contingency), alternative=direction)

        # get overlap numbers
        nforeground = contingency[0, 0]
        nbackground = contingency[0, 1]
        tforeground = contingency[0, 0] + contingency[1, 0]
        tbackground = contingency[0, 1] + contingency[1, 1]

        # convert back to python object...
        f = [list(x) for x in np.array(f)]

        # fisher.test returns pval, conf intervals, estimate of odds ratio
        pvalue, ci_low, ci_hi, OR = f[0][0], f[1][0], f[1][1], f[2][0]
        pvalues.append(pvalue)

        # if every background and every foreground, or no background
        # and no foreground interals are hit, make OR = 1 by definition
        if (nforeground + nbackground == 0 or
                (nforeground == tforeground and nbackground == tbackground)):
            OR = 1

        results.append([tf, OR, ci_low, ci_hi, pvalue, nforeground,
                        nbackground, tforeground, tbackground])

    # correct for multiple comparisons
    qvalues = padjpy(pvalues)
    for i in range(len(results)):
        outf.write("\t".join(map(str, results[i] + [qvalues[i]])) + "\n")
    outf.close()


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--foreground", dest="foreground", type="string",
                      help="foreground file")

    parser.add_option("--background", dest="background", type="string",
                      help="matched background file")

    parser.add_option("--database", dest="database", type="string",
                      help="PATH to sqlite database containing Transfac Match"
                      " output.")

    parser.add_option("--match-table", dest="match_table", type="string",
                      help="tablename containing Transfac Match output")

    parser.add_option("--outfile", dest="outfile", type="string",
                      help="output file name")

    parser.add_option("--direction", dest="direction", type="string",
                      help="direction to test for significance in Fisher's "
                      "exact test. Default=2-tailed")

    parser.add_option("--geneset-header", dest="genesets", type="string",
                      help="geneset header")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    PipelineTFM.testSignificanceOfMatrices(options.background,
                                           options.foreground,
                                           options.database,
                                           options.match_table,
                                           options.outfile,
                                           options.genesets,
                                           options.direction)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
