'''
PipelineTransfacMatch.py
========================
 - Classes and functions used in pipeline_trasfacmatch.py
=========================================================
'''

################
# import modules
################
import re
import sys
import os
import CGAT.IOTools as IOTools
import sqlite3
import collections
import pandas as pd
import pandas.io.sql as pdsql
import CGAT.Pipeline as P
import numpy as np
from rpy2.robjects import r as R
from rpy2.robjects.packages import importr
from rpy2.robjects.numpy2ri import numpy2ri
import rpy2.robjects as robjects
import glob
import random
import pickle
import CGAT.FastaIterator as FastaIterator
import CGAT.Experiment as E
from math import factorial

########################################
# helper function
########################################


def filenameToTablename(filename):
    '''
    converts filename containing "." to tablename where "." converted to "_"
    '''
    return filename.replace(".", "_")


####################################
# Dealing with Match output
####################################
class Match:

    '''
    class with methods to parse the output from match.
    Each instance of this class will be an identifiable
    sequence
    '''

    def __init__(self):

        seq_id = None
        matrix_ids = None
        position = None
        strand = None
        core_score = 0
        matrix_score = 0
        sequence = None

    def parseInfo(self, infile):
        inf = open(infile)
        matrix_library = inf.readline()
        sequence_file = inf.readline()
        profile = inf.readline()
        inf.readline()
        inf.readline()

        # initialise with the first sequence identifier
        self.seq_id = inf.readline()[:-1].split(" ")
        self.seq_id = self.seq_id[len(self.seq_id) - 1:][0]

        # iterate over file and yield information for
        # each sequence identifier
        for line in inf.readlines():
            if line.find("Inspecting sequence ID") != -1:
                self.seq_id = line[:-1].split(" ")
                self.seq_id = self.seq_id[len(self.seq_id) - 1:][0]
            else:
                if line.strip() == "":
                    continue
                elif line.find("Total") != -1 or line.find("Frequency") != -1:
                    continue
                self.matrix_id, position_strand, self.core_score, \
                    self.matrix_score, self.sequence = [
                        x.strip() for x in line[:-1].split("|")]
                self.position = position_strand.split(" ")[0]
                self.strand = position_strand.split(
                    " ")[1].split("(")[1].split(")")[0]
                match = Match()
                yield self


def match_iterator(infile):
    '''
    return a generator object of class Match
    '''
    match = Match()
    return match.parseInfo(infile)

####################################
# Match metrics
####################################


def frequencyMetrics(database, match_table, outfile):
    '''
    builds metrics using Match data loaded into an sqlite database
    '''
    tf_number = collections.defaultdict(set)
    tf_max = collections.defaultdict(list)

    # connect to database
    dbh = sqlite3.connect(database)
    cc = dbh.cursor()

    # open outfile
    outf = open(outfile, "w")
    outf.write("seq_id\ttotal\tmax\n")

    # get the total number of tfs that have motifs per sequence
    for data in cc.execute("""SELECT seq_id, matrix_id FROM %s""" %
                           match_table).fetchall():
        seq_id, tfid = data[0], data[1]
        tf_number[seq_id].add(tfid)
        tf_max[seq_id].append(tfid)

    # get max number of hits for any transcription factor
    for seq_id, tfs in tf_max.iteritems():
        result = collections.defaultdict(int)
        for tf in tfs:
            result[tf] += 1

        # output the results
        outf.write("\t".join(map(str, [seq_id,
                                       len(tf_number[seq_id]),
                                       max(result.values())])) + "\n")
    outf.close()

###################################################
# Significance testing
###################################################


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
    background_file = open(background)
    foreground_file = open(foreground)
    interval_sets = {"foreground": set(), "background": set()}
    for inf_name in ["background_file", "foreground_file"]:
        name = P.snip(inf_name, "_file")
        if header_line:
            eval(inf_name).readline()
        for line in eval(inf_name).readlines():
            interval_sets[name].add(line[:-1])

    # connect to database
    dbh = sqlite3.connect(database)
    cc = dbh.cursor()

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

    outf = open(outfile, "w")
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
        temp = P.getTempFile("/ifs/scratch")
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

######################################################
# CpG matching
######################################################


def calculateSequenceComposition(interval_names,
                                 sequence_file,
                                 outfile,
                                 header_line=True):
    '''
    given a set of interval names that are present in a
    fasta file, return CpG content file
    '''
    interval_file = open(interval_names)
    if header_line:
        interval_file.readline()
    sequence_file = open(sequence_file)

    interval_set = set()
    for line in interval_file.readlines():
        interval_set.add(line[:-1])

    temp = P.getTempFile("/ifs/scratch")
    for record in FastaIterator.iterate(sequence_file):
        seq_id = record.title.split(" ")[0]
        if seq_id in interval_set:
            temp.write(">%s\n%s\n" % (record.title, record.sequence))
    temp.close()

    inf = temp.name
    statement = '''
    cat %(inf)s | python %(scriptsdir)s/fasta2table.py
    -s na -s cpg -s length
    --log=%(outfile)s.log > %(outfile)s'''
    P.run()

####################################################
####################################################


def matchBgSequenceComposition(gc_load_files,
                               background,
                               foreground,
                               fasta_file,
                               outfile,
                               database="csvdb",
                               header_line=True,
                               bg_stat="pCpG",
                               stat="fisher"):
    '''
    take the background set and subset it for intervals with a sequence
    composition distribution that is the same as the foreground set.
    Subsetting is done without replacement.
    This requires that the background set is sufficiently large, if the
    returned matched_background set is <90% size of foreground set, the
    pipeline will crash.
    '''
    background_file = open(background)
    foreground_file = open(foreground)

    if header_line:
        background_file.readline()
        foreground_file.readline()

    background_set = set()
    foreground_set = set()
    for interval_id in background_file.readlines():
        background_set.add(interval_id[:-1])
    for interval_id in foreground_file.readlines():
        foreground_set.add(interval_id[:-1])

    dbh = sqlite3.connect(database)
    cc = dbh.cursor()
    tablenames = [filenameToTablename(
        P.snip(os.path.basename(x), ".load")) for x in gc_load_files]

    # jj: cpg scores rounded to three dp.
    # background dict - key <cpg score>: val <set of gene_ids with that score>
    # foreground dict - key <gene_id>: val <cpg score>

    gc = {"background": collections.defaultdict(set), "foreground": {}}
    for tablename in tablenames:

        # MM: need to make sure `-` in filenames don't break the sql statement

        tablename = tablename.replace("-", "_")
        for data in cc.execute("""SELECT * FROM %s;""" % tablename):
            interval_id = data[3].split(" ")[0]
            cpg = data[2]
            print cpg

            # jj: store background in 1 percent bins
            cpg_str = "%.3f" % cpg
            if re.search("background", tablename):
                if interval_id in background_set:
                    gc["background"][cpg_str].add(interval_id)

            elif re.search("foreground", tablename):
                if interval_id in foreground_set:
                    gc["foreground"][interval_id] = cpg_str

            else:
                raise ValueError("Unrecognized table name %s. Should contain"
                                 "'foreground' or 'background'" % tablename)

    # debug: pickle and dump the gc dict
    pickle_file = P.snip(foreground, ".foreground.tsv") + ".p"
    pickle.dump(gc, open(pickle_file, "wb"))

    # match the background set to the foreground set by taking a random
    # background interval with the the same sequence composition as each
    # foreground interval.
    outf = open(outfile, "w")
    if header_line:
        outf.write("gene_id\n")

    # jj: sample background gene_ids without replacement
    matched_background = set()
    X = 0
    for interval, cpg in gc["foreground"].iteritems():

        # print("Finding background for foreground gene: %s (%s)" %
        # (interval, cpg))
        if cpg in gc["background"].keys():

            # get set of bg gene_ids with relevant cpg score
            bg_gene_ids = gc["background"][cpg]

            # print "There are %i background genes in total" % len(bg_gene_ids)
            # remove foreground genes from background set
            bg_gene_ids = bg_gene_ids - foreground_set

            # print("There are %i background genes after removing foreground" %
            # len(bg_gene_ids))

            if bg_gene_ids:
                # select one gene_id to add to matched_background

                bg_id = random.sample(gc["background"][cpg], 1)[0]
                matched_background.add(bg_id)
                # remove selected background gene_id from set

                gc["background"][cpg].remove(bg_id)
            else:
                X += 1
                E.warn("Missing background gene for %s %s, no gene with"
                       " matching %s" % (foreground_file.name,
                                         interval, bg_stat))

        else:
            X += 1
            E.warn("Missing background gene for %s %s, no gene with"
                   " matching %s" % (foreground_file.name, interval, bg_stat))

    # Hack
    # jj: check that background gene_list is <10% shorter than foregroung
    # hack
    # MM: only need to check sufficient background size for Fisher's exact test
    if stat == "fisher":
        assert len(matched_background) > 0.9*len(foreground_set), (
            "There are insufficient genes with matched background to perform"
            " test for sample %s" % foreground_file)
    else:
        pass
    print "Number of genes with no available background: %i" % X
    print "Foreground set: %i" % len(foreground_set)
    print "Backfround set: %i" % len(matched_background)
    outf.write("\n".join(matched_background) + "\n")
    outf.close()


def matchGenesByComposition(bg_gc, fg_gc, bg_stat):
    '''
    bg_gc and fg_gc are pandas DataFrames.

    For each gene in the test (foreground) gene set, generate
    a set of matched background genes based on nucleotide composition
    statistic.  Matching genes are done so on a 2% interval.
    '''

    gc_dict = {}
    match_dict = {}

    props = (x/100.0 for x in xrange(0, 101, 1))

    bg_genes = bg_gc.index
    bg_cpg = bg_gc[bg_stat]
    fg_cpg = fg_gc[bg_stat]

    # bin all background genes into 2% intervals

    for i in props:
        matches = []
        for x in bg_genes:
            poss_ = round(bg_cpg.loc[x], 2)
            if (poss_ >= i) and (poss_ <= (i + 0.02)):
                matches.append(x)
        gc_dict[i] = matches

    for gene in fg_gc.index.tolist():
        fore_pCpG = fg_cpg[gene]
        match_dict[gene] = gc_dict[round(fore_pCpG, 2)]

    return match_dict


def TFBSgeneset(tfbs_id, tfbs_table):
    '''
    Pull out genes for each tfbs from the TFBS table
    to test for enrichment.
    '''

    if len(tfbs_table.loc[tfbs_id]) > 1:
        tfbs_genes = set(tfbs_table.loc[tfbs_id]['seq_id'].tolist())
    elif len(tfbs_table.loc[tfbs_id]) == 1:
        tfbs_genes = set(tfbs_table.loc[tfbs_id]['seq_id'])

    return tfbs_genes


def countTFBSEnrichment(tfbs_genes, gene_set):
    '''
    Count occurances of genes with TFBS
    '''

    in_counter = len(set(gene_set).intersection(tfbs_genes))

    return in_counter, len(tfbs_genes)


def genNullGeneSet(match_background):
    '''
    Generate a null gene set matched to the test gene set
    based on pCpG and GC content.  The null gene set is
    the same size as the test gene set.
    '''

    null_list = set()
    for gene in match_background.keys():
        matches = match_background[gene]
        if len(matches) > 1:
            rand = random.randint(0, len(matches)-1)
            match = matches[rand]
        elif len(matches) == 0:
            pass
        else:
            match = matches[0]

        if match in null_list and len(matches) > 1:
            matches.remove(match)
            rand = random.randint(0, len(matches)-1)
            match = matches[rand]
            null_list.add(match)
        else:
            null_list.add(match)

    return list(null_list)


def nullDistPermutations(tfbs_genes, matched_genes, nPerms=1000):
    '''
    Randomly generate a null distribution of genes from the background
    set, match for pCpG and total number of genes.  Calculate
    enrichment of these genes for the given TFBS.  Permuate nPerm times.
    Return a median enrichment and p-value from the empirical
    cumulative frequency distribution.
    '''

    null_dist = []

    for i in xrange(0, nPerms):
        null_genes = genNullGeneSet(matched_genes)
        null_counts = countTFBSEnrichment(tfbs_genes, null_genes)
        null_enrich = null_counts[0]/float(null_counts[1])
        null_dist.append(null_enrich)

    # null_dist = robjects.FloatVector([float(x) for x in null_dist])
    # null_dist_r = robjects.FloatVector([float(x) for x in null_dist])

    # null_ecdf = ecdf(null_dist)
    out_dict = {}
    out_dict['null'] = null_dist
    out_dict['median'] = np.median(null_enrich)

    return out_dict


def permuteTFBSEnrich(tfbs_table,
                      fg_gc,
                      bg_gc,
                      nPerms,
                      bg_stat):
    '''
    Generate p-value from empirical cumulative frequency distribution
    for all TFBS by permutation.
    '''
    results_dict = {}
    fg_geneset = fg_gc.index.tolist()
    matched_genes = matchGenesByComposition(bg_gc, fg_gc, bg_stat)
    tfbs_all = set(tfbs_table.index)

    ecdf_py = R['ecdf']
    for tfbs in tfbs_all:
        tfbs_res = {}
        tfbs_genes = TFBSgeneset(tfbs, tfbs_table)
        fg_genes_counters = countTFBSEnrichment(tfbs_genes, fg_geneset)

        n_foregenes = fg_genes_counters[0]
        n_tfbsgenes = fg_genes_counters[1]
        total_fg = len(fg_geneset)
        fore_enrich = n_foregenes/float(n_tfbsgenes)
        if fore_enrich > 0:
            null_perms = nullDistPermutations(tfbs_genes=tfbs_genes,
                                              matched_genes=matched_genes,
                                              nPerms=nPerms)
            null_dist = null_perms['null']
            null_median = null_perms['median']
            null_dist_r = robjects.FloatVector([f for f in null_dist])
            null_ecdf = ecdf_py(null_dist_r)
            null_func = null_ecdf(fore_enrich)
            pvalue = 1.00 - null_func[0]

        else:
            null_median = 0
            pvalue = 1.0

        tfbs_res['nForegenes'] = n_foregenes
        tfbs_res['tForegenes'] = total_fg
        tfbs_res['nTFBSGenes'] = n_tfbsgenes
        tfbs_res['propForeInTFBS'] = fore_enrich
        tfbs_res['null_median'] = null_median
        tfbs_res['pvalue'] = pvalue

        results_dict[tfbs] = tfbs_res

    return results_dict


def nCr(n, r):
    '''
    n choose r
    '''

    n_fac = factorial(n)
    r_fac = factorial(r)
    n_r_fac = factorial(n-r)

    ncr = n_fac/r_fac/n_r_fac

    return ncr
