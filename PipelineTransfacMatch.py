#############################################################
#############################################################
# Classes and functions used in pipeline_trasfacmatch.py
#############################################################
#############################################################

################
# import modules
################
import re, sys, os
import IOTools
import sqlite3
import collections
import Pipeline as P
import numpy as np
from rpy2.robjects import r as R
from rpy2.robjects.packages import importr
from rpy2.robjects.numpy2ri import numpy2ri
import rpy2.robjects as robjects
import glob, random
import FastaIterator

########################################
# helper function
########################################

def filenameToTablename(filename):
    '''
    converts filename containing "." to tablename where "." converted to "_"
    '''
    return filename.replace(".", "_" )



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
        self.seq_id = self.seq_id[len(self.seq_id)-1:][0]

        # iterate over file and yield information for 
        # each sequence identifier
        for line in inf.readlines():
            if line.find("Inspecting sequence ID") != -1:
                self.seq_id = line[:-1].split(" ")
                self.seq_id = self.seq_id[len(self.seq_id)-1:][0]
            else:
                if line.strip() == "": continue
                elif line.find("Total") != -1 or line.find("Frequency") != -1: continue
                self.matrix_id, position_strand, self.core_score, self.matrix_score, self.sequence = [x.strip() for x in line[:-1].split("|")]
                self.position = position_strand.split(" ")[0]
                self.strand = position_strand.split(" ")[1].split("(")[1].split(")")[0]
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
    builds metrics using Match data loaded
    into an sqlite database
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
    for data in cc.execute("""SELECT seq_id, matrix_id FROM %s""" % match_table).fetchall():
        seq_id, tfid = data[0], data[1]
        tf_number[seq_id].add(tfid)
        tf_max[seq_id].append(tfid)

    # get max number of hits for any transcription factor
    for seq_id, tfs in tf_max.iteritems():
        result = collections.defaultdict(int)
        for tf in tfs:
            result[tf] += 1
                    
        # output the results
        outf.write("\t".join(map(str,[seq_id, len(tf_number[seq_id]), max(result.values())])) + "\n")
    outf.close()
            
###################################################
# Significance testing
###################################################

def testSignificanceOfMatrices( background
                               , foreground
                               , database
                               , match_table
                               , outfile ):
    '''
    uses the fishers exact test to estimate the significance of 
    transcription factor motifs in a foregound set of intervals
    compared to a background set
    '''
    # get background and foreground sets                                                                                                                                                                                                     
    background_file = open(background)
    foreground_file = open(foreground)
    interval_sets = {"foreground":set(), "background": set()}
    for inf_name in ["background_file", "foreground_file"]:
        name = P.snip(inf_name, "_file")
        eval(inf_name).readline()
        for line in eval(inf_name).readlines():
            interval_sets[name].add(line[:-1])
 
    # connect to database                                                                                                                                                                                                                    
    dbh = sqlite3.connect(database)
    cc = dbh.cursor()

    # retrieve a set of sequence ids and the transcription                                                                                                                                                                                   
    # factors that they bind                                                                                                                                                                                                                 
    # NB this will only add 1 count regardless of the number of 
    # motifs per sequence
    ids_tfs = collections.defaultdict(set)
    all_tfs = set()
    for data in cc.execute("""SELECT seq_id, matrix_id FROM %s""" % match_table).fetchall():
        seq_id, tf = data[0], data[1]
        ids_tfs[seq_id].add(tf)
        all_tfs.add(tf)

    outf = open(outfile, "w")
    outf.write("matrix_id\toddsRatio\t95CI_low\t95_CI_hi\tpvalue\tnforeground\tnbackground\ttforeground\ttbackground\tqvalue\n")
    
    # iterate over factors and calculate the significance                                                                                                                                                                                    
    # using fishers exact test                                                                                                                                                                                                               
    fisherpy = R["fisher.test"]
    padjpy = R["p.adjust"]

    results = []
    pvalues = []
    for tf in all_tfs:
        # initialise contingency table                                                                                                                                                                                                       
        contingency = np.zeros((2,2))
        temp = P.getTempFile()
        for seq_id in ids_tfs.keys():
            if seq_id in interval_sets["foreground"]:
                if tf in ids_tfs[seq_id]:
                    contingency[0,0] += 1
                elif tf not in ids_tfs[seq_id]:
                    contingency[1,0] += 1
            elif seq_id in interval_sets["background"]:
                if tf in ids_tfs[seq_id]:
                    contingency[0,1] += 1
                elif tf not in ids_tfs[seq_id]:
                    contingency[1,1] += 1
        f =  fisherpy(numpy2ri(contingency))
        
        # get overlap numbers
        nforeground = contingency[0,0]
        nbackground = contingency[0,1]
        tforeground = len(interval_sets["foreground"])
        tbackground = len(interval_sets["background"])
        
        # convert back to python object                                                                                                                                                                                                      
        f = [list(x) for x in np.array(f)]
        pvalue, ci_low, ci_hi, OR =  f[0][0], f[1][0], f[1][1], f[2][0]
        pvalues.append(pvalue)
        results.append([tf, OR, ci_low, ci_hi, pvalue, nforeground, nbackground, tforeground, tbackground])

    # correct for multiple comparisons
    qvalues = padjpy(pvalues)
    for i in range(len(results)):
        outf.write("\t".join( map(str, results[i] + [qvalues[i]]) ) + "\n")
    outf.close()
                   
######################################################
# CpG matching                   
######################################################

def calculateCpGComposition(interval_names, sequence_file, outfile):
    '''
    given a set of interval names that are present in a 
    fasta file, return CpG content file
    '''
    interval_file = open(interval_names)
    interval_file.readline()
    sequence_file = open(sequence_file)
    
    interval_set = set()
    for line in interval_file.readlines():
        interval_set.add(line[:-1])
        
    temp = P.getTempFile()
    for record in FastaIterator.iterate(sequence_file):
        seq_id = record.title.split(" ")[0]
        if seq_id in interval_set:
            temp.write(">%s\n%s\n" % (record.title, record.sequence))
    temp.close()

    inf = temp.name 
    statement = '''cat %(inf)s | python %(scriptsdir)s/fasta2table.py
                   -s cpg -s length 
                   --log=%(outfile)s.log > %(outfile)s'''
    P.run()

####################################################
####################################################
def matchBackgroundForCpGComposition( gc_load_files
                                      , background
                                      , foreground
                                      , database
                                      , outfile ):
        '''
        take the background set and subset it for intervals with
        a CpG distribution that is the same as the foreground set
        - this requires that the background set is sufficiently
        large
        '''
        background_file = open(background)
        background_file.readline()
        foreground_file = open(foreground)
        foreground_file.readline()
        fasta_file = open(glob.glob("fasta.dir/*.fasta")[0]) # NB this is specific to this pipeline
        
        background_set = set()
        foreground_set = set()
        for interval_id in background_file.readlines():
            background_set.add(interval_id[:-1])
        for interval_id in foreground_file.readlines():
            foreground_set.add(interval_id[:-1])

        dbh = sqlite3.connect(database)
        cc = dbh.cursor()
        tablenames = [filenameToTablename(P.snip(os.path.basename(x), ".load")) for x in gc_load_files]

        # NB background set will be keyed on CpG proportion
        gc = {"background":collections.defaultdict( set ), "foreground":{}}
        for tablename in tablenames:
            for data in cc.execute("""SELECT id, pCpG FROM %s""" % tablename):
                interval_id, cpg = data[0].split(" ")[0], data[1]
                if interval_id in background_set:
                    gc["background"][cpg].add( interval_id )
                elif interval_id in foreground_set:
                    gc["foreground"][interval_id] = cpg
            
        # match the background set to the foreground set 
        # by taking a random background interval with the 
        # the same CpG content as each foreground interval
        # this will match the length of each list
        outf = open(outfile, "w")
        outf.write("gene_id\n")
        matched_background = set()
        for interval, cpg in gc["foreground"].iteritems():
            if gc["foreground"][interval] in gc["background"].keys():
                outf.write( random.sample( gc["background"][cpg], 1 )[0] + "\n") 
        outf.close()



