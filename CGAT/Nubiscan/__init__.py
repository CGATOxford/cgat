import sys, os, re, collections, math, random
import numpy
import numpy.random

import CGAT.Stats as Stats
import CGAT.Experiment as E

from cnubiscan import *

def iterate_overlaps(matches):

    if len(matches) == 0: return
    #matches.sort( key = lambda x: (x.sequence, x.start) )
    matches.sort()
    lseq, lstart, lend = matches[0].sequence, matches[0].start, matches[0].end
    r = [matches[0]]
    for m in matches[1:]:
        if lseq != m.sequence or m.start > lend:
            yield (lseq, lstart, lend, r)
            lseq, lstart, lend = m.sequence, m.start, m.end
            r = []
        r.append( m )
        lend = max( m.end, lend )
    yield (lseq, lstart, lend, r) 

def combineMotifs( matches ):
    '''merge overlapping motifs.

    In the case of several overlapping motif arrangements, take the one with the
    highest zscore.
    '''
    filtered_results = []
    for seq, start, end, r in iterate_overlaps( matches ):
        r.sort( key = lambda x: x.zscore )
        if len(r) > 0:
            take = r[-1]._replace( alternatives = r[:-1])
        else:
            take = r[-1]
        filtered_results.append( take )
    return filtered_results

def runNubiscanAnalysis( sequences,
                        control_sequences,
                        arrangements,
                        sense_matrix ):

    matcher = Nubiscan( sense_matrix )
    
    all_sequences = sequences + control_sequences

    # find motifs in both foreground and control together
    # result = findNubiscanNuclearReceptors( all_sequences, arrangements )
    result = matcher.run( all_sequences,
                          arrangements,
                          zscore_threshold = 6.5 )
    
    nsequences = len(sequences)
    # separate foreground (fg) and control (background = bg)
    fg_result = [ x for x in result if x.sequence < nsequences ]
    bg_result = [ x for x in result if x.sequence >= nsequences ]

    fg_counter = E.Counter()
    bg_counter = E.Counter()
    fg_seqs = set()
    bg_seqs = set()
    
    for x in fg_result:
        fg_counter[x.arrangement] += 1
        fg_seqs.add( x.sequence )
    for x in bg_result:
        bg_counter[x.arrangement] += 1
        bg_seqs.add( x.sequence )

    print len(fg_result), len(bg_result) / 2
    print len(fg_seqs), len(bg_seqs) / 2

    for x in arrangements:
        print x, fg_counter[x], bg_counter[x] / 2

    fg_filtered = combineMotifs( fg_result )
    bg_filtered = combineMotifs( bg_result )
    
    fg_counter = E.Counter()
    bg_counter = E.Counter()
    fg_seqs = set()
    bg_seqs = set()
    
    for x in fg_filtered:
        fg_counter[x.arrangement] += 1
        fg_seqs.add( x.sequence )
    for x in bg_filtered:
        bg_counter[x.arrangement] += 1
        bg_seqs.add( x.sequence )

    print "###########################################"
    print len(fg_result), len(bg_result) / 2
    print len(fg_seqs), len(bg_seqs) / 2

    for x in arrangements:
        print x, fg_counter[x], bg_counter[x] / 2
    
    return []

def runAnalysis( sequences,
                 arrangements,
                 matrix = 'nr',
                 qvalue_threshold = 0.05 ):

    if matrix == 'nr':
        sense_matrix = NR
    elif matrix == "rxrvdr":
        sense_matrix = RXRVDR
    else:
        raise ValueError("unknown matrix")

    matcher = MatcherRandomisationSequence( sense_matrix )
    
    # find motifs in both foreground and control together
    results = []
    for x, sequence in enumerate(sequences):
        result = matcher.run( sequence,
                              arrangements,
                              qvalue_threshold = qvalue_threshold )
        
        for r in result:
            results.append( r._replace( sequence = x ) )
            
    nsequences = len(sequences)

    fg_filtered = combineMotifs( results )
    fg_counter = E.Counter()
    fg_seqs = set()
    co_counter = E.Counter()
    co_seqs = set()
    
    for x in results:
        fg_counter[x.arrangement] += 1
        fg_seqs.add( x.sequence )
        
    for x in fg_filtered:
        co_counter[x.arrangement] += 1
        co_seqs.add( x.sequence )

    for x in arrangements:
        print x, fg_counter[x], co_counter[x]

    print len(fg_seqs), len(co_seqs)

    return fg_filtered

if __name__ == "__main__":

    sequences = (
        'caaaggtcacagagttcattcatagacctgatgttcaagctctggctccagtaactgccatctgggtcgcttacccagcctgaatctcagtgttctcaac',
        'aaccatctgtgagactggttttgaacccagtgacctctgcctctgggggcagttctgtttccagtggttttgccacagcccagacagggttgaggttgct',
        'CTTCAAATTTGGTAGTAGTTATTTCTGTGAGGTCACTGGGTTCAGACCTTAACCTATCAACTGGAGTGTCTGCTCTAGAGCTAAGTGCTTTGCTGTGGAC',
        "ACATCCTGAACCCTGTGACCTAGGCCCAGGGCTGCTGCGCGGACGGTAGCTCCCCCTGCAGGAAGCAAGGTTCCTCCGGGCCCCCAGACTGCTGCTGGAC",
        "ATGAACTGTATGACCTCCTTTGAGGTCAGTATCCCTGCAGACAAAAATAGGGGCACCACCATTTACATGGCAGATCTGCCATGATGGTGAAAGCACCTGC")

    sequences = open("positive.fasta", "r" ).readlines()

    import random
    ss = []
    for x in sequences:
        xx = list(x)
        random.shuffle(xx)
        ss.append("".join(xx))
    sequences = ss

    sequences = ('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGTAAGTCATATCTCTGAACTTTAACTGTTGTTCCTCTAGCCAGCACCCCACCCCACCAGCTCCCATGTTTGTAAAAGCCTATTTTCCCCTTCCATTCCCTAGAGGGTTGGGTACTCAGGATAGGCCAGAAAGGGCAGGAAGAGAGAGACTCATCCTTACTAAAAGGAAAAGAAAGGGCAACAAGTTCAGCCATGCAGTTCATTTTCACTTTTGAATCCATTCTAAACATGACTAGCCATTATGTTAAAAGCAGGGGAGGTGAGGGGGGCAGGACTAGAAATTAGCATTTCTATATTTATGCTATACTTCAAATTAAACCAGTAATTTCAGTGCCAAAACTATGGAGAATAGGGCTGGGCATGGTGGCTTACACATGTAATCCCAGCATTTGGGGAGGTGGAAGCAGGCTGATC',)
    
    arrangements = [ "DR%s" % x for x in range(0,15) ] + [ "ER%s" % x for x in range(0,15) ]
    for r in runAnalysis( sequences, arrangements, matrix = 'rxrvdr' ):
        print str(r)




    
