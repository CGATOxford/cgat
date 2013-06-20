#!/home/belgardt/bin/python2.6

# sorts patterns by relative standard deviation

from sys import argv
from math import sqrt

rpkm_file = open(argv[1],'r')
association_file = open(argv[2],'r')
min_rpkm = float(argv[3])
min_genes = float(argv[4])

def calc_rel_std_dev(values):
    variance = 0.0
    mean = float(sum(values))/len(values)
    if mean == 0.0:
        return 0.0
    for val in values:
        variance += (val-mean)**2
    return sqrt(variance/mean)/mean

genes_dict = {}
header = rpkm_file.readline().rstrip('\n').split('\t')[2::]
for line in rpkm_file:
    la = line.rstrip('\n').split('\t')
    if max( map(float, la[2::]) ) < min_rpkm:
        continue
    sum_tot = sum( map(float, la[2::]) )
    genes_dict[la[0]] = map(lambda x: float(x)/sum_tot, la[2::])

association_file.readline()
associations_dict = {}
for line in association_file:
    la = line.rstrip('\n').split('\t')
    associations_dict[ la[1] ] = associations_dict.get( la[1], set() )
    associations_dict[ la[1] ].add( la[0] )

by_variance = []
for term in associations_dict.keys():
    term_dist = [0.0] * len( header )
    if len( associations_dict[term].intersection(genes_dict.keys()) ) < min_genes:
        continue
    for gene in associations_dict[term]:
        try:
            term_dist = map(sum, zip(term_dist, genes_dict[gene]) )
        except:
            pass
    by_variance.append( (calc_rel_std_dev( term_dist ), term, len( associations_dict[term].intersection(genes_dict.keys()) ), associations_dict[term].intersection(genes_dict.keys()), term_dist ) )

by_variance.sort()
by_variance.reverse()

for entry in by_variance:
    print "\t".join(map(str,entry[0:4])) + "\t" + "\t".join(map(str,entry[4]))
