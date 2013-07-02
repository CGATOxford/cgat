#!/home/belgardt/bin/python2.6

# sort genes by relative standard deviation

from sys import argv
from math import sqrt

rpkm_file = open(argv[1],'r')
min_rpkm = float(argv[2])

def calc_rel_std_dev(values):
    variance = 0.0
    mean = float(sum(values))/len(values)
    for val in values:
        variance += (val-mean)**2
    return sqrt(variance/mean)/mean

by_variance=[]
rpkm_file.readline()
for line in rpkm_file:
    la = line.rstrip('\n').split('\t')
    if max( map(float, la[2::]) ) < min_rpkm:
        continue
    by_variance.append( (calc_rel_std_dev( map(float, la[2::]) ), la) )

by_variance.sort()
by_variance.reverse()

for entry in by_variance:
    print str(entry[0]) + "\t" + "\t".join(entry[1])
