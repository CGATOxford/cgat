#!/home/belgardt/bin/python2.6

from sys import stdin

chrom=""
for line in stdin:
    la = line.strip().rstrip('\n').split('\t')
    if not la[0] == chrom:
        chrom = la[0]
        print "variableStep chrom=%s" % chrom
    print "\t".join( la[1:3] )

