#!/home/belgardt/bin/python2.6

# sort genes by mean expression level

from sys import argv

rpkm_file = open(argv[1],'r')

by_mean_expn=[]
rpkm_file.readline()
for line in rpkm_file:
    la = line.rstrip('\n').split('\t')
    by_mean_expn.append( (float(sum( map(float, la[2::]) ))/len(la[2::]), la) )

by_mean_expn.sort()
by_mean_expn.reverse()

for entry in by_mean_expn:
    print str(entry[0]) + "\t" + "\t".join(entry[1])
