#!/home/belgardt/bin/python2.6

from sys import argv
import pysam
from commands import getoutput

coords_file=open(argv[1],'r')
bamfile=pysam.Samfile( argv[2], 'rb' )  # bamfile

oldgene=''
genelen=0
num_reads=0
first=True
for line in coords_file:
    la=line.rstrip('\n').split('\t')
    if len(la)<4:
        break
    chrom=la[0]
    if la[3]==oldgene:  # if still on the same gene
        genelen+=int(la[2])-int(la[1])+1
        coords=chrom+':'+str(la[1])+'-'+str(la[2])
        for alignedread in bamfile.fetch(coords):
            if not (alignedread.qname, alignedread.is_read1) in anames:
                anames.add((alignedread.qname, alignedread.is_read1))
                num_reads+=1
    else:   # if on a new gene
        if first:
            first=False
        else:
            print oldgene + '\t' + str(num_reads) + '\t' + str(genelen)
            num_reads=0
        oldgene=la[3]
        genelen=int(la[2])-int(la[1])+1
        coords=chrom+':'+la[1]+'-'+la[2]
        anames=set([])
        for alignedread in bamfile.fetch(coords):
            if not (alignedread.qname, alignedread.is_read1) in anames:
                anames.add((alignedread.qname, alignedread.is_read1))
                num_reads+=1
print oldgene + '\t' + str(num_reads) + '\t' + str(genelen)
