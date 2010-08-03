#!/home/belgardt/bin/python2.6

from sys import argv

infile = open(argv[1],'r')

header = infile.readline().rstrip('\n')
print header
num_fields = len(header.split('\t')) - 2

total_reads = [0] * num_fields
genes_list = []
for line in infile:
    la = line.rstrip('\n').split('\t')
    if len(la)<3:
        continue
    genes_list.append(la)
    total_reads = map(lambda x, y: float(x)+float(y), total_reads, la[2::])
total_reads = map(lambda x: float(x)/1000000, total_reads)

for gene in genes_list:
    my_str_list = gene[0:2]
    my_str_list.extend( map(lambda x, y: str(float(x)/float(y)/(float(gene[1])/1000)), gene[2::], total_reads) )
    print "\t".join(my_str_list)
