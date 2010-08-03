#!/home/belgardt/bin/python2.6

# first add up counts for each gene from all samples
# then find the median ratio of sample X : sum of all samples
# leave the smallest median ratio alone; scale the counts of all others based on that ratio

from numpy import median
from sys import argv

ratios_file=open(argv[1],'w')
samples=argv[2::]

# first read in all values
gene_dict={}
for sample in samples:
    numgenes=0
    gene_dict[sample]=[]
    other_gene_info=[]
    file=open(sample,'r')
    for line in file:
        la=line.rstrip('\n').split('\t')
        if not len(la)==3:  # for example, when qsub spits something dumb out
            continue
        other_gene_info.append((la[0],la[2]))
        gene_dict[sample].append(int(la[1]))
        numgenes+=1

# then sum all samples and calculate ratios
ratio_dict={}
for sample in samples:  # initialize ratios dicts
    ratio_dict[sample]=[]
for gene_idx in range(numgenes):
    totreads=0
    for sample in samples:
        totreads+=gene_dict[sample][gene_idx]
    for sample in samples:
        if totreads==0:
            continue
        ratio=gene_dict[sample][gene_idx]/float(totreads)
        ratio_dict[sample].append(ratio)

# find median ratio for each and print the header
med_ratio=[]
my_str="Name\tLength\t"
my_str2=""
for sample in samples:
    my_str+=sample.rsplit('/',1)[1].split('.')[0]+'\t'
    my_str2 += sample.rsplit('/',1)[1].split('.')[0] + '\t'
    ratio_dict[sample].sort()
    med_ratio.append(float(median(ratio_dict[sample])))
print my_str.rstrip('\t')
ratios_file.write( my_str2.rstrip('\t') + '\n' )

smallest_ratio=min(med_ratio)
my_str2=""
for index, sample in enumerate(samples):
    my_str2 += str( smallest_ratio/med_ratio[index] ) + '\t'
ratios_file.write( my_str2.rstrip('\t') )

# identify the smallest median ratio; correct all counts and prints out data from all samples...
for gene_idx in range(numgenes):
    my_str=other_gene_info[gene_idx][0]+'\t'+other_gene_info[gene_idx][1]+'\t'
    for index, sample in enumerate(samples):
        my_str+=str(gene_dict[sample][gene_idx]*smallest_ratio/med_ratio[index]) + '\t'
    print my_str.rstrip('\t')

"""
files=[]
for index, sample in samples:
    files.append(open(sample,'r'))
for line in file[0]:
    la=line.rstrip('\n').split('\t')
    gene_name=la[0]
    my_str=gene_name + '\t' + str( float(la[1])*smallest_ratio/med_ratio[0] )
    for each in files[1::]:
        each.readline()
    output.write(my_str)
"""
