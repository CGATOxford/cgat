# first add up counts for each gene from all samples
# then find the median ratio of sample X : sum of all samples
# leave the smallest median ratio alone; scale the counts of all others based on that ratio

import sys, optparse, re, os
import numpy

import Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    options, args = E.Start()

    ratios_file=open(args[0],'w')
    samples=args[1:]

    # first read in all values
    gene_dict={}
    for sample in samples:
        numgenes=0
        gene_dict[sample]=[]
        other_gene_info=[]
        for line in open(sample,'r'):
            if line.startswith("#"): continue
            if line.startswith("gene_id"): continue
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
            if totreads==0: continue
            ratio=gene_dict[sample][gene_idx]/float(totreads)
            ratio_dict[sample].append(ratio)

    # find median ratio for each and print the header
    med_ratio=[]
    my_str="gene_id\tlength\t"
    my_str2=""
    for sample in samples:
        my_str+=sample.rsplit('/',1)[1].split('.')[0]+'\t'
        my_str2 += sample.rsplit('/',1)[1].split('.')[0] + '\t'
        ratio_dict[sample].sort()
        med = float(numpy.median(ratio_dict[sample]))
        if med == 0.0:
            E.warn( "median ratio is 0 in %s - added small amount" % sample )
            med += 0.00001

        med_ratio.append( med )

    print my_str.rstrip('\t')
    ratios_file.write( my_str2.rstrip('\t') + '\n' )

    smallest_ratio=min(med_ratio)
    my_str2=""
    for index, sample in enumerate(samples):
        try:
            my_str2 += str( smallest_ratio/med_ratio[index] ) + '\t'
        except ZeroDivisionError:
            my_str2 += 'na\t'
    ratios_file.write( my_str2.rstrip('\t') )

    # identify the smallest median ratio; correct all counts and prints out data from all samples...
    for gene_idx in range(numgenes):
        my_str=other_gene_info[gene_idx][0]+'\t'+other_gene_info[gene_idx][1]+'\t'
        for index, sample in enumerate(samples):
            try:
                my_str+=str(gene_dict[sample][gene_idx]*smallest_ratio/med_ratio[index]) + '\t'
            except ZeroDivisionError:
                my_str += 'na\t'

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

    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
