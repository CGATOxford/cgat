#!/home/belgardt/bin/python2.6

import sys, optparse

import Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )


    infile = open(args[0],'r')
    genes_list = []
    header = None

    for line in infile:
        if line.startswith("#"): continue
        if line.startswith("gene_id"):
            header = line.rstrip('\n')
            num_fields = len(header.split('\t')) - 2
            total_reads = [0] * num_fields
            continue

        la = line.rstrip('\n').split('\t')
        if len(la)<3:
            continue
        genes_list.append(la)
        total_reads = map(lambda x, y: float(x)+float(y), total_reads, la[2::])

    total_reads = map(lambda x: float(x)/1000000, total_reads)

    print header

    for gene in genes_list:
        my_str_list = gene[0:2]
        vals = map(lambda x, y: float(x)/float(y)/(float(gene[1])/1000.0), gene[2::], total_reads)
        my_str_list.extend( map(str, vals ) )
        print "\t".join(my_str_list)
        
    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
