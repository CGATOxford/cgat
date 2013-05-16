#!/home/belgardt/bin/python2.6
# computes enrichments in gene categories

from sys import argv
from rpy import *
from math import log

foreground_file = open( argv[1], 'r' )  # gene identifiers, one per line, of genes in the foreground set
background_file = open( argv[2], 'r' )  # gene identifiers, one per line, of genes in the background set
association_file = open( argv[3], 'r' ) # gene identifiers and associated term identifiers, tab separated, one pair per line
summary_file = open( argv[4], 'w' )     # summary results
expanded_file = open( argv[5], 'w' )    # results with all associated gene identifiers

foreground_genes = set([])
for line in foreground_file:
    if line.rstrip('\n') == '':
        continue
    foreground_genes.add(line.rstrip('\n'))

background_genes = set([])
for line in background_file:
    if line.rstrip('\n') == '':
        continue
    background_genes.add(line.rstrip('\n'))

association_to_genes = {}
associations = {}
association_file.readline() # get rid of header line
for line in association_file:
    la = line.rstrip('\n').split('\t')
    if len(la) < 2:
        continue
    associations[ la[0] ] = associations.get( la[0], set([]) )
    associations[ la[0] ].add( la[1] )
    association_to_genes[ la[1] ] = association_to_genes.get( la[1], set([]) )
    association_to_genes[ la[1] ].add( la[0] )

genes_in_fore = 0; genes_in_back = 0
foreground_associations = []
for gene in foreground_genes:
    genes_in_fore += 1
    try:
        foreground_associations.extend( list( associations[ gene ] ) )
    except:
        pass
background_associations = []
for gene in background_genes:
    genes_in_back += 1
    try:
        background_associations.extend( list( associations[ gene ] ) )
    except:
        pass
background_set = set(background_associations)

preresults = []
for association in background_set:
    fore_genes_associated = foreground_associations.count( association )
    back_genes_associated = background_associations.count( association )
    Fisher_table = r.matrix(r.c( back_genes_associated, genes_in_back - back_genes_associated, fore_genes_associated, genes_in_fore - fore_genes_associated),nr=2)
    Fisher_results = r.fisher_test( Fisher_table )
    try:
        log_fold_diff = log( (float(fore_genes_associated)/float(genes_in_fore)) / (float(back_genes_associated)/float(genes_in_back)), 2 )
    except:
        log_fold_diff = "TINY"
    preresults.append( ( float(Fisher_results['p.value']), association, log_fold_diff, (back_genes_associated, genes_in_back, fore_genes_associated, genes_in_fore), association_to_genes[association].intersection(foreground_genes) ) )
preresults.sort()

summary_file.write( "Name\tChange\tLog2(Fold Difference)\tFDR\tSingle test p-value\n" )
expanded_file.write( "Name\tChange\tLog2(Fold Difference)\tFDR\tSingle test p-value\tGenes in this association\n" )
for num_terms, result in enumerate( preresults ):
    FDR = result[0]*len(set(background_associations))/(1+num_terms)
    summ_str = "%s\t(%s/%s)->(%s/%s)\t" % ( result[1], result[3][0], result[3][1], result[3][2], result[3][3] ) + "\t".join( [str(result[2]), str(FDR), str(result[0]) ] ) + '\n'
    summary_file.write( summ_str )
    exp_str = summ_str.rstrip('\n') + '\t%s' % ( ",".join( list(result[4]) ) ) + '\n'
    expanded_file.write( exp_str )

