'''utility tasks for importing data from biomart.
'''

from rpy import r as R
import rpy

def importFromBiomart( outfile, columns, 
                       biomart = "ensembl", 
                       dataset = "hsapiens_gene_ensembl" ):
    '''download a dataset from biomart and output as a 
    tab-separated table.

    *columns* is a dictionary mapping biomart columns to columns
    in the output tables. *biomart* and *dataset* denote the 
    database and dataset to get the data from.

    '''

    R.library("biomaRt")
    
    keys = columns.keys()

    mart = R.useMart(biomart=biomart, dataset= dataset )
    result = R.getBM( attributes=keys, mart=mart )
    
    outf = open( outfile, "w" )
    outf.write( "\t".join( [columns[x] for x in keys ] ) + "\n" )
    
    # for x in ("mim_gene_accession", "mim_morbid_accession"):
    #     result[x] = [ ("", y)[y >= 0] for y in result[x] ]

    for data in zip( *[ result[x] for x in keys] ):
        outf.write( "\t".join( map(str, data) ) + "\n" )

    outf.close()

def biomart_iterator( columns, 
                      biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl" ):
    '''download a dataset from biomart and output as a 
    tab-separated table.

    *columns* is a list with field to obtain.
    
    *biomart* and *dataset* denote the 
    database and dataset to get the data from.

    returns a iterator over rows.
    '''

    R.library("biomaRt")

    mart = R.useMart(biomart=biomart, dataset= dataset )
    result = R.getBM( attributes=columns, mart=mart )
    
    for data in zip( *[ result[x] for x in columns] ):
        yield dict( zip(columns, data) )
        

