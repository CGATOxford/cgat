'''
Biomart.py - utility tasks for importing data from biomart
==================================================================

In order to find out about datasets, filters and attributes, use the

listDatasets(mart)
listFilters(mart)
listAttributes(mart)

functions.
'''

import CGAT.IOTools as IOTools

from rpy2.robjects import r as R
import rpy2.robjects.numpy2ri


def importFromBiomart(outfile,
                      columns,
                      biomart="ensembl",
                      dataset="hsapiens_gene_ensembl",
                      host='www.biomart.org'):
    '''download a dataset from biomart and output as a
    tab-separated table.

    *columns* is a dictionary mapping biomart columns to columns in
    the output tables. *biomart* and *dataset* denote the database and
    dataset to get the data from.

    '''

    R.library("biomaRt")

    keys = columns.keys()

    mart = R.useMart(biomart=biomart, dataset=dataset, host=host)
    result = R.getBM(attributes=keys, mart=mart)

    outf = IOTools.openFile(outfile, "w")
    outf.write("\t".join([columns[x] for x in keys]) + "\n")

    # for x in ("mim_gene_accession", "mim_morbid_accession"):
    #     result[x] = [ ("", y)[y >= 0] for y in result[x] ]

    for data in zip(*[result[x] for x in keys]):
        outf.write("\t".join(map(str, data)) + "\n")

    outf.close()


def biomart_iterator(columns,
                     biomart="ensembl",
                     dataset="hsapiens_gene_ensembl",
                     host='www.biomart.org',
                     path="/biomart/martservice",
                     filters=None,
                     values=None,
                     archive=False):
    '''download a dataset from biomart and output as a
    tab-separated table.

    *columns* is a list with field to obtain.

    *biomart* and *dataset* denote the
    database and dataset to get the data from.

    returns a iterator over rows.
    '''

    R.library("biomaRt")

    mart = R.useMart(biomart=biomart,
                     dataset=dataset,
                     host=host,
                     path=path,
                     archive=archive)

    if filters is not None:
        filter_names = rpy2.robjects.vectors.StrVector(filters)
    else:
        filter_names = ""

    if values is not None:
        filter_values = values
    else:
        filter_values = ""

    # result is a dataframe
    result = R.getBM(
        attributes=rpy2.robjects.vectors.StrVector(columns),
        filters=filter_names,
        values=filter_values,
        mart=mart)

    # access via result.rx was broken in rpy2 2.4.2, thus try
    # numeric access
    assert tuple(result.colnames) == tuple(columns),\
        "colnames in dataframe: %s different from expected: %s" % \
        (str(tuple(result.colnames)), tuple(columns))

    for data in zip(*[result[x] for x in range(len(columns))]):
        yield dict(zip(columns, data))
