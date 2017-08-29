'''
Biomart.py - Utilities for importing data from biomart
==================================================================

This module provides access to biomart servers using the R biomaRt module.

In order to find out about datasets, filters and attributes, use R to query
the server::

   listDatasets(mart)
   listFilters(mart)
   listAttributes(mart)

Reference
---------

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

    Arguments
    ---------
    outfile : string
       Filename of output file
    columns : dict
       Dictionary mapping biomart columns to columns in the output table.
    biomart : string
       Biomart name
    dataset : string
       Biomart dataset
    host : string
       Biomart host
    '''

    R.library("biomaRt")

    keys = list(columns.keys())

    # The default value for host in the biomaRt package is
    # www.biomart.org but for some reason R errors if you specify
    # host manually but then use the default - but it is fine if
    # host is anything valid apart from www.biomart.org.  So I have
    # changed this to only specify a value if the value you
    # are specifying is different to the default KB

    if host == 'www.biomart.org':
        mart = R.useMart(biomart=biomart, dataset=dataset)
    else:
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


    Arguments
    ---------
    columns : dict
       List of fields to obtain.
    biomart : string
       Biomart name
    dataset : string
       Biomart dataset
    host : string
       Biomart host
    filters : list
       List of filter to use
    values : list
       Values of the filters
    archive : bool
       If True, use archived version

    Returns
    -------
    iterator
       Iterator over rows in biomart database. Each row
       is dictionary mapping column names to values.

    '''

    R.library("biomaRt")

    # The default value for host in the biomaRt package is
    # www.biomart.org but for some reason R errors if you specify
    # host manually but then use the default - but it is fine if
    # host is anything valid apart from www.biomart.org.  So I have
    # changed this to only specify a value if the value you
    # are specifying is different to the default KB

    if host == 'www.biomart.org':
        mart = R.useMart(biomart=biomart, dataset=dataset, path=path,
                         archive=archive)
    else:
        mart = R.useMart(biomart=biomart, dataset=dataset, host=host,
                         path=path, archive=archive)

    if filters is not None:
        filter_names = rpy2.robjects.vectors.StrVector(list(filters))
    else:
        filter_names = ""

    if values is not None:
        filter_values = values
    else:
        filter_values = ""

    # result is a dataframe
    result = R.getBM(
        attributes=rpy2.robjects.vectors.StrVector(list(columns)),
        filters=filter_names,
        values=filter_values,
        mart=mart)

    # access via result.rx was broken in rpy2 2.4.2, thus try
    # numeric access
    assert tuple(result.colnames) == tuple(columns),\
        "colnames in dataframe: %s different from expected: %s" % \
        (str(tuple(result.colnames)), tuple(columns))

    for data in zip(*[result[x] for x in range(len(columns))]):
        yield dict(list(zip(columns, data)))
