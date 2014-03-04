'''
PipelineKegg.py - Pipeline components related to KEGG analysis 
==============================================================
'''

import CGAT.Experiment as E
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import CGAT.IOTools as IOTools
import re


def importKEGGAssignments(outfile, mart, host, biomart_dataset):
    ''' import the KEGG annotations from the R KEGG.db 
    annotations package. Note that since KEGG is no longer
    publically availible, this is not up-to-date and maybe removed
    from bioconductor in future releases '''

    R.library("KEGG.db")
    R.library("biomaRt")

    E.info("getting entrez to ensembl mapping ...")
    mart = R.useMart(biomart=mart,
                     host=host,
                     path="/biomart/martservice",
                     dataset=biomart_dataset)

    entrez2ensembl = R.getBM(attributes=ro.StrVector(["ensembl_gene_id", "entrezgene"]),
                             mart=mart)

    entrez = entrez2ensembl.rx2("entrezgene")
    ensembl = entrez2ensembl.rx2("ensembl_gene_id")
    entrez2ensembl = dict(zip(entrez, ensembl))

    E.info("Done")

    E.info("getting entrez to kegg mapping ... ")
    entrez2path = R('as.list(KEGGEXTID2PATHID)')
    E.info("Done")

    E.info("Getting KEGG names")
    pathnames = R('as.list(KEGGPATHID2NAME)')
    pathid2name = dict(zip(pathnames.names, R.unlist(pathnames)))
    E.info("Done")

    outf = IOTools.openFile(outfile, "w")
    outf.write("ontology\tgene_id\tkegg_ID\tkegg_name\tevidence\n")

    for gene in entrez2path.names:

        try:
            gene = int(gene)
        except ValueError:
            continue

        if gene in entrez2ensembl:
            ensid = entrez2ensembl[gene]

        else:
            continue

        for pathway in entrez2path.rx2(str(gene)):
            pathid = re.match("[a-z]+([0-9]+)", pathway).groups()[0]
            pathname = pathid2name[pathid]
            outf.write(
                "\t".join(["kegg", ensid, str(pathway), pathname, "NA"]) + "\n")
