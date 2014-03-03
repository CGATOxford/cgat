"""

:Author: Andreas Heger
:Release: $Id: PipelineGO.py 2877 2010-03-27 17:42:26Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

Pipeline components - GO analysis

Tasks related to gene set GO analysis.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----


"""
import sys
import tempfile
import optparse
import shutil
import itertools
import csv
import math
import random
import re
import glob
import os
import shutil
import collections
import sqlite3

import CGAT.Experiment as E
import CGAT.Pipeline as P
import CGAT.Stats as Stats
import CGAT.IOTools as IOTools
import CGAT.CSV as CSV

try:
    PARAMS = P.getParameters()
except IOError:
    pass

############################################################
############################################################
############################################################
# get GO assignments


def createGOFromENSEMBL(infile, outfile):
    '''get GO assignments from ENSEMBL'''

    job_options = "-l mem_free=5G"
    statement = '''
        python %(scriptsdir)s/runGO.py 
                     --filename-dump=%(outfile)s 
                     --host=%(go_host)s 
                     --user=anonymous 
                     --database=%(go_database)s 
                     --port=%(go_port)i > %(outfile)s.log
        '''

    P.run()

############################################################
############################################################
############################################################
# get go assignments


def createGOFromGeneOntology(infile, outfile):
    '''get GO assignments from Geneontology.org

    GO terms are mapped to ensembl gene names via uniprot identifiers.
    '''

    filename = "geneontology.goa.gz"
    if not os.path.exists(filename):
        statement = '''
    wget -O %(filename)s http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/%(go_geneontology_file)s?rev=HEAD
    '''

        P.run()

    # see http://www.geneontology.org/gene-associations/readme/goa.README
    Data = collections.namedtuple("Data", "db db_object_id db_object_symbol qualifier goid dbreference evidence "
                                  " with_id aspect "
                                  " db_object_name synonym db_object_type "
                                  " taxon_id date assigned_by "
                                  " annotation_extension"
                                  " gene_product_form_id")

    dbh = sqlite3.connect(PARAMS["database"])
    cc = dbh.cursor()
    map_uniprot2ensembl = dict(
        cc.execute("SELECT DISTINCT gene_name, gene_id FROM transcript_info").fetchall())
    map_goid2description = dict(
        cc.execute("SELECT DISTINCT go_id, description FROM go_assignments").fetchall())

    aspect2name = {"P": "biol_process",
                   "F": "mol_function",
                   "C": "cell_location"}

    c = E.Counter()
    found_uniprot, found_genes, notfound_uniprot = set(), set(), set()
    outf = IOTools.openFile(outfile, "w")
    outf.write("go_type\tgene_id\tgo_id\tdescription\tevidence\n")
    for line in IOTools.openFile(filename):
        if line.startswith("!"):
            continue
        c.input += 1
        data = Data._make(line[:-1].split("\t"))

        if data.db_object_symbol in map_uniprot2ensembl:
            gene_id = map_uniprot2ensembl[data.db_object_symbol]
            found_uniprot.add(data.db_object_symbol)
            found_genes.add(gene_id)
            outf.write("%s\t%s\t%s\t%s\t%s\n" %
                       (aspect2name[data.aspect],
                        gene_id,
                        data.goid,
                        map_goid2description.get(data.goid, ""),
                        data.evidence))
            c.output += 1

        else:
            c.notfound += 1
            notfound_uniprot.add(data.db_object_symbol)

    c.found_genes = len(found_genes)
    c.found_uniprot = len(found_uniprot)
    c.notfound_uniprot = len(notfound_uniprot)

    E.info("%s" % str(c))
    E.info("not found=%s" % str(notfound_uniprot))
    outf.close()

############################################################
############################################################
############################################################
# get GO descriptions
############################################################


def imputeGO(infile_go, infile_paths, outfile):
    '''impute GO accessions.

    Infile is a file with GO assocations and a file
    with paths of term to ancester (see go2fmt.pl).
    '''

    c = E.Counter()

    term2ancestors = collections.defaultdict(set)
    with IOTools.openFile(infile_paths) as inf:
        for line in inf:
            parts = line[:-1].split()
            term = parts[0]
            ancestors = [parts[x] for x in range(2, len(parts), 2)]
            # there can be multiple paths
            term2ancestors[term].update(ancestors)

    goid2description = {}
    gene2goids = collections.defaultdict(list)
    goid2type = {}
    with IOTools.openFile(infile_go) as inf:
        for line in inf:
            if line.startswith("go_type"):
                continue
            go_type, gene_id, goid, description, evidence = line[
                :-1].split("\t")
            gene2goids[gene_id].append(goid)
            goid2description[goid] = description
            goid2type[goid] = go_type

    outf = IOTools.openFile(outfile, "w ")
    for gene_id, in_goids in gene2goids.iteritems():
        c.genes += 1
        out_goids = set(in_goids)
        for goid in in_goids:
            out_goids.update(term2ancestors[goid])
        if len(in_goids) != len(out_goids):
            c.increased += 1
        else:
            c.complete += 1

        for goid in out_goids:
            outf.write("\t".join(
                (goid2type.get(goid, ""), gene_id, goid, goid2description.get(goid, ""), "NA")) + "\n")
            c.assocations += 1

    outf.close()

    E.info("%s" % str(c))

############################################################


def buildGOPaths(infile, outfile):
    '''output file with paths of terms to root.

    infile is an ontology obo file.
    '''
    use_cluster = True
    statement = '''
    go2fmt.pl -w pathlist %(infile)s > %(outfile)s
    '''
    P.run()

############################################################


def buildGOTable(infile, outfile):
    '''output file with paths of terms to root.

    infile is an ontology obo file.
    '''
    use_cluster = True
    statement = '''
    echo -e "go_id\\tdescription\\tlong_description\\ttext\\n" > %(outfile)s;
    go2fmt.pl -w tbl %(infile)s >> %(outfile)s
    '''
    P.run()

############################################################
############################################################
############################################################
# get GO descriptions
############################################################


def getGODescriptions(infile):
    '''return dictionary mapping GO category to description
    and namespace.
    '''

    with IOTools.openFile(infile) as inf:
        fields, table = CSV.ReadTable(inf, as_rows=False)

    return dict([(y, (x, z)) for x, y, z in zip(table[fields.index("go_type")],
                                                table[fields.index("go_id")],
                                                table[fields.index("description")])])


############################################################
############################################################
############################################################
# get GO Slim assignments
############################################################
def createGOSlimFromENSEMBL(infile, outfile):
    '''get GO assignments from ENSEMBL'''

    statement = '''wget %(go_url_goslim)s --output-document=goslim.obo'''
    P.run()

    statement = '''wget %(go_url_ontology)s --output-document=go_ontology.obo'''
    P.run()

    to_cluster = True
    job_options = "-l mem_free=5G"
    statement = '''
        map2slim -outmap %(outfile)s.map goslim.obo go_ontology.obo
    '''
    P.run()

    job_options = "-l mem_free=5G"
    statement = '''
        zcat < %(infile)s
        | python %(scriptsdir)s/runGO.py 
                --go2goslim 
                --filename-ontology=go_ontology.obo 
                --slims=%(outfile)s.map 
                --log=%(outfile)s.log 
        | gzip
        > %(outfile)s
        '''
    P.run()

############################################################


def runGOFromFiles(outfile,
                   outdir,
                   fg_file,
                   bg_file=None,
                   go_file=None,
                   ontology_file=None,
                   samples=None,
                   minimum_counts=0,
                   pairs=False,
                   gene2name=None):
    '''check for GO enrichment within a gene list.

    The gene list is given in ``fg_file``. It is compared
    against ``bg_file`` using the GO assignments from
    ``go_file``. Results are saved in ``outfile`` and
    ``outdir``.

    if *bg_file* is None, the all genes with GO annotations
    will be used.

    If *gene2name* is given, it will be supplied to the runGO.py script.

    If *pairs* is set, each category for each pair of gene sets will be tested for 
    differential enrichment.
    '''

    to_cluster = True

    if ontology_file is None:
        ontology_file = PARAMS.get("go_ontology", None)

    options = []
    if ontology_file:
        options.append("--filename-ontology=%(ontology_file)s" % locals())

    if bg_file is not None:
        options.append("--background=%(bg_file)s" % locals())

    if samples is not None:
        options.append("--fdr")
        options.append("--sample=%(samples)i" % locals())
        options.append("--qvalue-method=empirical")
    else:
        options.append("--fdr")
        options.append("--qvalue-method=BH")

    if pairs:
        options.append("--pairwise")

    if gene2name:
        options.append("--filename-gene2name=%s" % gene2name)

    options = " ".join(options)
    statement = '''
    python %(scriptsdir)s/runGO.py 
        --filename-input=%(go_file)s 
        --genes=%(fg_file)s 
        --output-filename-pattern='%(outdir)s/%%(set)s.%%(go)s.%%(section)s' 
        --minimum-counts=%(minimum_counts)i 
        --log=%(outfile)s.log
        %(options)s
    > %(outfile)s'''

    P.run()

    dbhandle = sqlite3.connect(PARAMS["database"])

############################################################


def runGOFromDatabase(outfile, outdir,
                      statement_fg,
                      statement_bg,
                      go_file,
                      ontology_file=None,
                      samples=1000):
    '''Take gene lists from the SQL database using
    ``statement_foreground`` and ``statement_background``
    '''

    dbhandle = sqlite3.connect(PARAMS["database"])

    cc = dbhandle.cursor()
    fg = set([x[0] for x in cc.execute(statement_fg).fetchall()])
    bg = set([x[0] for x in cc.execute(statement_bg).fetchall()])

    if len(fg) == 0:
        P.touch(outfile)
        return

    fg_file = os.path.join(outdir, "foreground")
    bg_file = os.path.join(outdir, "background")
    outf = open(fg_file, "w")
    outf.write("\n".join(map(str, fg)) + "\n")
    outf.close()
    outf = open(bg_file, "w")
    outf.write("\n".join(map(str, bg)) + "\n")
    outf.close()

    runGOFromFiles(outfile, outdir,
                   fg_file, bg_file,
                   go_file,
                   ontology_file=ontology_file,
                   samples=samples)

############################################################
############################################################
############################################################
##
############################################################


def loadGO(infile, outfile, tablename):
    '''import GO results into individual tables.'''

    indir = infile + ".dir"

    if not os.path.exists(indir):
        P.touch(outfile)
        return

    statement = '''
    python %(toolsdir)s/cat_tables.py %(indir)s/*.overall |\
    python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --allow-empty \
              --index=category \
              --index=goid \
              --table=%(tablename)s \
    > %(outfile)s
    '''
    P.run()


############################################################
############################################################
############################################################
##
############################################################
def loadGOs(infiles, outfile, tablename):
    '''import GO results into a single table.

    This method also computes a global QValue over all
    tracks, genesets and annotation sets.
    '''

    header = False

    tempf1 = P.getTempFile()

    pvalues = []

    for infile in infiles:
        indir = infile + ".dir"

        if not os.path.exists(indir):
            continue

        track, geneset, annotationset = re.search(
            "^(\S+)_vs_(\S+)\.(\S+)", infile).groups()

        for filename in glob.glob(os.path.join(indir, "*.overall")):
            for line in open(filename, "r"):
                if line.startswith("#"):
                    continue
                data = line[:-1].split("\t")
                if line.startswith("code"):
                    if header:
                        continue
                    tempf1.write("track\tgeneset\tannotationset\t%s" % line)
                    header = True
                    assert data[10] == "pover" and data[
                        11] == "punder", "format error, expected pover-punder, got %s-%s" % (data[10], data[11])
                    continue
                tempf1.write("%s\t%s\t%s\t%s" %
                             (track, geneset, annotationset, line))
                pvalues.append(min(float(data[10]), float(data[11])))

    tempf1.close()

    E.info("analysing %i pvalues" % len(pvalues))
    fdr = Stats.doFDR(pvalues)
    E.info("got %i qvalues" % len(fdr.mQValues))
    qvalues = ["global_qvalue"] + fdr.mQValues

    tempf2 = P.getTempFile()

    for line, qvalue in zip(open(tempf1.name, "r"), qvalues):
        tempf2.write("%s\t%s\n" % (line[:-1], str(qvalue)))

    tempf2.close()
    tempfilename = tempf2.name
    print tempf1.name
    print tempf2.name

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              --allow-empty 
              --index=category 
              --index=track,geneset,annotationset
              --index=geneset
              --index=annotationset
              --index=goid 
              --table=%(tablename)s 
    < %(tempfilename)s
    > %(outfile)s
    '''
    P.run()

    #os.unlink( tempf1.name )
    #os.unlink( tempf2.name )
