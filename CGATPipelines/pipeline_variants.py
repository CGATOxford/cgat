"""
===========================
Variant annotation pipeline
===========================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

The Variants pipeline attempts to annotate variants in
a :term:`vcf` formatted file. It computes
 
   1. the effects of SNPs on transcripts and genes

This pipeline works on a single genome.

Overview
========

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineReporting`). To start with, use the files supplied with the
Example_ data.

Input
-----

Variants
++++++++

Variants are read from a :term:`vcf` formatted file called :file:`variants.vcf.gz`. 
The file is assumed to have been compressed with :file:`bgzip` and compressed with
tabix.

The tracks are taken from the headers in the :term:`vcf` file. Please avoid any special
characters like ``_][*.+-``  within strain names.

The pipeline expects the following information within the genotype field in 
the :term:`vcf` file:

GT
   The genotype

DP
   The read depth




Optional inputs
+++++++++++++++

Requirements
------------

The pipeline requires the results from :doc:`pipeline_annotations`. Set the configuration variable 
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|polyphen_           |>=2.0.23           |prediction of deleterious substitutions         |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The major output is in the database file :file:`csvdb`.

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_variants.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_variants.tgz
   tar -xvzf pipeline_variants.tgz
   cd pipeline_variants
   python <srcdir>/pipeline_variants.py make full

.. note:: 
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::

   polyphen
       polyphen_ - a program to predict the deleteriousness of substitutions

   
.. _polyphen: http://genetics.bwh.harvard.edu/pph2/dokuwiki/start

Code
====

"""
from ruffus import *
import sys
import glob
import gzip
import os
import itertools
import CGAT.CSV as CSV
import re
import math
import collections
import numpy
import sqlite3
import CGAT.GTF as GTF
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Database as Database
import CGAT.FastaIterator as FastaIterator
import PipelineEnrichment as PEnrichment
import PipelineBiomart as PBiomart
import PipelineDatabase as PipelineDatabase
import CGATPipelines.PipelineUCSC as PipelineUCSC
import scipy.stats
import CGAT.Stats as Stats
import pysam

# only update R if called as pipeline
# otherwise - failure with sphinx
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# load options from the config file
import CGAT.Pipeline as P
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    on_error_raise=__name__ == '__main__')


PipelineDatabase.PARAMS = PARAMS
PipelineUCSC.PARAMS = PARAMS

###################################################################
###################################################################
# Helper functions mapping tracks to conditions, etc
###################################################################
import CGATPipelines.PipelineTracks as PipelineTracks

# need to be refactored

SEPARATOR = "|"

###################################################################
###################################################################
# Helper functions mapping tracks to conditions, etc
###################################################################
import CGATPipelines.PipelineTracks as PipelineTracks


class TracksVCF (PipelineTracks.Tracks):

    def load(self, filename, exclude=None):
        '''load tracks from a vcf file.'''
        tracks = []
        v = pysam.VCF()
        v.setversion(40)

        if not os.path.exists(filename):
            self.tracks = tracks
            return self
        v.connect(filename)

        if exclude:
            to_exclude = [re.compile(x) for x in exclude]

        for sample in v.getsamples():
            if exclude:
                for x in to_exclude:
                    if x.search(sample):
                        skip = True
                        break
                if skip:
                    continue

            tracks.append(self.factory(sample))

        self.tracks = tracks
        return self

TRACKS = TracksVCF(PipelineTracks.Sample).load("variants.vcf.gz")

###################################################################
###################################################################
###################################################################


def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh

# @transform( buildPileups, suffix(".pileup.gz"), ".pileup.stats")
# def countPileups( infile, outfile ):
#     '''get some basic counts from the pileup files.'''

#     to_cluster = True

#     statement = '''gunzip < %(infile)s
#     | python %(scriptsdir)s/snp2counts.py
#             --genome-file=genome
#             --module=contig-counts
#     > %(outfile)s
#     '''

#     P.run()

###################################################################
###################################################################
###################################################################
# Geneset
###################################################################
if "refseq_filename_gtf" in PARAMS:
    @split((PARAMS["refseq_filename_gtf"],
            PARAMS["refseq_filename_pep"],
            PARAMS["refseq_filename_cdna"],
            PARAMS["refseq_filename_map"],
            PARAMS["refseq_filename_ensembl"],
            ),
           (PARAMS["ensembl_filename_gtf"],
            PARAMS["ensembl_filename_pep"],
            PARAMS["ensembl_filename_cdna"]))
    def importRefseq(infiles, outfiles):
        '''convert a refseq gtf formatted file into an ensembl like
        gtf file.

        The refseq files should have been downloaded by USCS's
        table browser.

        Only unique refseq entries are used - all duplicates are
        removed. 

        This method imports the following files:
        gtf.gz, pep.fa.gz, cdna.fa.gz from the UCSC

        It also requires:

        * link.tsv.gz from the UCSC (table refLink)
           to add peptide identifiers and gene numbers.
        * ccdsinfo.tsv.gz from the UCSC (table ccdsInfo)
           to add a map from transcripts to ENSEMBL genes
        * refgene.tsv.gz from the UCSC (table refgene )
          to add gene_name. The refgene table contains
          most of the fields required for the gtf file,
          but unfortunately, the UCSC parser does
          not add it.
        '''

        infile_gtf, infile_pep, infile_cdna, infile_map, infile_ensembl = infiles
        outfile_gtf, outfile_pep, outfile_cdna = outfiles

        # build map between mrna and prot
        tmpfilename1 = P.getTempFilename()
        statement = '''gunzip < %(infile_map)s 
        | python %(scriptsdir)s/csv_cut.py mrnaAcc protAcc 
        | perl -p -e "s/\.\d+//g"
        > %(tmpfilename1)s
        '''
        P.run()

        # build map between mrna and gene - use ccds gene
        tmpfilename2 = P.getTempFilename()
        statement = '''gunzip < %(infile_map)s 
        | python %(scriptsdir)s/csv_cut.py mrnaAcc geneName
        | perl -p -e "s/\.\d+//g"
        > %(tmpfilename2)s
        '''
        P.run()

        statement = '''gunzip < %(infile_gtf)s
        | awk -v FS="\\t" -v OFS="\\t" '
              { $2 = "protein_coding"; print } '
        | python %(scriptsdir)s/gtf2gtf.py
           --method=remove-duplicates --duplicate-feature=ucsc
           --log=%(outfile_gtf)s.log
           --verbose=2
        | python %(scriptsdir)s/gtf2gtf.py
           --method=add-protein-id --map-tsv-file=%(tmpfilename1)s
           --log=%(outfile_gtf)s.log
           --verbose=2
        | python %(scriptsdir)s/gtf2gtf.py
           --method=rename-genes=gene
           --map-tsv-file=%(tmpfilename2)s
           --log=%(outfile_gtf)s.log
           --verbose=2
        | python %(scriptsdir)s/gtf2gtf.py
           --method=sort --sort-order=gene
        | gzip
        > %(outfile_gtf)s'''
        if not os.path.exists(outfile_gtf):
            P.run()

        for infile, outfile in ((infile_pep, outfile_pep),
                                (infile_cdna, outfile_cdna)):
            # remove numerical suffixes from identifiers
            statement = '''gunzip < %(infile)s
            | perl -p -e "s/\.\d+//g" 
            | gzip 
            > %(outfile)s'''
            if not os.path.exists(outfile):
                P.run()

        table = "ensembl2refseq"
        # use ENSEMBL mapping
        if 0:
            outf = open(tmpfilename1, "w")

            reader = CSV.DictReader(
                IOTools.openFile(infile_ensembl), dialect="excel-tab")
            c = E.Counter()

            outf.write(
                "gene_id\ttranscript_id\trefseq_transcript_id\trefseq_protein_id\tccds_id\n")

            for row in reader:
                c.input += 1
                gene_id, transcript_id, refseq_transcript_id, refseq_protein_id, ccds_id = \
                    [x.strip() for x in
                     (row["Ensembl Gene ID"],
                      row["Ensembl Transcript ID"],
                      row["RefSeq DNA ID"],
                      row["RefSeq Protein ID"],
                      row["CCDS ID"],
                      )]

                if not (transcript_id and gene_id and refseq_transcript_id and refseq_protein_id):
                    c.skipped += 1
                    continue

                c.output += 1
                outf.write("%s\t%s\t%s\t%s\t%s\n" %
                           (gene_id, transcript_id, refseq_transcript_id, refseq_protein_id, ccds_id))
            outf.close()

            statement = '''cat < %(tmpfilename1)s
            |python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
               --add-index=gene_id 
               --add-index=transcript_id 
               --add-index=refseq_transcript_id 
               --add-index=refseq_protein_id 
               --add-index=ccds_id 
               --table=%(table)s
            > refseq.load'''

            P.run()
            E.info("%s" % str(c))

        # use UCSC mapping
        statement = '''gunzip < %(infile_map)s
            | perl -p -i -e "s/\.\d+//g"
            | awk 'BEGIN {printf("ccds_id\\tsrc_db\\tttranscript_id\\tprotein_id\\n")} 
                   /^ccds/ {next} {print}'
            |python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
               --add-index=ccds_id 
               --add-index=transcript_id 
               --add-index=protein_id 
               --table=%(table)s
            > refseq.load'''

        P.run()

        os.unlink(tmpfilename1)
        os.unlink(tmpfilename2)

if "refseq_filename_gtf" in PARAMS:
    @split((PARAMS["refseq_filename_gtf"],
            PARAMS["refseq_filename_pep"],
            PARAMS["refseq_filename_cdna"],
            PARAMS["refseq_filename_map"],
            ),
           (PARAMS["ensembl_filename_gtf"],
            PARAMS["ensembl_filename_pep"],
            PARAMS["ensembl_filename_cdna"],
            "refseq.load"))
    def importRefseqFromUCSC(infiles, outfiles):
        '''convert a refseq gtf formatted file into an ensembl like
        gtf file.

        The refseq files should have been downloaded by USCS's
        table browser.

        Only unique refseq entries are used - all duplicates are
        removed. 

        This method imports the following files:
        gtf.gz, pep.fa.gz, cdna.fa.gz from the UCSC

        It also requires:

        * link.tsv.gz from the UCSC (table refLink)
           to add peptide identifiers and gene numbers.
        * ccdsinfo.tsv.gz from the UCSC (table ccdsInfo)
           to add a map from transcripts to ENSEMBL genes
        * refgene.tsv.gz from the UCSC (table refgene )
          to add gene_name. The refgene table contains
          most of the fields required for the gtf file,
          but unfortunately, the UCSC parser does
          not add it.
        '''

        infile_gtf, infile_pep, infile_cdna, infile_map = infiles
        outfile_gtf, outfile_pep, outfile_cdna, outfile_load = outfiles

        if not os.path.exists(outfile_gtf):
            PipelineUCSC.importRefSeqFromUCSC(
                infile_gtf, outfile_gtf, remove_duplicates=True)

        for infile, outfile in ((infile_pep, outfile_pep),
                                (infile_cdna, outfile_cdna)):
            # remove numerical suffixes from identifiers
            statement = '''gunzip < %(infile)s
            | perl -p -e "s/\.\d+//g" 
            | gzip 
            > %(outfile)s'''
            if not os.path.exists(outfile):
                P.run()

        # table = "ensembl2refseq"

        # use UCSC mapping
        # statement = '''gunzip < %(infile_map)s
        #     | perl -p -i -e "s/\.\d+//g"
        #     | awk 'BEGIN {printf("ccds_id\\tsrc_db\\tttranscript_id\\tprotein_id\\n")}
        #            /^ccds/ {next} {print}'
        #     |python %(scriptsdir)s/csv2db.py %(csv2db_options)s
        #        --add-index=ccds_id
        #        --add-index=transcript_id
        #        --add-index=protein_id
        #        --table=%(table)s
        #     > %(outfile_load)s'''
        # P.run()


@files("%s.fasta" % PARAMS["genome"], "%s.fa" % PARAMS["genome"])
def indexGenome(infile, outfile):
    '''index the genome for samtools.

    Samtools does not like long lines, so create a new file
    with split lines (what a waste).
    '''

    # statement = '''fold %(infile)s | perl -p -e "s/chr//" > %(outfile)s'''
    statement = '''fold %(infile)s > %(outfile)s'''
    P.run()

    pysam.faidx(outfile)

######################################################################
######################################################################
######################################################################


@files(((None, "pseudogenes.load"),))
def importPseudogenes(infile, outfile):
    '''import pseudogene data from pseudogenes.org'''

    tmpfile = "pseudogenes.tsv"

    # download
    if not os.path.exists(tmpfile + ".gz"):
        statement = '''
        wget -O %(tmpfile)s http://tables.pseudogene.org/dump.cgi?table=Mouse56;
        gzip %(tmpfile)s
        ''' % locals()

        P.run()

    tablename = P.toTable(outfile)

    statement = '''
    zcat %(tmpfile)s.gz
    | perl -p -i -e "s/Parent Protein/protein_id/; s/Chromosome/contig/; s/Start Coordinate/start/; s/Stop Coordiante/end/"
    |python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                     --table=%(tablename)s
                     --add-index=protein_id
    > %(outfile)s
    '''

    P.run()

######################################################################
######################################################################
######################################################################


@files(((None, "mgi.import"),))
def importMGI(infile, outfile):
    '''create via BIOMART'''

    filename = "mgi_biomart.tsv"

    if False:
        R.library("biomaRt")

        columns = {

            "marker_symbol_107": "marker_symbol",
            "marker_name_107": "marker_name",
            "mgi_allele_id_att": "allele_id",
            "allele_symbol_101": "allele_symbol",
            "allele_name_101": "allele_name",
            "allele_type_101": "allele_type",
            "phenotype_id_106_att": "phenotype_id",
            "ensembl_gene_id_103": "gene_id"}

        keys = columns.keys()

        mgi = R.useMart(biomart="biomart", dataset="markers")
        result = R.getBM(attributes=keys, mart=mgi)

        outf = open(filename, "w")
        outf.write("\t".join([columns[x] for x in keys]) + "\n")

        for data in zip(*[result[x] for x in keys]):
            outf.write("\t".join(map(str, data)) + "\n")

        outf.close()

    if not os.path.exists(filename):

        # associations need to be downloaded individually

        R.library("biomaRt")

        columns = {
            "mgi_marker_id_att": "marker_id",
            "marker_name_107": "marker_name",
            "mgi_allele_id_att": "allele_id",
            "allele_symbol_101": "allele_symbol",
            "allele_name_101": "allele_name",
            "allele_type_101": "allele_type",
            "phenotype_id_106_att": "phenotype_id",
            "ensembl_gene_id_103": "gene_id"}

        def downloadData(filename, columns):
            '''download data via biomart into filename.
               translate column headers.'''

            if os.path.exists(filename):
                return

            E.info("downloading data for %s" % filename)

            keys = columns.keys()

            mgi = R.useMart(biomart="biomart", dataset="markers")
            result = R.getBM(attributes=keys,
                             mart=mgi)

            if len(result.rx(keys[0])) == 0:
                raise ValueError(
                    "no data for %s: using keys=%s" % (filename, keys))

            outf = open(filename, "w")
            outf.write("\t".join([columns[x] for x in keys]) + "\n")

            for data in zip(*[result.rx(x) for x in keys]):
                outf.write("\t".join(map(str, data)) + "\n")

            outf.close()

        downloadData("mgi_marker2allele.tsv",
                     {"mgi_marker_id_att": "marker_id",
                      "mgi_allele_id_att": "allele_id"})

        downloadData("mgi_allele2phenotype.tsv",
                     {"mgi_allele_id_att": "allele_id",
                      "phenotype_id_106_att": "phenotype_id"})

        downloadData("mgi_marker2gene.tsv",
                     {"mgi_marker_id_att": "marker_id",
                      "ensembl_gene_id_103": "gene_id"})

        downloadData("mgi_markers.tsv",
                     {"mgi_marker_id_att": "marker_id",
                      "marker_symbol_107": "symbol",
                      "marker_name_107": "name",
                      "marker_type_107": "type",
                      })

        downloadData("mgi_alleles.tsv",
                     {"mgi_allele_id_att": "allele_id",
                      "allele_symbol_101": "symbol",
                      "allele_name_101": "name",
                      "allele_type_101": "type"
                      })

        downloadData("mgi_phenotypes.tsv",
                     {"phenotype_id_106_att": "phenotype_id",
                      "term_106": "term"})

        for filename in glob.glob("mgi_*.tsv"):
            tablename = filename[:-len(".tsv")]
            E.info("loading %s" % tablename)

            # remove duplicate rows
            # remove rows where only the first field is set
            statement = '''
            perl -p -e "s/\\s+\\n/\\n/" < %(filename)s
            | %(scriptsdir)s/hsort 1
            | uniq
            | awk '{ for (x=2; x<=NF; x++) { if ($x != "") { print; break;} } }'
            |python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                     --table=%(tablename)s
                     --add-index=marker_id
                     --add-index=allele_id
                     --add-index=gene_id
                     --map=allele_name:str
                     --map=symbol:str
                     --add-index=phenotype_id
           >> %(outfile)s
           '''

            P.run()


@merge(None, "gene2omim.load")
def loadGene2Omim(infile, outfile):
    '''download gene id - OMIM associations via BIOMART.

    Missing numerical entries are set to 0 (from -2147483648). 
    '''

    tablename = P.toTable(outfile)

    columns = {
        "ensembl_gene_id": "gene_id",
        "mim_gene_accession": "mim_gene_id",
        "mim_morbid_accession": "mim_morbid_id",
        "mim_morbid_description": "mim_morbid_description",
    }

    data = PBiomart.biomart_iterator(
        columns.keys(), biomart="ensembl", dataset="hsapiens_gene_ensembl")

    def transform_data(data):
        for result in data:
            for x in ("mim_gene_accession", "mim_morbid_accession"):
                result[x] = ("", result[x])[result[x] >= 0]
            result["mim_morbid_description"] = result[
                "mim_morbid_description"].strip()
            yield result

    PipelineDatabase.importFromIterator(
        outfile, tablename, transform_data(data),
        columns=columns, indices=("gene_id", ))


@merge(None, "orthologs.load")
def loadHumanOrthologs(infile, outfile):
    '''download human2mouse orthologs
    '''

    tablename = P.toTable(outfile)

    if "orthology_species" in PARAMS:
        # we are within human, create a set of default 1:1 orthologs
        species = PARAMS["orthology_species"]
        columns = {
            "ensembl_gene_id": "hs_gene_id",
            "%s_ensembl_gene" % species: "gene_id",
            "%s_orthology_type" % species: "orthology_type",
            "%s_homolog_ds" % species: "ds",
        }

        data = PBiomart.biomart_iterator(
            columns.keys(), biomart="ensembl", dataset="hsapiens_gene_ensembl")

        PipelineDatabase.importFromIterator(
            outfile, tablename, data,
            columns=columns, indices=("hs_gene_id", "gene_id", ))

    else:
        # we are within human, create a set of default 1:1 orthologs
        columns = {
            "ensembl_gene_id": "hs_gene_id",
            "gene_id": "gene_id",
            "orthology_type": "orthology_type",
            "ds": "ds",
        }

        dbh = connect()

        statement = '''CREATE TABLE %(tablename)s AS 
                       SELECT gene_id AS hs_gene_id,
                              gene_id AS gene_id,
                              "ortholog_one2one" AS orthology_type,
                              0 AS ds 
                       FROM annotations.gene_info'''

        Database.executewait(
            dbh, "DROP TABLE IF EXISTS %(tablename)s" % locals())
        Database.executewait(dbh, statement % locals())
        Database.executewait(
            dbh, "CREATE INDEX %(tablename)s_index1 ON %(tablename)s (hs_gene_id)" % locals())
        Database.executewait(
            dbh, "CREATE INDEX %(tablename)s_index2 ON %(tablename)s (hs_gene_id)" % locals())

    P.touch(outfile)

###################################################################
###################################################################
###################################################################
# MAIN PIPELINE
###################################################################
###################################################################
###################################################################

############################################################
############################################################
############################################################


@merge(None, "seleno.list")
def buildSelenoList(infile, outfile):
    '''export a list of seleno cysteine transcripts.'''

    dbhandle = connect()
    cc = dbhandle.cursor()
    statement = '''
    SELECT DISTINCT transcript_id
    FROM annotations.transcript_info as t,
         annotations.protein_stats as p
    WHERE p.protein_id = t.protein_id AND
         p.nU > 0
    '''
    outf = open(outfile, "w")
    outf.write("transcript_id\n")
    outf.write("\n".join([x[0] for x in cc.execute(statement)]) + "\n")
    outf.close()

###################################################################
###################################################################
###################################################################
# Targets for prepare
###################################################################

###################################################################
###################################################################
###################################################################


@files(((PARAMS["ensembl_ensembl2uniprot"], "ensembl2uniprot.load"), ))
def loadEnsembl2Uniprot(infile, outfile):
    '''load mapping from ENSEMBL transcripts ids to
    uniprot ids.

    This method expects an BioMart output file with the following 
    five columns: 
    Ensembl gene id, 
    Ensembl transcript id, 
    Uniprot Swissprot Id,
    Uniprot Accession
    Uniport/Trembl Accession
    '''

    table = P.toTable(outfile)

    statement = '''gunzip 
    < %(infile)s
    | perl -p -e 
          "s/Ensembl Gene ID/gene_id/; 
           s/Ensembl Transcript ID/transcript_id/; 
           s/UniProt\/SwissProt ID/swissprot_id/;       
           s/UniProt\/SwissProt Accession/swissprot_acc/;
           s/UniProt\/TrEMBL Accession/trembl_acc/"
    |python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
        --add-index=gene_id \
        --add-index=transcript_id \
        --add-index=trembl_acc \
        --table=%(table)s
    > %(outfile)s
    '''
    P.run()


###################################################################
###################################################################
###################################################################
@files(os.path.join(PARAMS["annotations_dir"],
                    PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
       "annotations_bases.fasta")
def buildBaseAnnotations(infile, outfile):
    """build base annotations"""

    to_cluster = True

    dbname = outfile[:-len(".fasta")]
    statement = """
        gunzip < %(infile)s 
        | python %(scriptsdir)s/gtf2fasta.py 
                --force-output 
                --genome=%(genome_dir)s/%(genome)s 
                --output-filename-pattern=annotations_bases.%%s 
                --log=%(outfile)s.log 
        | python %(toolsdir)s/index_fasta.py 
        --log=%(outfile)s.log 
        %(dbname)s - > %(outfile)s.log
    """

    P.run()

###################################################################
###################################################################
###################################################################


@files(os.path.join(PARAMS["annotations_dir"],
                    PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
       "annotations_exons.gtf")
def buildExonAnnotations(infile, outfile):
    """build exon annotations"""

    to_cluster = True

    statement = """
        gunzip < %(infile)s 
        | awk '$3 == "CDS"' 
        | python %(scriptsdir)s/gff2gff.py --method=sanitize=genome --skip-missing 
               --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gtf2gff.py 
                --method=exons 
                --restrict-source=protein_coding 
                --log=%(outfile)s.log 
        > %(outfile)s
    """

    P.run()

###################################################################
###################################################################
###################################################################


@files(os.path.join(PARAMS["annotations_dir"],
                    PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
       "annotations_genes.gtf")
def buildGeneAnnotations(infile, outfile):
    """build gene annotations.

    Merge exons per gene within the reference set. The
    output includes the UTR and non-coding genes.
    """
    statement = """
        gunzip < %(infile)s |\
        python %(scriptsdir)s/gtf2gtf.py --method=merge-exons --with-utr --log=%(outfile)s.log |\
        python %(scriptsdir)s/gtf2gtf.py --method=set-transcript-to-gene --log=%(outfile)s.log |\
        python %(scriptsdir)s/gff2gff.py --skip-missing --method=sanitize=genome --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log |\
        %(scriptsdir)s/gff_sort gene-pos \
        > %(outfile)s
    """
    queue = "server"
    P.run()

###################################################################
###################################################################
###################################################################


@files(buildGeneAnnotations, "annotations_genes.counts")
def makeGeneCounts(infile, outfile):
    """coun gene exon statistics.
    """

    statement = """
    cat < %(infile)s |\
    python %(scriptsdir)s/gtf2table.py \
        --genome-file=genome \
        --counter=length \
        --log=%(outfile)s.log \
    > %(outfile)s
    """
    P.run()

###################################################################
###################################################################
###################################################################


@follows(buildBaseAnnotations, buildExonAnnotations)
@files([("variants.vcf.gz", "%s.annotations.gz" % x, x) for x in TRACKS])
def buildAnnotations(infile, outfile, sample):
    """annotate snps with gene set."""

    to_cluster = True

    bases = "annotations_bases"

    statement = """
    python %(scriptsdir)s/snp2table.py 
        --input-format=vcf
        --vcf-sample=%(sample)s
        --vcf-file=%(infile)s
        --genome-file=%(genome_dir)s/%(genome)s 
        --annotations-tsv-file=%(bases)s 
        --log=%(outfile)s.log 
    | gzip > %(outfile)s
    """
    P.run()

###################################################################
###################################################################
###################################################################


@transform(buildAnnotations,
           suffix('.annotations.gz'),
           '_annotations.load')
def loadAnnotations(infile, outfile):
    '''load annotations'''

    tablename = P.toTable(outfile)

    statement = '''gunzip 
    < %(infile)s
    |python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              --quick
              --map=gene_id:str 
              --add-index=gene_id 
              --table=%(tablename)s
              --map=base_qualities:text 
    > %(outfile)s
    '''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(buildAnnotations,
           suffix('.annotations.gz'),
           '_annotations.summary')
def summarizeAnnotations(infile, outfile):
    '''compute summary stats for annotation files.'''

    to_cluster = True

    # count substitutions for each category
    statement = '''gunzip 
    < %(infile)s
    | python %(scriptsdir)s/csv_cut.py code reference_base genotype variant_type 
    | awk '$4 == "variant_type" { printf("%%s-%%s-%%s\\tcounts\\n", $1,$2,$3); } 
           $4 == "E" || $4 == "O" {printf("%%s-%%s-%%s\\t1\\n", $1,$2,$3)}'
    | sort 
    | uniq -c 
    | awk 'BEGIN{ printf("code-reference_base-genotype\\tcounts\\n" ); } \
                  $2 !~ /^code/ {printf("%%s\\t%%i\\n",$2,$1);}'
    | perl -p -i -e "s/-/\\t/g unless (/^#/)"
    > %(outfile)s
    '''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(summarizeAnnotations,
           suffix('_annotations.summary'),
           '_annotations_summary.load')
def loadAnnotationsSummary(infile, outfile):
    '''load annotations'''

    tablename = P.toTable(outfile)

    statement = '''cat
    < %(infile)s
    |python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              --add-index=code
              --table=%(tablename)s
    > %(outfile)s
    '''
    P.run()

###################################################################
###################################################################
###################################################################


@follows(buildSelenoList)
@files([("variants.vcf.gz", "%s.effects.gz" % x, x) for x in TRACKS])
def buildEffects(infile, outfile, sample):
    """annotate snps with gene set."""

    to_cluster = True

    seleno = "seleno.list"

    transcripts = os.path.join(PARAMS["annotations_dir"],
                               PARAMS_ANNOTATIONS["interface_geneset_cds_gtf"])

    statement = """
    python %(scriptsdir)s/snp2counts.py 
        --genome-file=%(genome_dir)s/%(genome)s
        --input-format=vcf
        --vcf-sample=%(sample)s
        --module=transcript-effects 
        --seleno-tsv-file=%(seleno)s 
        --vcf-file=%(infile)s
        --exons-file=%(transcripts)s 
        --output-filename-pattern=%(outfile)s.%%s.gz
        --log=%(outfile)s.log 
    | gzip 
    > %(outfile)s
    """
    P.run()

###################################################################
###################################################################
###################################################################


@transform(buildEffects,
           suffix(".effects.gz"),
           "_effects.load")
def loadEffects(infile, outfile):
    '''load transcript effects into tables.'''

    root = infile[:-len(".effects.gz")]

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --from-zipped \
              --add-index=transcript_id \
              --table=%(root)s_effects \
    < %(infile)s > %(outfile)s
    '''
    P.run()

    for suffix in ("cds", "intron", "splicing", "translation"):

        statement = '''
        gunzip 
        < %(infile)s.%(suffix)s.gz
        |python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
        --allow-empty-file
        --add-index=transcript_id 
        --table=%(root)s_effects_%(suffix)s 
        --ignore-column=seq_na
        --ignore-column=seq_aa
        >> %(outfile)s
        '''
        P.run()

###################################################################
###################################################################
###################################################################


@merge(buildEffects, "effects.load")
def mergeEffects(infiles, outfile):
    '''load transcript effects into single table.'''

    tablename = P.toTable(outfile)
    outf = IOTools.openFile('effects.tsv.gz', 'w')
    first = True
    for f in infiles:
        track = P.snip(os.path.basename(f), ".effects.gz")
        if not os.path.exists(f):
            E.warn("File %s missing" % f)
            continue
        lines = [x for x in gzip.open(f, "r").readlines()]
        if first:
            outf.write("%s\t%s" % ("track", lines[0]))
        first = False
        for i in range(1, len(lines)):
            outf.write("%s\t%s" % (track, lines[i]))
    outf.close()
    tmpfilename = outf.name

    statement = '''zcat %(tmpfilename)s 
                   | python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                       --add-index=transcript_id 
                       --table=%(tablename)s 
                   > %(outfile)s'''
    P.run()

    for suffix in ("cds", "intron", "splicing", "translation"):

        outf = IOTools.openFile('effects_' + suffix + '.tsv.gz', 'w')
        first = True
        for f in infiles:
            track = P.snip(os.path.basename(f), ".effects.gz")
            statfile = P.snip(f, ".gz") + suffix + ".gz"
            if not os.path.exists(statfile):
                E.warn("File %s missing" % statfile)
                continue
            lines = [x for x in gzip.open(statfile, "r").readlines()]
            if first:
                outf.write("%s\t%s" % ("track", lines[0]))
            first = False
            for i in range(1, len(lines)):
                outf.write("%s\t%s" % (track, lines[i]))
        outf.close()
        tmpfilename = outf.name

        statement = '''zcat %(tmpfilename)s 
                       | python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                           --allow-empty-file
                           --add-index=transcript_id 
                           --table=%(tablename)s_%(suffix)s 
                           --ignore-column=seq_na
                           --ignore-column=seq_aa
                       >> %(outfile)s'''
        P.run()

###################################################################
###################################################################
###################################################################


@follows(buildSelenoList)
@files([("variants.vcf.gz", "%s.alleles" % x, x) for x in TRACKS])
def buildAlleles(infile, outfile, sample):
    """annotate snps with gene set."""

    to_cluster = True

    seleno = "seleno.list"

    transcripts = os.path.join(PARAMS["annotations_dir"],
                               PARAMS_ANNOTATIONS["interface_geneset_cds_gtf"])

    statement = """zcat %(transcripts)s 
    | python %(scriptsdir)s/gtf2alleles.py 
        --genome-file=%(genome_dir)s/%(genome)s
        --seleno-tsv-file=%(seleno)s 
        --output-filename-pattern=%(outfile)s.%%s.gz
        --vcf-file=%(infile)s
        --vcf-sample=%(sample)s
    > %(outfile)s
    """
    P.run()

###################################################################
###################################################################
###################################################################


@transform(buildAlleles,
           suffix(".alleles"),
           "_alleles.load")
def loadAlleles(infile, outfile):
    '''load allele.'''

    tablename = outfile[:-len(".load")]

    statement = '''gunzip
    < %(infile)s.table.gz
    | perl -p -e "s/False/0/g; s/True/1/g;"
    |python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              --add-index=gene_id 
              --add-index=transcript_id 
              --ignore-column=cds
              --ignore-column=peptide
              --table=%(tablename)s 
    > %(outfile)s
    '''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(loadAlleles,
           suffix("_alleles.load"),
           "_alleles_transcripts.load")
def summarizeAllelesPerTranscript(infile, outfile):
    '''summarize effects on a per-gene level.

    The following fields are exclusive:

    is_wildtype
       both alleles are wildtype

    is_knockout
       both alleles knocked out

    is_truncated
       both alleles truncated or truncated and knocked out

    is_affected
       one allele is truncated or knocked out

    The other fields are not necessarily exclusive, for example there
    could be transcripts with one knocked out allele and one wildtype
    allele, such that ``is_nmd_affected``, ``is_affected`` and ``has_wildtype`` 
    are all true.
    '''

    tablename = outfile[:-len(".load")]
    track = infile[:-len("_alleles.load")]

    dbhandle = connect()

    statement = '''
    CREATE TABLE %(tablename)s AS
    SELECT DISTINCT 
           transcript_id,
           COUNT(DISTINCT allele_id) AS nalleles,
           CASE WHEN SUM( is_nmd_knockout) = 2 THEN 1 ELSE 0 END AS is_nmd_knockout,
           CASE WHEN SUM( is_nmd_knockout) >= 1 THEN 1 ELSE 0 END AS is_nmd_affected,
           CASE WHEN SUM( is_splice_truncated) = 2 THEN 1 ELSE 0 END AS is_splice_truncated,
           CASE WHEN SUM( is_splice_truncated) >= 1 THEN 1 ELSE 0 END AS is_splice_affected,
           CASE WHEN SUM( is_stop_truncated) = 2 THEN 1 ELSE 0 END AS is_stop_truncated,
           CASE WHEN SUM( is_stop_truncated) >= 1 THEN 1 ELSE 0 END AS is_stop_affected,
           CASE WHEN SUM( is_wildtype ) = 2 THEN 1 ELSE 0 END AS is_wildtype, 
           CASE WHEN SUM( is_wildtype ) >= 1 THEN 1 ELSE 0 END AS has_wildtype, 
           contig AS contig, 
           strand AS strand, 
           GROUP_CONCAT( reference_first_stop_start ) AS stop_codons_start,
           GROUP_CONCAT( reference_first_stop_end ) AS stop_codons_end,
           0 AS is_knockout,
           0 AS is_truncated,
           0 AS is_affected
    FROM %(track)s_alleles AS a
    GROUP BY transcript_id
    ''' % locals()

    Database.executewait(
        dbhandle, "DROP TABLE IF EXISTS %(tablename)s" % locals())
    Database.executewait(dbhandle, statement)
    Database.executewait(
        dbhandle, "CREATE INDEX %(tablename)s_transcript_id ON %(tablename)s (transcript_id)" % locals())
    Database.executewait(
        dbhandle, "UPDATE %(tablename)s SET is_knockout = is_nmd_knockout" % locals())
    Database.executewait( dbhandle, '''UPDATE %(tablename)s SET is_truncated = 
                                       is_splice_truncated OR is_stop_truncated OR 
                                       (is_splice_affected AND is_stop_affected) OR 
                                       (is_splice_affected AND is_nmd_affected) OR 
                                       (is_stop_affected AND is_nmd_affected)
                                       ''' % locals())
    Database.executewait(dbhandle, 'UPDATE %(tablename)s SET is_affected ='
                         '(is_nmd_affected OR is_splice_affected OR is_stop_affected) AND NOT'
                         '(is_knockout or is_truncated)' % locals())
    dbhandle.commit()

    P.touch(outfile)

###################################################################
###################################################################
###################################################################


@transform(summarizeAllelesPerTranscript,
           suffix("_alleles_transcripts.load"),
           "_alleles_genes.load")
def summarizeAllelesPerGene(infile, outfile):
    '''summarize effects on a per-gene level.'''

    tablename = outfile[:-len(".load")]
    track = infile[:-len(".load")]

    dbhandle = connect()

    statement = '''
    CREATE TABLE %(tablename)s AS
    SELECT DISTINCT 
           i.gene_id AS gene_id,
           COUNT( DISTINCT a.transcript_id ) AS ntranscripts,
           CASE WHEN SUM( is_nmd_knockout ) = COUNT(DISTINCT a.transcript_id) THEN 1 ELSE 0 END AS is_nmd_knockout,
           SUM( is_nmd_knockout ) AS is_nmd_affected,
           CASE WHEN SUM( is_splice_truncated) = COUNT(DISTINCT a.transcript_id) THEN 1 ELSE 0 END AS is_splice_truncated,
           SUM( is_splice_truncated ) AS is_splice_affected,
           CASE WHEN SUM( is_stop_truncated ) = COUNT(DISTINCT a.transcript_id) THEN 1 ELSE 0 END AS is_stop_truncated,
           SUM( is_stop_truncated ) AS is_stop_affected,
           CASE WHEN SUM( is_wildtype ) = COUNT(DISTINCT a.transcript_id) THEN 1 ELSE 0 END AS is_wildtype, 
           SUM( is_wildtype ) AS has_wildtype, 
           contig AS contig, 
           strand AS strand, 
           GROUP_CONCAT( stop_codons_start ) AS stop_codons_start,
           GROUP_CONCAT( stop_codons_end ) AS stop_codons_end,
           0 AS is_knockout,
           0 AS is_truncated,
           0 AS is_affected
    FROM %(track)s AS a, annotations.transcript_info AS i
    WHERE i.transcript_id = a.transcript_id
    GROUP BY i.gene_id
    ''' % locals()

    Database.executewait(
        dbhandle, "DROP TABLE IF EXISTS %(tablename)s" % locals())
    Database.executewait(dbhandle, statement)
    Database.executewait(
        dbhandle, "CREATE INDEX %(tablename)s_gene_id ON %(tablename)s (gene_id)" % locals())
    Database.executewait(
        dbhandle, "UPDATE %(tablename)s SET is_knockout = is_nmd_knockout" % locals())
    Database.executewait( dbhandle, '''UPDATE %(tablename)s SET is_truncated = 
                                       is_splice_truncated OR is_stop_truncated OR 
                                       (is_splice_affected + is_stop_affected >= ntranscripts)  
                          ''' % locals())
    Database.executewait(dbhandle, 'UPDATE %(tablename)s SET is_affected ='
                         '(is_nmd_affected OR is_splice_affected OR is_stop_affected) AND NOT'
                         '(is_knockout or is_truncated)' % locals())

    dbhandle.commit()

    P.touch(outfile)


@merge(summarizeAllelesPerGene,
       "summary_alleles_genes.load")
def combineSummaryAllelesPerGene(infiles, outfile):

    dbhandle = connect()

    tracks = [x[:-len("_alleles_genes.load")] for x in infiles]

    tablename_prefix = P.toTable(outfile)

    fields = ", ".join(["%s INT" % x for x in tracks])

    statement_create = '''
    CREATE TABLE %(tablename)s 
           ( gene_id TEXT,
             total INT,
             %(fields)s )'''

    statement_insert = '''
    INSERT INTO %(tablename)s
             VALUES( '%(gene_id)s', %(total)i, %(matrix)s )
    '''

    statement_allgenes = "SELECT DISTINCT gene_id FROM gene_info"

    for field in ("is_knockout", "is_truncated"):

        tablename = "%s_%s" % (tablename_prefix, field)
        E.info("creating %s" % tablename)

        all_genes = dict([(x[0], set())
                          for x in Database.executewait(dbhandle, statement_allgenes % locals())])

        Database.executewait(
            dbhandle, "DROP TABLE IF EXISTS %(tablename)s" % locals())
        Database.executewait(dbhandle, statement_create % locals())

        for track in tracks:
            statement = """SELECT gene_id 
                           FROM %(track)s_alleles_genes
                           WHERE %(field)s"""

            genes = [x[0]
                     for x in Database.executewait(dbhandle, statement % locals())]

            for gene in genes:
                all_genes[gene].add(track)

        for gene_id, data in all_genes.iteritems():
            matrix = [0] * len(tracks)
            for x, track in enumerate(tracks):
                if track in data:
                    matrix[x] = 1
            total = sum(matrix)
            matrix = ",".join([str(x) for x in matrix])
            Database.executewait(dbhandle, statement_insert % locals())

        Database.executewait(dbhandle,
                             "CREATE INDEX %(tablename)s_index1 on %(tablename)s (gene_id)" % locals())

    P.touch(outfile)

###################################################################
###################################################################
###################################################################


@transform(loadEffects,
           suffix("_effects.load"),
           "_effects_genes.load")
def summarizeEffectsPerGene(infile, outfile):
    '''summarize effects on a per-gene level.'''

    tablename = outfile[:-len(".load")]
    track = infile[:-len("_effects.load")]

    dbhandle = connect()

    statement = '''
    CREATE TABLE %(tablename)s AS
    SELECT DISTINCT 
           gene_id, 
           COUNT(*) AS ntranscripts,
           MIN(e.nalleles) AS min_nalleles,
           MAX(e.nalleles) AS max_nalleles,
           MIN(e.stop_min) AS min_stop_min,
           MAX(e.stop_min) AS max_stop_min,
           MIN(e.stop_max) AS min_stop_max,
           MAX(e.stop_max) AS max_stop_max,
           SUM( CASE WHEN stop_min > 0 AND cds_len - stop_min * 3 < last_exon_start THEN 1  
                     ELSE 0 END) AS nmd_knockout,
           SUM( CASE WHEN stop_max > 0 AND cds_len - stop_max * 3 < last_exon_start THEN 1  
                     ELSE 0 END) AS nmd_affected
    FROM annotations.transcript_info as i,
         %(track)s_effects AS e
    WHERE i.transcript_id = e.transcript_id
    GROUP BY i.gene_id
    ''' % locals()

    Database.executewait(
        dbhandle, "DROP TABLE IF EXISTS %(tablename)s" % locals())
    Database.executewait(dbhandle, statement)
    Database.executewait(
        dbhandle, "CREATE INDEX %(tablename)s_gene_id ON %(tablename)s (gene_id)" % locals())
    dbhandle.commit()

    P.touch(outfile)


###################################################################
@follows(buildGeneAnnotations)
@files_re(glob.glob('*.pileup.gz'),
          '(.*).pileup.gz',
          [r'\1.pileup.gz', "annotations_genes.gtf"],
          r'\1.genecounts.gz')
def makeSNPCountsPerGene(infiles, outfile):
    """count snps within genes"""

    infile_snps, infile_genes = infiles

    statement = """
    gunzip < %(infile_snps)s |\
    grep -v "^NT" |\
    python %(scriptsdir)s/snp2counts.py \
        --genome-file=genome \
        --exons-file=%(ensembl_filename_gtf)s \
        --log=%(outfile)s.log |\
    gzip > %(outfile)s
    """
    P.run()

############################################################################


@follows(mkdir(os.path.join(PARAMS["scratchdir"], "malis.dir")))
@merge(buildAlleles, "malis.map")
def setupMultipleAlignment(infiles, outfile):
    '''prepare input files for multiple alignment computations.

    This script does some id-mapping to resolve coordinates.

    Basically, each genome is separated into two alleles. 
    Gene_id's will be suffixed with the allele_id. This ensures
    that exons of a gene with multiple transcipts will be resolved
    correctly with consistent coordinates. 

    From an alignment point of view, the two alleles of the genes will be treated
    independently, but transcripts within a gene  will be merged correctly at exon 
    boundaries, again on a per-allele basis.

    Later, when collecting the results, the allele id is moved from the gene to
    the transcript.
    '''

    targetdir = os.path.join(PARAMS["scratchdir"], "malis.dir")

    filepool_gtf = IOTools.FilePoolMemory(
        "%(targetdir)s/cluster_%%s.dir/cluster_%%s.gtf" % locals())
    filepool_pep = IOTools.FilePoolMemory(
        "%(targetdir)s/cluster_%%s.dir/cluster_%%s_pep.fasta" % locals())
    filepool_cds = IOTools.FilePoolMemory(
        "%(targetdir)s/cluster_%%s.dir/cluster_%%s_cds.fasta" % locals())

    outf = open(outfile, "w")
    outf.write("id\tgroup_id\n")

    map_gene2group = {}
    map_seqid2code = {}
    x = 0
    counts = E.Counter()
    for infile in infiles:
        track = infile[:-len(".alleles")]
        E.info("adding track %s" % track)

        reader = CSV.DictReader(
            open(infile + ".table", "rU"), dialect="excel-tab")
        for row in reader:
            counts.input += 1
            gene_id, allele_id, transcript_id = row[
                "gene_id"], row["allele_id"], row["transcript_id"]
            if gene_id not in map_gene2group:
                map_gene2group[gene_id] = len(map_gene2group)
            group_id = map_gene2group[gene_id]
            new_gene_id = "-".join((gene_id, allele_id))
            if row["is_wildtype"] == "1":
                code = "WT"
            if row["is_nmd_knockout"] == "1":
                counts.nmd_knockouts += 1
                continue

            else:
                code = "VA"
            seq_id = SEPARATOR.join((track, transcript_id, new_gene_id))
            map_seqid2code[seq_id] = code
            seq_id = SEPARATOR.join((seq_id, code))
            outf.write("%s\t%i\n" % (seq_id, group_id))
            filepool_pep.write(str(group_id), ">%s\n%s\n" %
                               (seq_id, row["peptide"]))
            filepool_cds.write(str(group_id), ">%s\n%s\n" %
                               (seq_id, row["cds"]))
            counts.written += 1

        with open(infile + ".gtf") as inf:
            for gtf in GTF.iterator(inf):
                group_id = map_gene2group[gtf.gene_id]
                new_gene_id = "-".join((gtf.gene_id, gtf["allele_id"]))
                seq_id = SEPARATOR.join(
                    (track, gtf.transcript_id, new_gene_id))
                seq_id = SEPARATOR.join((seq_id, map_seqid2code[seq_id]))
                gtf.transcript_id = seq_id
                filepool_gtf.write(group_id, str(gtf) + "\n")

        x += 1
        # if x > 2: break

    E.info("writing data")
    filepool_gtf.close()
    filepool_pep.close()
    filepool_cds.close()
    outf.close()
    counts.ngroups = len(map_gene2group)
    counts.nsequences = len(map_seqid2code)

    E.info("%s\n" % (str(counts)))


@transform(os.path.join(PARAMS["scratchdir"], "malis.dir", "*", "*.gtf"),
           suffix(".gtf"),
           ".mali")
def buildMultipleAlignments(infile, outfile):
    '''build multiple alignments.'''

    track = infile[:-len(".gtf")]
    filename_cds = track + "_cds.fasta"
    filename_pep = track + "_pep.fasta"

    to_cluster = True

    statement = '''
	python %(scriptsdir)s/align_transcripts.py \
		--gtf-file=%(infile)s \
		--cds-gtf-file=%(filename_cds)s \
		--force-map \
		--verbose=2 \
		--output-filename-pattern=%(track)s_%%s.fasta \
		--output=final_aa \
		--output=final_na \
		--output=aligned_aa \
		--output=aligned_na \
		--output-format="plain-fasta" \
	< %(filename_pep)s > %(outfile)s
      '''

    P.run()


@merge(buildMultipleAlignments, "variants")
def buildMultipleAlignmentVariantColumns(infile, outfile):
    '''build multiple alignments.'''

    track = infile[:-len(".gtf")]
    filename_cds = track + "_cds.fasta"
    filename_pep = track + "_pep.fasta"

    to_cluster = True

    statement = '''
	python %(scriptsdir)s/malis2mali.py \
		--gtf-file=%(infile)s \
		--cds-gtf-file=%(filename_cds)s \
		--force-map \
		--verbose=2 \
		--output-filename-pattern=%(track)s_%%s.fasta \
		--output=final_aa \
		--output=final_na \
		--output=aligned_aa \
		--output=aligned_na \
		--output-format="plain-fasta" \
	< %(filename_pep)s > %(outfile)s
      '''

    P.run()


@merge(buildMultipleAlignments, "malis.result")
def mergeMultipleAlignments(infiles, outfile):
    '''collect multiple alignment results into files that
    are compatible with OPTIC.
    '''

    for section in ("final_aa", "final_na", "aligned_aa", "aligned_na"):
        outfilename = outfile + "." + section + ".gz"

        counter = E.Counter()

        E.info("processing %s into %s" % (section, outfilename))
        outf = gzip.open(outfilename, "w")
        outf.write(
            "cluster_id\tspecies\ttranscript_id\tgene_id\tcode\tsequence\n")
        for infile in infiles:
            counter.input += 1
            dirname, filename = os.path.split(infile)
            cluster_id = re.match("cluster_(\d+).mali", filename).groups()[0]
            infilename = os.path.join(
                dirname, "cluster_%s_%s.fasta" % (cluster_id, section))
            # E.debug( "adding %s - %s from %s" % (filename, cluster_id, infilename) )
            if not os.path.exists(infilename):
                counter.missing += 1
                E.warn("multiple alignment %s missing" % infilename)
                continue
            for entry in FastaIterator.FastaIterator(open(infilename, "r")):
                parts = entry.title.split(SEPARATOR)
                if len(parts) == 4:
                    species, transcript_id, gene_id, code = entry.title.split(
                        SEPARATOR)
                elif len(parts) == 2:
                    species, gene_id = entry.title.split(SEPARATOR)
                    transcipt_id = gene_id
                    code = "CG"
                # transfer the allele_id from the gene to the transcript
                gene_id, allele_id = gene_id.split("-")
                transcript_id += "-" + allele_id

                outf.write("\t".join(map(str,
                                         (cluster_id,
                                          species,
                                          transcript_id,
                                          gene_id,
                                          code,
                                          entry.sequence))) + "\n")
            counter.output += 1

        outf.close()
        E.info("%s: %s" % (outfilename, str(counter)))

    P.touch(outfile)


@merge('*_pileup.load',
       "genome.maf.gz")
def buildMAF(infiles, outfile):

    tracks = " ".join(["--track=%s" % x[:-len(".load")] for x in infiles])

    statement = '''
    gunzip 
    < transcripts.gtf.gz 
    | python %(scriptsdir)s/gtf2gtf.py
           --method=merge-transcripts --with-utr 
    | %(cmd-farm)s --split-at-lines=100 --log=%(outfile)s.log --is-binary -v 10 
    "python %(scriptsdir)s/snp2maf.py 
          --genome=genome 
          %(tracks)s 
          --reference=mm9 
          --is-gtf 
          --pattern='\(\\\\\\S+\)_pileup'
          --log=%(outfile)s.log" | gzip
    > %(outfile)s
    '''

    P.run()

###################################################################
###################################################################
###################################################################


@merge(summarizeAllelesPerGene,
       ("export/nmd_knockouts.tsv.gz",
        "export/nmd_knockouts_summary.tsv.gz",
        ))
def exportKnockoutLists(infiles, outfiles):

    dbhandle = connect()

    outf = gzip.open(outfiles[0], "w")

    headers = ("strain",
               "gene_id",
               "gene_name",
               "ntranscripts",
               "contig",
               "strand",
               "stops-start",
               "stops-end")

    outf.write("%s\n" % "\t".join(headers))

    for infile in infiles:

        track = infile[:-len(".load")]
        strain = track[:-len("_alleles_genes")]

        statement = '''
                     SELECT DISTINCT '%(strain)s',
                            g.gene_id, 
                            i.gene_name, 
                                    g.ntranscripts,
                                    g.contig, g.strand,
                                    g.stop_codons_start,
                                    g.stop_codons_end
                         FROM %(track)s as g,
                            annotations.transcript_info AS i
                        WHERE g.gene_id = i.gene_id AND g.is_nmd_knockout
        ''' % locals()

        outf.write("\n".join(["\t".join(map(str, x))
                              for x in Database.executewait(dbhandle, statement).fetchall()]) + "\n")

    outf.close()

    headers = ("gene_id",
               "gene_name",
               "nmd_knockout_total",
               "strains")

    outf = gzip.open(outfiles[1], "w")
    outf.write("%s\n" % "\t".join(headers))

    columns = ["%s_nmd_knockout" % t for t in TRACKS]
    fields = ",".join(columns)

    statement = '''
        SELECT DISTINCT gene_id, gene_name, nmd_knockout_total, %(fields)s 
        FROM view_genes WHERE nmd_knockout_total > 0
        ''' % locals()

    data = list(dbhandle.execute(statement))

    d = dict(
        zip(["gene_id", "gene_name", "nmd_knockout_total"] + columns, zip(*data)))

    c = []
    for x in range(len(d["gene_id"])):
        s = []
        for t in TRACKS:
            if d["%s_nmd_knockout" % t][x] != 0:
                s.append(t)
        c.append(",".join(s))

    for t in TRACKS:
        del d["%s_nmd_knockout" % t]

    for d, strains in zip(data, c):
        outf.write("\t".join(map(str, d[:3])) + "\t%s\n" % strains)
    outf.close()

###################################################################
###################################################################
###################################################################


@merge(loadEffects, "polyphen.input")
def buildPolyphenInput(infiles, outfile):
    '''build polyphen input file.

    SNPS across all species are aggregated into a single
    file to avoid multiple submissions for the same variant.

    Mapping to Uniprot ids was not successful - 40% of the
    SNPs would have been lost. Hence I map to ensembl protein
    identifiers. Note that the sequence file is then to be 
    submitted to POLYPHEN as well.

    Note that this method outputs 1-based coordinates for polyphen,
    while the coordinates in the .map file are still 0-based.

    SNPs are assigned a snp_id and a locus_id. The snp_id refers
    to the SNP within a peptide sequence while the locus_id refers
    to the genomic location. If there are alternative
    transcripts overlapping a SNP, the same SNP will get two
    snp_ids, but the same locus_id. As the peptide background might
    be different for the same SNP depending on the transcript,
    its effect needs to be predicted twice.
    '''

    statement = '''SELECT
        transcript_id,
        cds_start,
        cds_end,
        orig_codons,
        variant_codons,
        orig_na,
        variant_na,
        contig,
        snp_position
    FROM %(table)s_cds
    WHERE variant_code = '=' AND code = 'N'
    '''

    dbhandle = connect()
    cc = dbhandle.cursor()

    infiles.sort()

    # uniprot mapping:
    # map_transcript2id = dict(
    # cc.execute( "SELECT transcript_id, trembl_acc FROM ensembl2uniprot WHERE
    # trembl_acc IS NOT NULL").fetchall() )

    # ensembl mapping
    map_transcript2id = dict(
        cc.execute("SELECT transcript_id, protein_id FROM annotations.transcript_info WHERE protein_id IS NOT NULL").fetchall())

    total_counts = E.Counter()
    notfound, found = set(), set()

    outf_map = open(outfile + ".map", "w")
    outf_map.write(
        "snp_id\ttrack\ttranscript_id\tprotein_id\tprotein_pos\tlocus_id\tcontig\tpos\tphase\n")

    outf = open(outfile, "w")

    snps = {}
    locus_ids = {}

    for infile in infiles:

        table = P.toTable(infile)
        track = table[:-len("_effects")]
        cc.execute(statement % locals())

        counts = E.Counter()

        snp_id = 0
        for transcript_id, cds_start, cds_end, orig_codons, variant_codons, orig_na, variant_na, contig, pos in cc:

            counts.input += 1

            if transcript_id not in map_transcript2id:
                notfound.add(transcript_id)
                counts.not_found += 1
                continue

            if "," in variant_codons:
                counts.heterozygous += 1
                continue

            for phase in range(0, 3):
                if orig_na[phase].lower() != variant_na[phase].lower():
                    break

            pid = map_transcript2id[transcript_id]
            # one-based coordinates
            peptide_pos = int(math.floor(cds_start / 3.0)) + 1
            key = "%s-%i-%s" % (pid, peptide_pos, variant_codons)

            if key in snps:
                snp_id = snps[key]
            else:
                snp_id = len(snps)
                snps[key] = snp_id
                outf.write("snp%010i\t%s\t%i\t%s\t%s\n" %
                           (snp_id,
                            pid,
                            peptide_pos,
                            orig_codons,
                            variant_codons,
                            ))
                counts.output += 1

            locus_key = "%s-%i-%s" % (contig, pos, variant_codons)
            if locus_key not in locus_ids:
                locus_ids[locus_key] = len(locus_ids)

            # use 0-based coordinates throughout, including peptide pos
            outf_map.write("snp%010i\t%s\t%s\t%s\t%i\tloc%010i\t%s\t%i\t%i\n" %
                           (snp_id,
                            track,
                            transcript_id,
                            pid,
                            peptide_pos - 1,
                            locus_ids[locus_key],
                            contig,
                            pos,
                            phase))

            found.add(transcript_id)

        total_counts += counts

        E.info("%s: %s" % (table, str(counts)))

    outf.close()
    outf_map.close()

    E.info("%s: transcripts: %s found, %i not found" % (table,
                                                        len(found),
                                                        len(notfound)))

    E.info("total=%s, snp_ids=%i, locus_ids=%i" %
           (str(total_counts), len(snps), len(locus_ids)))
    if notfound:
        E.warn("%i transcripts had SNPS that were ignored because there was no uniprot accession" %
               len(notfound))
        E.warn("notfound: %s" % ",".join(notfound))

    statement = '''sort -k2,2 -k3,3n %(outfile)s > %(outfile)s.tmp; mv %(outfile)s.tmp %(outfile)s'''

    P.run()

###################################################################
###################################################################
###################################################################


@transform(buildPolyphenInput, suffix(".input"), ".features")
def buildPolyphenFeatures(infile, outfile):
    '''run polyphen on the cluster.

    To do this, first send uniref to all nodes:

    python ~/cgat/cluster_distribute.py 
           --collection=andreas 
           /net/cpp-group/tools/polyphen-2.0.18/nrdb/uniref100*.{pin,psd,psi,phr,psq,pal}
    '''

    nsnps = len([x for x in open(infile)])

    to_cluster = True
    stepsize = max(int(nsnps / 200000.0), 1000)
    job_array = (0, nsnps, stepsize)
    E.info("running array jobs on %i snps" % nsnps)

    scratchdir = os.path.join(os.path.abspath("."), "scratch")
    try:
        os.mkdir(scratchdir)
    except OSError:
        pass

    resultsdir = outfile + ".dir"
    try:
        os.mkdir(resultsdir)
    except OSError:
        pass

    filename_peptides = os.path.join(PARAMS["annotations_dir"],
                                     PARAMS_ANNOTATIONS["interface_peptides_fasta"])

    statement = '''
    %(polyphen_home)s/bin/run_pph.pl
       -s %(filename_peptides)s
       -b %(polyphen_blastdb)s
       -d %(scratchdir)s
       %(infile)s > %(resultsdir)s/%(outfile)s.$SGE_TASK_ID 2> %(resultsdir)s/%(outfile)s.err.$SGE_TASK_ID
    '''
    P.run()

    to_cluster = False
    job_array = None

    statement = '''find %(resultsdir)s -name "*.err.*" -exec cat {} \; > %(outfile)s.log'''
    P.run()

    statement = '''find %(resultsdir)s -not -name "*.err.*" -exec cat {} \; > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################
# do not run in parallel. run_weka.pl creates a $testfile
# that is not unique. run_weka.pl and pph2arff.pl could either
# be patched or the following jobs run in sequence.


@jobs_limit(1)
@files([(buildPolyphenFeatures, "polyphen_%s.output.gz" % x, x)
        for x in P.asList(PARAMS["polyphen_models"])])
def runPolyphen(infile, outfile, model):
    '''run POLYPHEN on feature tables to classify SNPs.
    '''

    to_cluster = True

    # options
    # -f: feature set, default is F11
    # -c: classifier, default is NBd (Naive Bayes with discretization)
    # -l: model name, default is HumDiv

    statement = '''
    %(polyphen_home)s/bin/run_weka.pl 
           -l %(polyphen_home)s/models/%(model)s.UniRef100.NBd.f11.model
           %(infile)s 
    | gzip 
    > %(outfile)s 
    2> %(outfile)s.log
    '''

    P.run()

###################################################################
###################################################################
###################################################################


@transform(buildPolyphenInput, suffix(".input"), "_map.load")
def loadPolyphenMap(infile, outfile):
    '''load polyphen input data.'''

    table = P.toTable(outfile)
    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --add-index=snp_id 
              --add-index=track,transcript_id
              --add-index=contig,pos
              --add-index=protein_id
              --add-index=transcript_id
              --table=%(table)s 
    < %(infile)s.map
    > %(outfile)s
    '''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(runPolyphen, suffix(".output.gz"), ".load")
def loadPolyphen(infile, outfile):
    '''load polyphen results.'''

    table = P.toTable(outfile)

    statement = '''
    gunzip 
    < %(infile)s
    | perl -p -e "s/o_acc/protein_id/; s/ +//g; s/^#//;"
    |python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --add-index=snp_id 
              --add-index=protein_id
              --table=%(table)s 
              --map=effect:str
    > %(outfile)s
    '''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(loadPolyphen, suffix(".load"), ".genestats")
def analysePolyphen(infile, outfile):
    '''compute enrichment of SNPs within genes
    and deleterious SNPs within SNPs within genes.

    del: enrichment of deleterious snps within snps per gene
    len: enrichment of snps within genes
    com: enrichment of deleterious snps within gene
    '''

    table = P.toTable(infile)
    tablename_map = "polyphen_map"

    dbhandle = connect()
    cc = dbhandle.cursor()

    statement = '''
        SELECT i.gene_id,
               COUNT(DISTINCT map.locus_id) as nsnps, 
               COUNT(DISTINCT case t.prediction when 'possiblydamaging' then map.locus_id when 'probablydamaging' then map.locus_id else NULL end) AS ndeleterious,
               MAX(s.length)
               FROM %(table)s as t, 
                    %(tablename_map)s as map, 
                    annotations.protein_stats as s,
                    annotations.transcript_info as i 
        WHERE map.snp_id = t.snp_id AND 
              i.transcript_id = map.transcript_id AND
              s.protein_id = map.protein_id
        GROUP BY i.gene_id
     ''' % locals()

    data = cc.execute(statement).fetchall()

    statement = '''SELECT DISTINCT i.gene_id, MAX(s.length) 
                   FROM annotations.transcript_info AS i, annotations.protein_stats AS s 
                   WHERE s.protein_id = i.protein_id 
                   GROUP BY i.gene_id'''
    gene_ids = cc.execute(statement).fetchall()

    total_nsnps = sum([x[1] for x in data])
    total_ndel = sum([x[2] for x in data])
    total_length = sum([x[1] for x in gene_ids])
    del_p = float(total_ndel) / total_nsnps
    len_p = float(total_nsnps) / total_length
    com_p = float(total_ndel) / total_length

    E.info("del: background probability: %i/%i = %f" %
           (total_ndel, total_nsnps, del_p))
    E.info("len: background probability: %i/%i = %f" %
           (total_nsnps, total_length, len_p))
    E.info("com: background probability: %i/%i = %f" %
           (total_ndel, total_length, com_p))

    outf = open(outfile, "w")
    outf.write("\t".join(("gene_id", "code",
                          "length", "nsnps", "ndel",
                          "del_p", "del_pvalue", "del_qvalue",
                          "len_p", "len_pvalue", "len_qvalue",
                          "com_p", "com_pvalue", "com_qvalue", )) + "\n")

    del_pvalues, len_pvalues, com_pvalues = [], [], []
    for gene_id, nsnps, ndel, length in data:

        # use -1, because I need P( x >= X)
        # sf = 1 - cdf and cdf = P( x <= X ), thus sf = 1 - P( x <= X ) = P (x
        # > X ).
        del_pvalues.append(scipy.stats.binom.sf(ndel - 1, nsnps, del_p))
        len_pvalues.append(
            scipy.stats.binom.sf(nsnps - 1, int(round(length)), len_p))
        com_pvalues.append(
            scipy.stats.binom.sf(ndel - 1, int(round(length)), com_p))

    if len(del_pvalues) > 10:
        del_qvalues = Stats.doFDR(del_pvalues).mQValues
    else:
        E.warn("no FDR computed for del")
        del_qvalues = del_pvalues

    if len(len_pvalues) > 10:
        len_qvalues = Stats.doFDR(len_pvalues).mQValues
    else:
        E.warn("no FDR computed for del")
        len_qvalues = len_pvalues

    if len(com_pvalues) > 10:
        com_q = Stats.doFDR(com_pvalues).mQValues
    else:
        E.warn("no FDR computed for com")
        com_qvalues = com_pvalues

    fdr = PARAMS["polyphen_fdr"]

    found = set()

    for a, del_pvalue, del_qvalue, len_pvalue, len_qvalue, com_pvalue, com_qvalue in \
            zip(data,
                del_pvalues, del_qvalues,
                len_pvalues, len_qvalues,
                com_pvalues, com_qvalues,
                ):
        gene_id, nsnps, ndel, length = a
        found.add(gene_id)

        del_p = float(ndel) / nsnps
        len_p = float(nsnps) / length

        code = "".join([str(int(x < fdr))
                       for x in (del_qvalue, len_qvalue, com_qvalue)])

        outf.write("\t".join((gene_id,
                              code,
                              "%i" % int(round(length)),
                              "%i" % int(nsnps),
                              "%i" % int(ndel),
                              "%6.4f" % del_p,
                              "%6.4g" % del_pvalue,
                              "%6.4g" % del_qvalue,
                              "%6.4f" % len_p,
                              "%6.4g" % len_pvalue,
                              "%6.4g" % len_qvalue,
                              "%6.4f" % com_p,
                              "%6.4g" % com_pvalue,
                              "%6.4g" % com_qvalue,
                              )) + "\n")

    # add missing genes:
    code = "---"
    for gene_id, length in gene_ids:
        if gene_id in found:
            continue
        outf.write("\t".join((gene_id,
                              code,
                              "%i" % int(round(length)),
                              "%i" % 0,
                              "%i" % 0,
                              "%6.4f" % 0,
                              "%6.4g" % 1,
                              "%6.4g" % 1,
                              "%6.4f" % 0,
                              "%6.4g" % 1,
                              "%6.4g" % 1,
                              "%6.4f" % 0,
                              "%6.4g" % 1,
                              "%6.4g" % 1,
                              )) + "\n")

    outf.close()

###################################################################
###################################################################
###################################################################


@transform(analysePolyphen, suffix(".genestats"), "_genestats.load")
def loadPolyphenAnalysis(infile, outfile):
    '''load polyphen analysis results.'''

    table = P.toTable(outfile)

    statement = '''
    cat < %(infile)s
    |python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --add-index=gene_id 
              --map=code:str
              --table=%(table)s 
    > %(outfile)s
    '''
    P.run()

###################################################################
###################################################################
###################################################################
# @files( ( ( buildPeptideFasta, "panther.scores"), ))
# def preparePanther( infile, outfile ):
#     '''lookup peptide sequences with panther.

#     The actual snps will get scored in the next step.
#     This step takes a while, but could be sped up
#     easily by parallelization.
#     '''

#     to_cluster = True

#     if type(infile) in (types.ListType, types.TupleType):
#         infile = infile[0]

#     tmpdir = P.getTempDir( "." )

#     statement = '''
#     (PERL5LIB=%(panther_home)s/lib:$PERL5LIB;
#      perl %(panther_home)s/pantherScore.pl
#                       -l %(panther_library)s
#                       -D B -V -n
#                       -i %(infile)s
#                       -o %(outfile)s
#                       -T %(tmpdir)s )
#     '''

#     P.run()

#     shutil.rmtree( tmpdir )

# ###################################################################
# ###################################################################
# ###################################################################
# @files( ( ( (buildPolyphenInput, preparePanther), "panther.output" ), ) )
# def runPanther( infiles, outfile):
#     '''run PANTHER analysis.
#     '''

# to_cluster = True

#     filename_snps, filename_scores = infiles

#     tmpdir = P.getTempDir( "." )

#     peptides = PARAMS["peptides"]
#     tmpfilename_snps = P.getTempFilename(".")

#     statement = '''
#     awk '{printf("%%s|%%s|%%s|%%s;%%s\\n",
#                  $1,$2,$3,$4,$5);}'
#     < %(filename_snps)s > %(tmpfilename_snps)s
#     '''
# P.run()

#     statement = '''
#     (PERL5LIB=%(panther_home)s/lib:$PERL5LIB;
#      PATH=%(panther_home)s:$PATH;
#      awk '{printf("%%s|%%s|%%s|%%s;%%s\\n",
#                  $1,$2,$3,$4,$5);}'
#     < %(filename_snps)s
#     | %(cmd-farm)s --split-at-lines=2000 --log=%(outfile)s.log -v 10 --output-header --env=PERL5LIB --env=PATH
#     "perl %(panther_home)s/snp_analysis.pl
#                       -l %(panther_library)s
#                       -c %(filename_scores)s
#                       -s %%STDIN%%
#                       -f %(peptides)s
#                       -b %(panther_home)s/BLOSUM62
#                       -V
#                       -p %(panther_home)s/uprior.9comp
#                       -o %%STDOUT%%
#                       -T %(tmpdir)s"
#     > %(outfile)s 2> %(outfile)s.log )
#     '''

#     P.run()
#     shutil.rmtree( tmpdir )
#     os.unlink( tmpfilename_snps )


# ###################################################################
# ###################################################################
# ###################################################################
# @transform( runPanther, suffix(".output"), ".load")
# def loadPanther( infile, outfile ):
#     '''load panther results.'''

#     table = P.toTable( outfile )

#     statement = '''
#     perl -p -e "s/snpId/snp_id/; s/seqId/protein_id/; s/HMM /hmm/g;"
#     < %(infile)s
#     |python %(scriptsdir)s/csv2db.py %(csv2db_options)s
#               --add-index=snp_id
#               --add-index=protein_id
#               --table=%(table)s
#     > %(outfile)s
#     '''
#     P.run()

###################################################################
###################################################################
###################################################################
@split(loadPolyphenMap, ("counts_shared.matrix", "counts_segregation.matrix", "counts_pid.matrix", "counts_distance.matrix", "counts.tree"
                         ))
def buildSharedSNPMatrix(infiles, outfiles):
    '''build matrix of shared coding nonsynonymous SNPs.

    Counts are per locus id.

    Percent identities are only within coding segregating loci
    and thus do not reflect the real divergence.

    '''

    dbhandle = connect()
    cc = dbhandle.cursor()

    segregating_sites = cc.execute(
        'SELECT COUNT( DISTINCT locus_id) FROM polyphen_map').fetchone()[0]

    statement = '''SELECT DISTINCT locus_id, track FROM polyphen_map ORDER BY locus_id'''
    cc.execute(statement)

    matrix = collections.defaultdict(int)
    for k, vals in itertools.groupby(cc, key=lambda x: x[0]):
        tracks = [x[1] for x in list(vals)]
        for t1 in tracks:
            matrix[(t1, t1)] += 1
        if len(tracks) > 1:
            for t1, t2 in itertools.combinations(tracks, 2):
                matrix[(t1, t2)] += 1
                matrix[(t2, t1)] += 1

    all_tracks = set([x[0] for x in matrix.keys()] + [x[1]
                     for x in matrix.keys()])

    # output matrix with shared SNPs.
    outf = open(outfiles[0], "w")
    outf.write("track\t%s\n" % "\t".join(all_tracks))
    for track1 in all_tracks:
        outf.write("%s" % track1)
        for track2 in all_tracks:
            outf.write("\t%i" % matrix[(track1, track2)])
        outf.write("\n")
    outf.close()

    # output matrix with shared segregating sites as
    # distance matrix
    outf = open(outfiles[1], "w")
    outf.write("track\t%s\n" % "\t".join(all_tracks))
    for track1 in all_tracks:
        outf.write("%s" % track1)
        for track2 in all_tracks:
            if track1 == track2:
                outf.write("\t%i" % 0)
            else:
                outf.write("\t%i" %
                           (segregating_sites - matrix[(track1, track2)]))
        outf.write("\n")
    outf.close()

    # output matrix as percent identity matrix
    # percent identity is given as
    # segregating sites - sites where strains differ = segregating_sites - (matrix[i,i] + matrix[j,j] - 2 * matrix[i,j])
    # simplifies to:
    # segsites - matrix[i,i] -matrix[j,j] +
    # divided by the total number of segregating sites
    outf = open(outfiles[2], "w")
    outf.write("track\t%s\n" % "\t".join(all_tracks))
    pids = {}
    for track1 in all_tracks:
        outf.write("%s" % track1)
        for track2 in all_tracks:
            a = segregating_sites - \
                (matrix[(track1, track1)] + matrix[(track2, track2)] -
                 2 * matrix[(track1, track2)])
            pid = 100.0 * a / segregating_sites
            outf.write("\t%6.4f" % pid)
            pids[(track1, track2)] = pid
        outf.write("\n")
    outf.close()

    # distance matrix
    outf = open(outfiles[3], "w")
    outf.write("track\t%s\n" % "\t".join(all_tracks))
    for track1 in all_tracks:
        outf.write("%s" % track1)
        for track2 in all_tracks:
            val = 100.0 - pids[(track1, track2)]
            outf.write("\t%6.4f" % val)
        outf.write("\n")
    outf.close()

    outfile_distance, outfile_tree = outfiles[3], outfiles[4]

    # build tree
    statement = '''python %(scriptsdir)s/matrix2matrix.py
       --output-format=phylip
    < %(outfile_distance)s
    | python %(scriptsdir)s/matrix2tree.py
       --method=nj
    > %(outfile_tree)s
    '''
    P.run()

###################################################################
###################################################################
###################################################################
# Enrichment analysis
###################################################################


@files(((None, "workspace_genomic.bed", "genomic"),
        (None, "workspace_cds.bed", "cds"),
        ))
def buildQTLWorkspaces(infile, outfile, workspace):
    PEnrichment.buildWorkSpace(outfile, workspace)


@files((("%s.fasta" % PARAMS["genome"], "workspace_isochores.bed.gz"), ))
def buildEnrichmentIsochores(infile, outfile):
    PEnrichment.buildIsochoresGC(infile, outfile)


@follows(mkdir("enrichment.dir"), loadPolyphen, loadPolyphenMap)
@transform("*_effects.load",
           regex("(.*)_effects.load"),
           r"enrichment.dir/\1.deleterious.bed.gz")
def buildDeleteriousSNPs(infile, outfile):

    track = infile[:-len("_effects.load")]

    outf = gzip.open(outfile, "w")
    outf.write("track name=%s.deleterious\n" % track)

    dbhandle = connect()
    cc = dbhandle.cursor()

    statement = '''SELECT DISTINCT map.contig, map.pos 
                          FROM polyphen_map AS map,
                          polyphen_HumDiv as result
                          WHERE map.track = '%(track)s'
                                AND map.snp_id = result.snp_id
                                AND (result.prediction = 'possiblydamaging'
                                    OR result.prediction = 'probablydamaging')
                          ''' % locals()

    cc.execute(statement)

    for contig, pos in cc:
        outf.write("%s\t%i\t%i\n" % (contig, pos, pos + 1))

    outf.close()


@follows(mkdir("enrichment.dir"), loadPolyphen, loadPolyphenMap)
@transform("*_effects.load",
           regex("(.*)_effects.load"),
           r"enrichment.dir/\1.benign.bed.gz")
def buildBenignSNPs(infile, outfile):

    track = infile[:-len("_effects.load")]

    outf = gzip.open(outfile, "w")
    outf.write("track name=%s.benign\n" % track)

    dbhandle = connect()
    cc = dbhandle.cursor()

    statement = '''SELECT DISTINCT map.contig, map.pos 
                          FROM polyphen_map AS map,
                          polyphen_HumDiv as result
                          WHERE map.track = '%(track)s'
                                AND map.snp_id = result.snp_id
                                AND NOT (result.prediction = 'possiblydamaging'
                                    OR result.prediction = 'probablydamaging')
                          ''' % locals()

    cc.execute(statement)

    for contig, pos in cc:
        outf.write("%s\t%i\t%i\n" % (contig, pos, pos + 1))

    outf.close()


@merge((buildBenignSNPs, buildDeleteriousSNPs),
       ("enrichment.dir/all.benign.bed.gz",
        "enrichment.dir/all.deleterious.bed.gz",
        "enrichment.dir/all.ambiguous.bed.gz", ),
       )
def mergeSNPs(infiles, outfiles):

    tmp1 = P.getTempFilename()
    tmp2 = P.getTempFilename()
    statement = '''zcat enrichment.dir/mouse*.benign.bed.gz 
                | grep -v "track" 
                | sort -k 1,1 -k2,2n 
                | uniq 
                > %(tmp1)s
    '''
    P.run()

    statement = '''zcat enrichment.dir/mouse*.deleterious.bed.gz 
                | grep -v "track" 
                | sort -k 1,1 -k2,2n 
                | uniq 
                > %(tmp2)s 
    '''
    P.run()

    statement = '''intersectBed -a %(tmp1)s -b %(tmp2)s 
                   | awk 'BEGIN {printf("track name=all.ambiguous\\n");} {print}' 
                   > enrichment.dir/all.ambiguous.bed'''
    P.run()

    statement = '''intersectBed -v -a %(tmp1)s -b enrichment.dir/all.ambiguous.bed
                   | awk 'BEGIN {printf("track name=all.benign\\n");} {print}' 
                   | gzip
                   > enrichment.dir/all.benign.bed.gz'''

    P.run()

    statement = '''intersectBed -v -a %(tmp2)s -b enrichment.dir/all.ambiguous.bed
                   | awk 'BEGIN {printf("track name=all.deleterious\\n");} {print}' 
                   | gzip
                   > enrichment.dir/all.deleterious.bed.gz'''

    P.run()

    statement = '''gzip enrichment.dir/all.ambiguous.bed'''
    P.run()

    os.unlink(tmp1)
    os.unlink(tmp2)


@merge((buildBenignSNPs, buildDeleteriousSNPs),
       "enrichment.dir/isochores.bed")
def buildSNPDensityIsochores(infile, outfile):
    '''build isochores with SNP density.'''

    statement = '''
           python %(scriptsdir)s/windows2gff.py 
                --genome=%(genome)s
                --fixed-width-windows=1000000
                --output-format=bed
           > tmp.bed'''
    P.run()

    statement = '''
           zcat enrichment.dir/mouse*.benign.bed.gz enrichment.dir/mouse*.deleterious.bed.gz
                | grep -v "track" 
                | sort -k 1,1 -k2,2n 
                | uniq > tmp2.bed
    '''
    P.run()


@merge(["/net/cpp-compute/backup/andreas/projects/mousestrains/data/qtl/jonathans/qtl_merged.bed",
        "/net/cpp-compute/backup/andreas/projects/mousestrains/data/qtl/jonathans/qtl_full.bed",
        "/net/cpp-compute/backup/andreas/projects/mousestrains/data/qtl/jonathans/qtl_rest.bed",
        ], "qtl.summary.tsv")
def QTLSummary(infiles, outfile):

    for infile in infiles:
        basename = os.path.basename(infile)
        statement = '''
        python %(scriptsdir)s/bed2gff.py
        < %(infile)s
        | python %(scriptsdir)s/gff2histogram.py
               --method=all
               --output-filename-pattern=%(outfile)s.%(basename)s
               --log=%(outfile)s.log
        > %(outfile)s
        '''
        P.run()


@follows(buildQTLWorkspaces)
@merge((buildDeleteriousSNPs, buildBenignSNPs, mergeSNPs), "qtl.table")
def runGATOnQTLs(infiles, outfile):
    '''run enrichment analysisusing the qtl definitions from
    Jonathan Flint's group.
    '''

    segments = IOTools.flatten(infiles)

    workspaces = ["workspace_cds.bed", ]

    annotations = ["/net/cpp-compute/backup/andreas/projects/mousestrains/data/qtl/jonathans/qtl_merged.bed",
                   "/net/cpp-compute/backup/andreas/projects/mousestrains/data/qtl/jonathans/qtl_full.bed",
                   "/net/cpp-compute/backup/andreas/projects/mousestrains/data/qtl/jonathans/qtl_rest.bed",
                   ]

    workspaces = " ".join(["--workspace-file=%s" % x for x in workspaces])
    annotations = " ".join(["--annotation-file=%s" % x for x in annotations])
    segments = " ".join(["--segment-file=%s" % x for x in segments])

    to_cluster = True
    job_options = "-l mem_free=8000M"

    statement = '''gatrun.py
                  %(workspaces)s
                  %(segments)s
                  %(annotations)s
                  --output-stats=annotations
                  --output-stats=workspaces
                  --output-filename-pattern=enrichment.dir/%%s.tsv
                  --force-output
                  --num-samples=10000
    > %(outfile)s
    '''
    P.run()


@follows(buildQTLWorkspaces)
@merge(mergeSNPs, "qtl_small.table")
def runGATOnQTLsSmall(infiles, outfile):
    '''run enrichment analysisusing the qtl definitions from
    Jonathan Flint's group.
    '''

    segments = IOTools.flatten(infiles)

    workspaces = ["workspace_cds.bed", ]

    annotations = ["/net/cpp-compute/backup/andreas/projects/mousestrains/data/qtl/jonathans/qtl_merged.bed",
                   "/net/cpp-compute/backup/andreas/projects/mousestrains/data/qtl/jonathans/qtl_full.bed",
                   "/net/cpp-compute/backup/andreas/projects/mousestrains/data/qtl/jonathans/qtl_rest.bed",
                   ]

    workspaces = " ".join(["--workspace-file=%s" % x for x in workspaces])
    annotations = " ".join(["--annotation-file=%s" % x for x in annotations])
    segments = " ".join(["--segment-file=%s" % x for x in segments])

    to_cluster = True
    job_options = "-l mem_free=8000M"

    statement = '''gatrun.py
                  %(workspaces)s
                  %(segments)s
                  %(annotations)s
                  --output-stats=annotations
                  --output-stats=workspaces
                  --output-filename-pattern=enrichment.dir/%%s.tsv
                  --force-output
                  --num-samples=10000
    > %(outfile)s
    '''
    P.run()


@transform((runGATOnQTLs, ), suffix(".table"), ".load")
def loadGATOnQTLs(infile, outfile):

    table = P.toTable(outfile)

    statement = '''
    cat < %(infile)s
    |python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --add-index=track
              --add-index=annotation
              --table=%(table)s
    > %(outfile)s
    '''
    P.run()

    stat_files = glob.glob("enrichment.dir/stats_*.tsv")

    for stat_file in stat_files:
        basename = os.path.basename(stat_file)
        table = os.path.splitext(basename)[0]
        statement = '''
        cat < %(stat_file)s
        |python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --add-index=track
              --add-index=contig
              --table=%(table)s
    >> %(outfile)s
    '''
        P.run()

###################################################################
###################################################################
###################################################################
# gene list analyses
###################################################################
############################################################


@files(((importMGI, "assignments.mgi"),))
def createMGI(infile, outfile):
    '''get GO assignments from MGI'''

    dbhandle = connect()

    statement = '''
                 SELECT DISTINCT 'MPheno.ontology', m2g.gene_id, a2p.phenotype_id, p.term, 'NA' 
                 FROM mgi_marker2gene as m2g, 
                      mgi_marker2allele as m2a, 
                      mgi_allele2phenotype as a2p, 
                      mgi_phenotypes as p 
                 WHERE m2g.marker_id = m2a.marker_id AND 
                       a2p.allele_id = m2a.allele_id AND 
                       p.phenotype_id = a2p.phenotype_id
                '''
    cc = dbhandle.cursor()
    data = cc.execute(statement).fetchall()

    outf = open(outfile, "w")
    outf.write("\n".join(["\t".join(x) for x in data]) + "\n")
    outf.close()


####################################################################
def buildGeneMatrix(tracks, analysis, statement, outfile):
    '''build a gene matrix.

    A gene matrix is an n x m matrix for n genes and m gene lists.
    Each column contains a 1 if a gene is present in a gene list,
    otherwise it is 0.
    '''

    dbhandle = connect()
    cc = dbhandle.cursor()

    all_genes = [x[0] for x in cc.execute(
        '''SELECT DISTINCT gene_id FROM annotations.gene_info''' % locals() )]

    gene2row = dict([(x[1], x[0]) for x in enumerate(all_genes)])
    matrix = numpy.zeros(
        (len(all_genes), len(analysis) * len(tracks)), numpy.int)

    col = 0
    for track in tracks:

        for label, field_where in analysis:
            genes = [x[0] for x in cc.execute(statement % locals())]
            for gene_id in genes:
                matrix[gene2row[gene_id]][col] = 1
            col += 1

    outf = open(outfile, "w")
    outf.write("gene_id\t%s\n" % "\t".join("%s_%s" %
               (x, y[0]) for x, y in itertools.product(tracks, analysis)))
    for gene_id in all_genes:
        outf.write("%s\t%s\n" %
                   (gene_id, "\t".join(map(str, matrix[gene2row[gene_id]]))))
    outf.close()

####################################################################


@merge(summarizeAllelesPerGene, "effects.genematrix")
def buildGeneMatrixEffects(infiles, outfile):
    '''build gene matrix with consequences data.

    Note that this analysis is confounded by gene length.
    '''

    analysis = (("benign", "benign"),
                ("probablydamaging", "probablydamaging"),
                ("possiblydamaging", "possiblydamaging"),
                ("unknown", "unknown"))

    tracks = [x[:-len("_alleles_genes.load")] for x in infiles]

    statement = '''
            SELECT DISTINCT i.gene_id 
            FROM
            transcript_info AS i,
            polyphen_map AS map,
            polyphen_HumVar as polyphen
            WHERE 
            polyphen.snp_id = map.snp_id AND
            map.track = '%(track)s' AND
            map.transcript_id = i.transcript_id AND
            prediction = '%(field_where)s' '''

    buildGeneMatrix(tracks, analysis, statement, outfile)

####################################################################


@merge(summarizeAllelesPerGene, "alleles.genematrix")
def buildGeneMatrixAlleles(infiles, outfile):
    '''build gene matrix from alleles results

    ``options`` is a tuple of (``track``, ``analysis``, ``ontology``)

    ``analysis`` can be:

    stoptruncated
        genes that are truncated due to stops
    nmdknockout
        genes that are knocked out due to NMD
    splicetruncated
        genes that are truncated due to deleted splice sites
    knockout
        any of the above

    Note that the analysis here needs to be background adjusted.
    NMD transcripts and splice truncated transcripts are only 
    multiple exon transcripts, while stoptruncated ones are
    only single exon ones.

    '''

    analysis = (("stoptruncated", "e.is_truncated"),
                ("nmdknockout", "e.is_nmd_knockout"),
                ("splicetruncated", "e.is_splice_truncated"),
                ("knockout", "(e.is_nmd_knockout or e.is_truncated or e.is_splice_truncated)"))

    tracks = [x[:-len("_alleles_genes.load")] for x in infiles]

    statement = '''
            SELECT DISTINCT e.gene_id 
            FROM
                  %(track)s_alleles_genes AS e
            WHERE 
                %(field_where)s
            '''

    buildGeneMatrix(tracks, analysis, statement, outfile)

####################################################################


@merge(summarizeEffectsPerGene, "consequences.genematrix")
def buildGeneMatrixConsequences(infiles, outfile):
    '''build gene matrix from effecs results

    nmdknockouttranscript
        genes for which one transcript has been knocked out
        due to NMD
    nmdaffectedtranscript
        genes in which one transcript is affected by NMD
    nmdknockoutgenes
        genes in which all transcripts have been knocked out 
        due to NMD

    Note that the analysis here needs to be background adjusted.
    For example, NMD transcripts are only multiple exon transcripts.
    '''

    analysis = (("nmdknockouttranscript", "e.nmd_knockout > 0"),
                ("nmdaffectedtranscript", "e.nmd_affected > 0"),
                ("nmdknockoutgenes", "e.nmd_knockout = e.ntranscripts"))

    tracks = [x[:-len("_effects_genes.load")] for x in infiles]

    statement = '''
            SELECT DISTINCT e.gene_id 
            FROM
                  %(track)s_effects_genes AS e
            WHERE 
                %(field_where)s
            '''

    buildGeneMatrix(tracks, analysis, statement, outfile)

####################################################################


@follows(buildGeneMatrixConsequences,
         buildGeneMatrixAlleles,
         buildGeneMatrixEffects)
@files([((x, y), "%s_vs_%s.gla" % (re.sub(".genematrix", "", x),
                                   re.sub("assignments.", "", y)))
        for x, y in
        itertools.product(
            glob.glob("*.genematrix"),
            glob.glob("assignments.*"))
        if not y.endswith(".log")])
def runGeneListAnalysis(infiles, outfile):
    '''run a gene list analysis.'''

    genematrix, assignments = infiles

    to_cluster = True

    try:
        options = "--fdr-lambda=%(genelist_analysis_qvalue_lambda)f" % PARAMS
    except TypeError:
        options = ""

    statement = '''
    python %(scriptsdir)s/genelist_analysis.py 
           --format=matrix 
           --assignments-tsv-file=%(assignments)s 
           --fdr 
           --fdr-method=%(genelist_analysis_qvalue_method)s
           --log=%(outfile)s.log
           %(options)s
    < %(genematrix)s
    > %(outfile)s
    '''

    P.run()

###########################################################################


@transform(runGeneListAnalysis, suffix(".gla"), "_gla.load")
def loadGeneListAnalysis(infile, outfile):
    '''load gene list analysis results.'''
    table = P.toTable(outfile)

    statement = '''
    cat < %(infile)s
    |python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --table=%(table)s
              --add-index=gene_list
              --add-index=pvalue
              --add-index=fdr
    > %(outfile)s
    '''
    P.run()

# ############################################################################
# @merge( runGOAnalysesOnAlleles, "alleles_go.load")
# def loadGOs( infile, outfile ):
#     '''load go results.'''
#     tablename = P.toTable( outfile )
#     PGO.loadGOs( infile, outfile, tablename )

# ############################################################################
# @transform( runGOAnalysesOnAlleles, suffix(".go"), "_go.load")
# def loadGO( infile, outfile ):
#     '''load go results.'''
#     tablename = P.toTable( outfile )
#     PGO.loadGO( infile, outfile, tablename )

# ############################################################################
# @transform( runGOAnalysesOnAlleles, suffix(".goslim"), "_goslim.load")
# def loadGOSlim( infile, outfile ):
#     '''load goslim results.'''
#     tablename = P.toTable( outfile )
#     PGO.loadGO( infile, outfile, tablename )

# ############################################################################
# @files( ((loadGOs, "goresults.table"),))
# def mergeGO( infile, outfile ):
#     '''merge all GO anlyses.

#     * collect all P-Values for all categories and experiments.
#     * compute stats on it
#     '''

#     dbhandle = sqlite3.connect( PARAMS["database"] )

#     statement = '''SELECT track, geneset, annotationset, category, min(pover,punder)
#                    FROM alleles_go'''
#     cc = dbhandle.cursor()
#     data = cc.execute(statement).fetchall()

#     pvalues = [ x[4] for x in data ]
#     E.info( "analysing %i pvalues" % len(pvalues ))

#     fdr = Stats.doFDR( pvalues )
#     E.info( "got %i qvalues" % len(fdr.mQValues ))

#     for d, qvalue in zip( data, fdr.mQValues ):
#         if qvalue > 0.05: continue
#         print data, qvalue

#     Database.executewait( dbhandle, '''ALTER TABLE %(table)s ADD COLUMN is_coding FLOAT''' % locals())

###################################################################
###################################################################
###################################################################
###################################################################
###################################################################
###################################################################


@follows(loadHumanOrthologs,
         loadGene2Omim)
def prepare():
    pass


@follows(buildEffects, loadEffects, mergeEffects, summarizeEffectsPerGene)
def consequences():
    pass


@follows(buildAlleles, loadAlleles,
         summarizeAllelesPerTranscript,
         summarizeAllelesPerGene,
         combineSummaryAllelesPerGene)
def alleles():
    pass


@follows(loadPolyphen, loadPolyphenMap, loadPolyphenAnalysis)
def effects():
    pass


@follows(loadAnnotations, loadAnnotationsSummary)
def annotations():
    pass


@follows(consequences, effects, alleles, annotations)
def full():
    pass


@follows(buildQTLWorkspaces,
         runGATOnQTLs,
         runGATOnQTLsSmall)
def qtl():
    pass


@follows(buildGeneMatrixConsequences,
         buildGeneMatrixAlleles,
         buildGeneMatrixEffects,
         runGeneListAnalysis,
         loadGeneListAnalysis,
         )
def go():
    pass


@files([(None, "clone.log"), ])
def clone(infile, outfile):
    '''clone a pipeline using symbolic links.'''

    src_dir, level = sys.argv[-2:]

    if not os.path.exists(src_dir):
        raise IOError("directory '%s' does not exist" % src_dir)

    if not os.path.exists(os.path.join(src_dir, "pipeline.ini")):
        raise IOError("directory '%s' is not a pipeline" % src_dir)

    if level in ("data", ):
        P.execute("ln -fs %(src_dir)s/*.pileup.* . ")
        P.execute("ln -fs %(src_dir)s/genome.* . ")


###################################################################
@merge((alleles, prepare),  "genes.views")
def createViewGenes(infile, outfile):
    '''create view in database for genes.

    This view aggregates all information on a per-gene
    basis. There is only a single entry per gene.
    '''

    # can not create views across multiple database, so use table
    view_type = "TABLE"

    dbhandle = connect()
    Database.executewait(
        dbhandle, "DROP %(view_type)s IF EXISTS view_genes" % locals())

    knockouts = ",".join(
        ["nmd.%s AS %s_nmd_knockout" % (track, track) for track in TRACKS])

    statement = '''
    CREATE %(view_type)s view_genes AS
    SELECT i.gene_id AS gene_id, 
           i.gene_name AS gene_name,
           nmd.total AS nmd_knockout_total, 
           %(knockouts)s,
           human_ortho.hs_gene_id AS hs_gene_id,
           human_ortho.ds AS hs_ds,
           omim.mim_gene_id as omim_gene_id,
           omim.mim_morbid_description as omim_description,
           omim.mim_morbid_id as omim_morbid_id
    FROM annotations.gene_info AS i,
         summary_alleles_genes_is_knockout AS nmd ON i.gene_id = nmd.gene_id
    LEFT JOIN orthologs AS human_ortho ON 
          human_ortho.gene_id = i.gene_id AND 
          human_ortho.orthology_type = "ortholog_one2one"
    LEFT JOIN gene2omim as omim ON omim.gene_id = human_ortho.hs_gene_id
    '''

    Database.executewait(dbhandle, statement % locals())


@follows(createViewGenes)
def views():
    pass


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''

    E.info("starting documentation build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.'''

    E.info("updating documentation")
    P.run_report(clean=False)

if __name__ == "__main__":

    # P.checkFiles( ("genome.fasta", "genome.idx" ) )
    sys.exit(P.main(sys.argv))
