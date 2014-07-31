"""
===========================
Variant annotation pipeline
===========================

:Author: Andreas Heger & David Sims
:Release: $Id$
:Date: |today|
:Tags: Python

The Variants pipeline attempts to annotate variants in
a :term:`vcf` formatted file. It computes
 
   1. The effects of SNPs on transcripts and genes
   2. The effects of indels on transcripts and genes

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
import gzip
import os
import itertools
import re
import math
import collections
import sqlite3
import CGAT.Experiment as E
import CGAT.Database as Database
import scipy.stats
import CGAT.Stats as Stats
import pysam

# only update R if called as pipeline
# otherwise - failure with sphinx
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri

###################################################################
###################################################################
###################################################################
# Pipeline configuration
import CGAT.Pipeline as P
P.getParameters(["%s/pipeline.ini" %
                os.path.splitext(__file__)[0], "../pipeline.ini",
                 "pipeline.ini"])

PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    on_error_raise=__name__ == "__main__")

SEPARATOR = "|"

###################################################################
###################################################################
###################################################################
# Helper functions mapping tracks to conditions, etc
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
# Database connectivity


def connect():
    '''connect to the database.
    This method also attaches to the annotation database.'''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh

###################################################################
###################################################################
# Annotations
###################################################################
###################################################################


@files([("variants.vcf.gz", "%s.annotations.gz" % x, x) for x in TRACKS])
def buildAnnotations(infile, outfile, sample):
    """annotate snps with gene set."""

    to_cluster = True
    bases = "annotations_bases"

    statement = """python %(scriptsdir)s/snp2table.py 
                       --input-format=vcf
                       --filename-vcf=%(infile)s
                       --vcf-sample=%(sample)s
                       --genome-file=%(genome_dir)s/%(genome)s 
                       --filename-annotations=%(bases)s 
                       --log=%(outfile)s.log 
                   | gzip > %(outfile)s """
    P.run()

###################################################################
###################################################################
###################################################################


@transform(buildAnnotations,
           suffix('.annotations.gz'),
           '_annotations.load')
def loadAnnotations(infile, outfile):
    '''load variant annotations into database'''

    tablename = P.toTable(outfile)

    statement = '''gunzip < %(infile)s
    |python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              --quick
              --map=gene_id:str 
              --index=gene_id 
              --table=%(tablename)s
              --map=base_qualities:text 
    > %(outfile)s
    '''
    P.run()

###################################################################
###################################################################
###################################################################


@merge(buildAnnotations, 'annotations.load')
def mergeAnnotations(infiles, outfile):
    '''load variant annotations into single database table'''

    tablename = P.toTable(outfile)
    outf = open('anno.txt', 'w')
    first = True
    for f in infiles:
        track = P.snip(os.path.basename(f), ".annotations.gz")
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

    statement = '''cat anno.txt |
                   python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                       --table=%(tablename)s 
                   > %(outfile)s'''
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
              --index=code
              --table=%(tablename)s
    > %(outfile)s
    '''
    P.run()

###################################################################
###################################################################
# Effects
###################################################################
###################################################################


@files([("variants.vcf.gz", "%s.effects.gz" % x, x) for x in TRACKS])
def buildEffects(infile, outfile, sample):
    """annotate snps with gene set."""

    to_cluster = True
    seleno = "seleno.list"
    transcripts = os.path.join(
        PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_geneset_cds_gtf"])

    statement = """python %(scriptsdir)s/snp2counts.py 
                       --genome-file=%(genome_dir)s/%(genome)s
                       --filename-vcf=%(infile)s
                       --input-format=vcf
                       --vcf-sample=%(sample)s
                       --module=transcript-effects 
                       --filename-seleno=%(seleno)s 
                       --filename-exons=%(transcripts)s 
                       --output-filename-pattern=%(outfile)s.%%s.gz
                       --log=%(outfile)s.log 
                   | gzip > %(outfile)s """
    P.run()

###################################################################
###################################################################
###################################################################


@transform(buildEffects, suffix(".effects.gz"), "_effects.load")
def loadEffects(infile, outfile):
    '''load transcript effects into tables.'''

    root = infile[:-len(".effects.gz")]

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --from-zipped \
              --index=transcript_id \
              --table=%(root)s_effects \
    < %(infile)s > %(outfile)s
    '''
    P.run()

    for suffix in ("cds", "intron", "splicing", "translation"):

        statement = '''
        gunzip < %(infile)s.%(suffix)s.gz
        | python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
        --allow-empty
        --index=transcript_id 
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
    outf = open('effects.txt', 'w')
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

    statement = '''cat effect.txt |
                   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
                       --index=transcript_id \
                       --table=%(tablename)s \
                   > %(outfile)s'''
    P.run()

    for suffix in ("cds", "intron", "splicing", "translation", "genes"):

        outf = open('effects.' + suffix + '.txt', 'w')
        first = True
        for f in infiles:
            track = P.snip(os.path.basename(f), ".effects.gz")
            statfile = f + "." + suffix + ".gz"
            print(statfile)
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

        statement = '''cat %(tmpfilename)s |
                       python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                           --allow-empty
                           --index=transcript_id 
                           --table=%(tablename)s_%(suffix)s 
                           --ignore-column=seq_na
                           --ignore-column=seq_aa
                       >> %(outfile)s'''
        P.run()

###################################################################
###################################################################
###################################################################


@transform(loadEffects, suffix("_effects.load"), "_effects_genes.load")
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
###################################################################
###################################################################


@merge(mergeEffects, "effects_genes.load")
def mergeEffectsPerGene(infile, outfile):
    '''summarize effects on a per-gene level.'''

    tablename = outfile[:-len(".load")]

    dbhandle = connect()

    statement = '''
    CREATE TABLE %(tablename)s AS
    SELECT DISTINCT 
           track,
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
    FROM annotations.transcript_info as i, effects AS e
    WHERE i.transcript_id = e.transcript_id
    GROUP BY i.gene_id, track
    ''' % locals()

    Database.executewait(
        dbhandle, "DROP TABLE IF EXISTS %(tablename)s" % locals())
    Database.executewait(dbhandle, statement)
    Database.executewait(
        dbhandle, "CREATE INDEX %(tablename)s_gene_id ON %(tablename)s (gene_id)" % locals())
    dbhandle.commit()

    P.touch(outfile)

###################################################################
###################################################################
# Polyphen
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
        print statement % locals()
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
              --index=snp_id 
              --index=track,transcript_id
              --index=contig,pos
              --index=protein_id
              --index=transcript_id
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
              --index=snp_id 
              --index=protein_id
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
              --index=gene_id 
              --map=code:str
              --table=%(table)s 
    > %(outfile)s
    '''
    P.run()

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


@follows(buildAnnotations, loadAnnotations, mergeAnnotations, summarizeAnnotations, loadAnnotationsSummary)
def annotations():
    pass


@follows(buildEffects, loadEffects, mergeEffects, summarizeEffectsPerGene)
def effects():
    pass


@follows(buildPolyphenInput, buildPolyphenFeatures, runPolyphen, loadPolyphen, loadPolyphenMap, analysePolyphen, loadPolyphenAnalysis)
def polyphen():
    pass


@follows(annotations, effects, polyphen)
def full():
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
