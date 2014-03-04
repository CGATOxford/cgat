"""
===========================
Ancestral repeats pipeline
===========================

:Author: Andreas Heger
:Release: $Id: pipeline_ancestral_repeats.py 2876 2010-03-27 17:42:11Z andreas $
:Date: |today|
:Tags: Python

The ancestral repeats pipeline defines ancestral repeats for a pair of genomes
and computes rates for these.

This pipeline performs the following actions:

   * collect repeatmasker annotation from external databases. Currently implemented are:
      * UCSC
      * Ensembl
   * build pairwise genomic alignment from axt or maf files
   * define ancestral repeats
   * compute rates of ancestral repeats


Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline expects a :term:`query` and :term:`target` genome. These should be set in the general section.
For each genome there should then be section on how to obtain the repeatmasker tracks. The default
configuration file gives an example.

Input
-----

The pipeline starts from an empty working directory. It will collect the input
data from directories specified in the configuration files.

The genomic alignment can both be build from :term:`axt` formatted pairwise alignments
and from :term:`maf` formatted multiple alignmentns. However, the latter currently 
only works if the :term:`query` genome is the reference species in the maf files. 

This is a consequence of :file:`maf2Axt` requiring that the strand of the reference species 
is always positive and I have not figured out how to invert maf alignments.

.. note::
   ENSEMBL import is not thoroughly tested.
   :term:`maf` formatted import is not thoroughly tested.

Type::

   python pipeline_ancestral_repeats.py --help

for command line help.

Requirements
------------



Output
======

The pipeline builds the following files:

aligned_repeats.psl.gz
   :term:`psl` formatted files of alignments between ancestral repeats
 
aligned_repeats.rates.gz
   rates between ancestral repeats

alignment.psl.gz
   :term:`psl` formatted genomic alignment between query and target.

<query>_rates.gff.gz
   :term:`gff` formatted file of ancestral repeats on the query. The score field is set
   to the estimated substitution rate of the repeat.

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_ancestral_repeats.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_ancestral_repeats.tgz
   tar -xvzf pipeline_ancestral_repeats.tgz
   cd pipeline_ancestral_repeats
   python <srcdir>/pipeline_ancestral_repeats.py make full

The example data builds ancestral repeats between human hg19:chr19 and mouse mm9:chr7.

Code
====


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

import CGAT.Experiment as E
import logging as L

from ruffus import *
import csv
import sqlite3
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IndexedGenome as IndexedGenome
import CGAT.FastaIterator as FastaIterator
import CGAT.Genomics as Genomics
import CGAT.GTF as GTF
import CGAT.Blat as Blat
import CGAT.IOTools as IOTools

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################
import CGAT.Pipeline as P
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"],
    defaults={
        'query': "",
        'target': ""})
PARAMS = P.PARAMS
USECLUSTER = True

#########################################################################
#########################################################################
#########################################################################
if os.path.exists("pipeline_conf.py"):
    L.info("reading additional configuration from pipeline_conf.py")
    execfile("pipeline_conf.py")


def getGenomes():
    '''return genome names of query and target.'''

    genome_query = os.path.join(PARAMS["genome_dir"], PARAMS["query"])
    genome_target = os.path.join(PARAMS["genome_dir"], PARAMS["target"])
    return genome_query, genome_target

#########################################################################
#########################################################################
#########################################################################


@files([("%s/%s.idx" % (PARAMS["genome_dir"], x), "%s.sizes" % x) for x in (PARAMS["query"], PARAMS["target"])])
def buildSizes(infile, outfile):
    '''extract size information from genomes.'''
    outf = open(outfile, "w")
    for line in open(infile):
        data = line[:-1].split("\t")
        if len(data) >= 4:
            contig = data[0]
            outf.write("%s\t%s\n" % (contig, data[3]))
    outf.close()

#########################################################################
#########################################################################
#########################################################################
if "axt_dir" in PARAMS:
    # build pairwise alignment from axt formatted data.'''
    @follows(buildSizes)
    @merge("%s/*.axt.gz" % PARAMS["axt_dir"], PARAMS["interface_alignment_psl"])
    def buildGenomeAlignment(infiles, outfile):
        '''build pairwise genomic aligment from axt files.'''

        to_cluster = USECLUSTER

        try:
            os.remove(outfile)
        except OSError:
            pass

        for infile in infiles:
            E.info("adding %s" % infile)
            statement = '''gunzip < %(infile)s 
                           | axtToPsl 
                               /dev/stdin
                               %(query)s.sizes 
                               %(target)s.sizes 
                               /dev/stdout 
                           | pslSwap /dev/stdin /dev/stdout 
                           | gzip >> %(outfile)s
                           '''
            P.run()


elif "maf_dir" in PARAMS:
    @follows(buildSizes)
    @merge("%s/*.maf.gz" % PARAMS["maf_dir"], "alignment.raw.psl.gz")
    def buildRawGenomeAlignment(infiles, outfile):
        '''build pairwise genomic aligment from maf files.
        '''

        try:
            os.remove(outfile)
        except OSError:
            pass

        to_cluster = USECLUSTER

        for infile in infiles:
            # skip maf files without Hsap on top.
            if "other" in infile or "supercontig" in infile:
                continue

            E.info("adding %s" % infile)

            genome_query, genome_target = getGenomes()

            statement = '''gunzip < %(infile)s 
             | python %(scriptsdir)s/maf2psl.py 
                  --query=%(maf_name_query)s
                  --target=%(maf_name_target)s
                  --log=%(outfile)s.log 
             | python %(scriptsdir)s/psl2psl.py 
                  --method=filter-fasta 
                  --method=sanitize
                  --filename-queries=%(genome_query)s
                  --filename-target=%(genome_target)s
                  --log=%(outfile)s.log 
             | gzip 
             >> %(outfile)s
             '''
            P.run()

    @transform(buildRawGenomeAlignment,
               suffix(".raw.psl.gz"),
               ".psl.gz")
    def buildGenomeAlignment(infile, outfile):
        '''remove non-unique alignments in genomic infile.'''

        to_cluster = True

        statement = '''gunzip < %(infile)s 
             | sort -k10,10 -k12,12n
             | python %(scriptsdir)s/psl2psl.py 
                  --method=remove-overlapping-query
                  --log=%(outfile)s.log 
             | sort -k14,14 -k16,16n
             | python %(scriptsdir)s/psl2psl.py 
                  --method=remove-overlapping-target
                  --log=%(outfile)s.log 
             | gzip
             >> %(outfile)s
             '''
        P.run()

    @follows(buildSizes)
    @merge("%s/*.maf.gz" % PARAMS["maf_dir"], PARAMS["interface_alignment_psl"])
    def buildGenomeAlignmentUCSCTools(infiles, outfile):
        '''build pairwise genomic aligment from maf files.'''

        try:
            os.remove(outfile)
        except OSError:
            pass

        to_cluster = USECLUSTER

        for infile in infiles:
            # skip maf files without Hsap on top.
            if "other" in infile or "supercontig" in infile:
                continue

            E.info("adding %s" % infile)

            genome_query, genome_target = getGenomes()

            statement = '''gunzip < %(infile)s 
            | mafToAxt
                  /dev/stdin
                  %(maf_name_target)s
                  %(maf_name_query)s
                  /dev/stdout 
                  -stripDb 
             | axtToPsl 
                  /dev/stdin 
                  %(target)s.sizes 
                  %(query)s.sizes 
                  /dev/stdout 
             | python %(scriptsdir)s/psl2psl.py 
                  --filename-queries=%(genome_query)s
                  --filename-target=%(genome_target)s
                  --method=sanitize
             | gzip 
             >> %(outfile)s
             '''
            P.run()
else:
    raise ValueError(
        "configuration error: please specify either maf_dir or axt_dir")

#########################################################################
#########################################################################
#########################################################################


def importRepeatsFromUCSC(infile, outfile, ucsc_database, repeattypes, genome):
    '''import repeats from a UCSC formatted file.

    The repeats are stored as a :term:`gff` formatted file.
    '''

    repclasses = "','".join(repeattypes.split(","))

    # Repeats are either stored in a single ``rmsk`` table (hg19) or in
    # individual ``rmsk`` tables (mm9) like chr1_rmsk, chr2_rmsk, ....
    # In order to do a single statement, the ucsc mysql database is
    # queried for tables that end in rmsk.

    import MySQLdb
    dbhandle = MySQLdb.Connect(host=PARAMS["ucsc_host"],
                               user=PARAMS["ucsc_user"])

    cc = dbhandle.cursor()
    cc.execute("USE %s " % ucsc_database)

    cc = dbhandle.cursor()
    cc.execute("SHOW TABLES LIKE '%rmsk'")
    tables = [x[0] for x in cc.fetchall()]
    if len(tables) == 0:
        raise ValueError("could not find any `rmsk` tables")

    tmpfile = P.getTempFile(".")

    for table in tables:
        E.info("loading repeats from %s" % table)
        cc = dbhandle.cursor()
        cc.execute("""SELECT genoName, 'repeat', 'exon', genoStart+1, genoEnd, '.', strand, '.', 
                      CONCAT('class \\"', repClass, '\\"; family \\"', repFamily, '\\";')
               FROM %(table)s
               WHERE repClass in ('%(repclasses)s') """ % locals() )
        for data in cc.fetchall():
            tmpfile.write("\t".join(map(str, data)) + "\n")

    tmpfile.close()
    tmpfilename = tmpfile.name

    to_cluster = USECLUSTER

    statement = '''cat %(tmpfilename)s
        | %(scriptsdir)s/gff_sort pos 
        | python %(scriptsdir)s/gff2gff.py 
            --sanitize=genome 
            --skip-missing 
            --genome-file=%(genome)s
            --log=%(outfile)s.log 
        | gzip
        > %(outfile)s
    '''
    P.run()

    os.unlink(tmpfilename)

#########################################################################
#########################################################################
########################################################################


def importRepeatsFromEnsembl(infile, outfile, ensembl_database, repeattypes, genome):
    '''import repeats from an ENSEMBL database.
    '''
    statement = '''
        perl %(scriptsdir)s/ensembl_repeats2gff.pl 
              -h %(ensembl_host)s 
              -u %(ensembl_user)s
              -p %(ensembl_password)s
              -d %(ensembl_database)s
              --repeattypes %(repeattypes)s 
	| %(scriptsdir)s/gff_sort pos 
        | python %(scriptsdir)s/gff2gff.py 
            --sanitize=genome
            --skip-missing 
            --genome-file=%(genome)s
            --log=%(outfile)s.log 
        | gzip
        > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
########################################################################


@files([(None, "%s_repeats.gff.gz" % x, x) for x in (PARAMS["query"], PARAMS["target"])])
def importRepeats(infile, outfile, track):
    '''import repeats from external sources.'''

    source = PARAMS["%s_source" % track]
    genome = os.path.join(PARAMS["genome_dir"], track)

    if source == "ensembl":
        importRepeatsFromEnsembl(infile, outfile,
                                 PARAMS["%s_database" % track],
                                 repeattypes=PARAMS["%s_repeattypes" % track],
                                 genome=genome)
    elif source == "ucsc":
        importRepeatsFromUCSC(infile, outfile,
                              PARAMS["%s_database" % track],
                              repeattypes=PARAMS["%s_repeattypes" % track],
                              genome=genome)

#########################################################################
#########################################################################
########################################################################


@transform(importRepeats,
           suffix("_repeats.gff.gz"),
           "_merged.gff.gz")
def mergeRepeats(infile, outfile):
    '''merge adjacent repeats.'''

    to_cluster = True

    statement = '''gunzip
    < %(infile)s 
    | python %(scriptsdir)s/gff2gff.py 
            --merge-features=0,10,0,0 
            --log=%(outfile)s.log 
    | gzip
    > %(outfile)s
    '''
    P.run()

########################################################
########################################################
########################################################


@follows(buildGenomeAlignment)
@merge(mergeRepeats, "aligned_repeats.psl.gz")
def buildAlignedRepeats(infiles, outfile):
    '''build alignment between repeats.
    '''

    infile_target = PARAMS["target"] + "_merged.gff.gz"
    infile_query = PARAMS["query"] + "_merged.gff.gz"

    # using farm.py to send to cluster
    # granularity should be set automatically.
    granularity = 5000

    # need to escape pipe symbols within farm.py command
    #to_cluster = False
    # statement = r'''
    #     gunzip < %(interface_alignment_psl)s
    #     | %(cmd-farm)s --split-at-lines=%(granularity)i --log=%(outfile)s.log --binary
    #          "python %(scriptsdir)s/psl2psl.py
    #             --method=test
    #     	--log=%(outfile)s.log
    #           | python %(scriptsdir)s/psl2psl.py
    #     	--method=map
    #     	--filter-query=%(infile_query)s
    #     	--filter-target=%(infile_target)s
    #     	--log=%(outfile)s.log "
    #      | gzip
    #      > %(outfile)s'''
    # P.run()

    statement = '''
        gunzip < %(interface_alignment_psl)s
        | python %(scriptsdir)s/psl2psl.py 
	        --method=test 
		--log=%(outfile)s.log 
	| python %(scriptsdir)s/psl2psl.py 
		--method=map 
		--filter-query=%(infile_query)s
		--filter-target=%(infile_target)s
		--log=%(outfile)s.log
         | gzip 
         > %(outfile)s'''
    P.run()

########################################################
########################################################
########################################################


@files(buildAlignedRepeats, "aligned_repeats.rates.gz")
def buildRepeatsRates(infile, outfile):
    '''compute rates for individual aligned repeats.'''

    to_cluster = False
    genome_query, genome_target = getGenomes()

    statement = '''gunzip < %(infile)s |
    sort -k10,10 -k14,14 -k9,9 -k12,12n |
    %(cmd-farm)s --split-at-lines=10000 --output-header --log=%(outfile)s.log
          "python %(scriptsdir)s/psl2psl.py 
		--log=%(outfile)s.log 
		--method=add-sequence 
		--filename-queries=%(genome_query)s
		--filename-target=%(genome_target)s |
	   python %(scriptsdir)s/psl2table.py 
                --method=query-counts 
                --method=baseml 
                --baseml-model=REV" |
    gzip > %(outfile)s
    '''
    P.run()


@transform((buildAlignedRepeats, buildGenomeAlignment),
           suffix(".psl.gz"),
           ".stats")
def computeAlignmentStats(infile, outfile):
    '''compute alignment coverage statistics'''

    to_cluster = USECLUSTER

    statement = '''
    gunzip < %(infile)s |
    python %(scriptsdir)s/psl2stats.py 
        --log=%(outfile)s.log 
    > %(outfile)s'''

    P.run()

########################################################
########################################################
########################################################


@transform(mergeRepeats, suffix(".gff.gz"), ".stats")
def computeRepeatsCounts(infile, outfile):
    '''count number and type of repeats.'''
    pass

# %_repeats_counts.stats: ucsc_%_repeats.table.gz
# 	$(PRELOG)
# @gunzip < $< | pe "s/#//" |\
# 	csv_cut genoName genoStart genoEnd repName repClass repFamily |\
# 	awk '/genoName/ {printf("%s\t%s\n", $$5, "length"); next;} {printf("%s\t%i\n", $$5, $$3-$$2); } ' |\
# 	t2t --group=1 --group-function=stats > $@
# 	$(EPILOG)

########################################################
########################################################
########################################################


@transform(mergeRepeats,
           suffix("_merged.gff.gz"),
           "_repeats_sizes.stats")
def buildRepeatDistribution(infile, outfile):
    '''count size and distance distribution of repeats.'''

    to_cluster = USECLUSTER

    statement = '''gunzip
    < %(infile)s 
    | python %(scriptsdir)s/gff2histogram.py 
        --output-filename-pattern="%(outfile)s.%%s" 
        --method=all 
    > %(outfile)s
    '''
    P.run()

########################################################
########################################################
########################################################


@files(buildRepeatsRates, PARAMS["interface_rates_query_gff"])
def exportRatesAsGFF(infile, outfile):
    '''export gff file with rate as score.'''

    to_cluster = USECLUSTER
    statement = '''gunzip
    < %(infile)s 
    | python %(toolsdir)s/csv_cut.py qName qStart qEnd distance converged 
    | awk '!/qName/ && $5 {printf("%%s\\tancestral_repeat\\texon\\t%%s\\t%%s\\t%%s\\t+\\t.\\t.\\n", $1, $2, $3, $4);}' 
    | gzip
    > %(outfile)s
    '''
    P.run()


@follows(importRepeats,
         mergeRepeats,
         buildAlignedRepeats,
         buildRepeatsRates,
         buildRepeatDistribution,
         computeAlignmentStats,
         computeRepeatsCounts,
         exportRatesAsGFF,
         )
def full():
    pass

###################################################################
###################################################################
###################################################################
# primary targets
###################################################################


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.'''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish_report():
    '''publish report.'''

    E.info("publishing report")
    P.publish_report()

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
