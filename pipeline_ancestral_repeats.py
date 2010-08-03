################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_ancestral_repeats.py 2876 2010-03-27 17:42:11Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
"""

:Author: Andreas Heger
:Release: $Id: pipeline_ancestral_repeats.py 2876 2010-03-27 17:42:11Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

This pipeline performs the following actions:

   * build pairwise genomic alignment from axt or maf files
   * define ancestral repeats
   * compute rates of ancestral repeats

The pairwise genomic alignment is stored in :term:`psl` format
with genome1 as the query and genome2 as the target.

.. note::
   The pipeline only works, if genome1 is the reference species
   in the maf files. This is a result of maf2Axt requiring that
   the strand of the reference species is always positive and 
   I have not figured out how to invert maf alignments.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----


"""
import sys, tempfile, optparse, shutil, itertools, csv, math, random, re, glob, os, shutil, collections

import Experiment as E
import Pipeline as P
from ruffus import *
import csv
import sqlite3
import IndexedFasta, IndexedGenome, FastaIterator, Genomics
import GFF, GTF, Blat
import IOTools

if not os.path.exists("conf.py"):
    raise IOError( "could not find configuration file conf.py" )

execfile("conf.py")

PARAMS = P.getParameters()

def getGenomes():
    '''return genome names of query and target.'''

    genome_query = PARAMS["%s_genome" % PARAMS["query"]]
    genome_target = PARAMS["%s_genome" % PARAMS["target"]]
    return genome_query, genome_target
        
#########################################################################
#########################################################################
#########################################################################
@transform( "*.idx",
            suffix(".idx"), 
            ".sizes" )
def buildSizes( infile, outfile ):
    '''extract size information from genomes.'''
    outf = open(outfile, "w")
    for line in open(infile):
        data = line[:-1].split( "\t" )
        if len(data) >= 4:
            contig = data[0]
            if contig.startswith("chr"): contig = contig[3:]
            outf.write("%s\t%s\n" % (contig, data[3]))
    outf.close()

#########################################################################
#########################################################################
#########################################################################
if PARAMS["axt_dir"]:
    ## build pairwise alignment from axt formatted data.'''
    @follows( buildSizes )
    @merge( "%s/*.axt.gz" % PARAMS["axt_dir"], "alignment.psl" )
    def buildGenomeAlignment(infiles, outfile):
        '''build pairwise genomic aligment from axt files.'''
    
        try:
            os.remove( outfile )
        except OSError:
            pass

        genome_query, genome_target = getGenomes()

        for infile in infiles:
            E.info( "adding %s" % infile )
            statement = '''gunzip < %(infile)s |
                           axtToPsl 
                               /dev/stdin
                               %(genome_query)s.sizes 
                               %(genome_target)s.sizes 
                               /dev/stdout |
                           pslSwap /dev/stdin /dev/stdout 
                           >> %(outfile)s
                           '''

            P.run( **locals() )

elif PARAMS["maf_dir"]:
    @follows( buildSizes )
    @merge( "%s/*.maf.gz" % PARAMS["maf_dir"], "alignment.raw.psl.gz" )
    def buildRawGenomeAlignment(infiles, outfile):
        '''build pairwise genomic aligment from maf files.'''
    
        try: os.remove( outfile )
        except OSError: pass

        for infile in infiles:
            # skip maf files without Hsap on top.
            if "other" in infile or "supercontig" in infile: continue

            E.info( "adding %s" % infile )

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

    @transform( buildRawGenomeAlignment, 
            suffix(".raw.psl.gz"),
            ".psl.gz" )
    def buildGenomeAlignment( infile, outfile ):
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

    @follows( buildSizes )
    @merge( "%s/*.maf.gz" % PARAMS["maf_dir"], "alignment.psl.gz" )
    def buildGenomeAlignmentUCSCTools(infiles, outfile):
        '''build pairwise genomic aligment from maf files.'''
    
        try: os.remove( outfile )
        except OSError: pass

        for infile in infiles:
            # skip maf files without Hsap on top.
            if "other" in infile or "supercontig" in infile: continue

            E.info( "adding %s" % infile )

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
                  %(genome_target)s.sizes 
                  %(genome_query)s.sizes 
                  /dev/stdout 
             | python %(scriptsdir)s/psl2psl.py 
                  --filename-queries=%(genome_query)s
                  --filename-target=%(genome_target)s
                  --method=sanitize
             | gzip 
             >> %(outfile)s
             '''
            P.run( **locals() )

#########################################################################
#########################################################################
#########################################################################
def importRepeatsFromUCSC( infile, outfile, ucsc_database, repeattypes, genome ):
    '''import repeats from a UCSC formatted file.

    The repeats are stored as a :term:`gff` formatted file.
    '''

    repclasses="','".join(repeattypes.split(","))

    # Repeats are either stored in a single ``rmsk`` table (hg19) or in
    # individual ``rmsk`` tables (mm9) like chr1_rmsk, chr2_rmsk, ....
    # In order to do a single statement, the ucsc mysql database is 
    # queried for tables that end in rmsk.

    import MySQLdb
    dbhandle = MySQLdb.Connect( host = PARAMS["ucsc_host"],
                                user = PARAMS["ucsc_user"] )

    cc = dbhandle.cursor()
    cc.execute( "USE %s " %  ucsc_database )
    


    cc = dbhandle.cursor()
    cc.execute("SHOW TABLES LIKE '%rmsk'")
    tables = [ x[0] for x in cc.fetchall()]
    if len(tables) == 0:
        raise ValueError( "could not find any `rmsk` tables" )

    tmpfile = P.getTempFile()
    
    for table in tables:
        E.info( "loading repeats from %s" % table )
        cc = dbhandle.cursor()
        cc.execute("""SELECT genoName, 'repeat', 'exon', genoStart+1, genoEnd, strand, '.', '.', 
                      CONCAT('class \\"', repClass, '\\"; family \\"', repFamily, '\\";')
               FROM %(table)s
               WHERE repClass in ('%(repclasses)s') """ % locals() )
        for data in cc.fetchall():
            tmpfile.write( "\t".join(map(str,data)) + "\n" )

    tmpfile.close()
    tmpfilename = tmpfile.name

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
    
    os.unlink( tmpfilename )

#########################################################################
#########################################################################
########################################################################
def importRepeatsFromEnsembl( infile, outfile, ensembl_database, repeattypes, genome ):
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
@transform( "*_genome.fasta", suffix("_genome.fasta"), "_repeats.gff.gz")
def importRepeats( infile, outfile ):
    
    track = infile[:-len("_genome.fasta")]

    source = PARAMS["%s_source" % track] 
    if source == "ensembl":
        importRepeatsFromEnsembl( infile, outfile, 
                                  PARAMS["%s_database" % track],
                                  repeattypes = PARAMS["%s_repeattypes" % track],
                                  genome = "%s_genome" % track )                                  
    elif source == "ucsc":
        importRepeatsFromUCSC( infile, outfile, 
                               PARAMS["%s_database" % track],
                               repeattypes = PARAMS["%s_repeattypes" % track],
                               genome = "%s_genome" % track )                                  
        
#########################################################################
#########################################################################
########################################################################
@transform( importRepeats,
            suffix("_repeats.gff.gz"), 
            "_merged.gff.gz") 
def mergeRepeats( infile, outfile):
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
@follows( buildGenomeAlignment )
@merge( mergeRepeats, "aligned_repeats.psl.gz" )
def buildAlignedRepeats( infiles, outfile ):
    '''build alignment between repeats.
    '''
    
    infile_target = PARAMS["target"] + "_merged.gff.gz"
    infile_query = PARAMS["query"] + "_merged.gff.gz" 

    to_cluster = False
    # need to escape pipe symbols within farm.py command
    statement = r'''
        gunzip < alignment.psl.gz 
        | %(cmd-farm)s --split-at-lines=100 --log=%(outfile)s.log --binary 
             "python %(scriptsdir)s/psl2psl.py 
	        --method=test 
		--log=%(outfile)s.log 
	      | python %(scriptsdir)s/psl2psl.py 
		--method=map 
		--filter-query=<(gunzip < %(infile_query)s )
		--filter-target=<(gunzip < %(infile_target)s )
		--log=%(outfile)s.log " 
         | gzip 
         > %(outfile)s'''
    P.run()

########################################################
########################################################
########################################################
@files( buildAlignedRepeats, "aligned_repeats.rates.gz" )
def buildRepeatsRates( infile, outfile ):
    '''compute rates for individual aligned repeats.'''
    # path ruffus bug:
    infile = infile[0]
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

@transform( (buildAlignedRepeats, buildGenomeAlignment), 
   suffix(".psl.gz"),
   ".stats")
def computeAlignmentStats( infile, outfile ):
    '''compute alignment coverage statistics'''
    
    to_cluster = True

    statement = '''
    gunzip < %(infile)s |
    python %(scriptsdir)s/psl2stats.py 
        --log=%(outfile)s.log 
    > %(outfile)s'''
    
    P.run()

########################################################
########################################################
########################################################
@transform( mergeRepeats, suffix(".gff.gz"), ".stats")
def computeRepeatsCounts( infile, outfile ):
    '''count number and type of repeats.'''
    pass

# %_repeats_counts.stats: ucsc_%_repeats.table.gz
# 	$(PRELOG)
# 	@gunzip < $< | pe "s/#//" |\
# 	csv_cut genoName genoStart genoEnd repName repClass repFamily |\
# 	awk '/genoName/ {printf("%s\t%s\n", $$5, "length"); next;} {printf("%s\t%i\n", $$5, $$3-$$2); } ' |\
# 	t2t --group=1 --group-function=stats > $@
# 	$(EPILOG)

########################################################
########################################################
########################################################
@transform( mergeRepeats,
            suffix( "_merged.gff.gz" ),
            "_repeats_sizes.stats" )
def buildRepeatDistribution( infile, outfile ):
    '''count size and distance distribution of repeats.'''

    to_cluster = True

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
@files( buildRepeatsRates, "%s_rates.gff.gz" % PARAMS["query"])
def exportRatesAsGFF( infile, outfile ):
    '''export gff file with rate as score.'''

    to_cluster = True
    statement = '''gunzip
    < %(infile)s 
    | python %(toolsdir)s/csv_cut.py qName qStart qEnd distance converged 
    | awk '!/qName/ && $5 {printf("%%s\\tancestral_repeat\\texon\\t%%s\\t%%s\\t%%s\\t+\\t.\\t.\\n", $1, $2, $3, $4);}' 
    | gzip
    > %(outfile)s
    '''
    P.run()

@follows( importRepeats,
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
    
if __name__== "__main__":
    sys.exit( P.main(sys.argv) )

