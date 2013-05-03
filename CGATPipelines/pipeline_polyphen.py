################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_snps.py 2870 2010-03-03 10:20:29Z andreas $
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
:Release: $Id: pipeline_snps.py 2870 2010-03-03 10:20:29Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

Polyphen prediction pipeline.

Input:

Indels in pileup format.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

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
import types
import collections
import optparse
import shutil
import sqlite3
import CGAT.GFF as GFF
import CGAT.GTF as GTF
import CGAT.Experiment as E
import CGAT.Pipeline as P
import CGAT.IOTools as IOTools
import CGAT.Genomics as Genomics
import CGAT.Database as Database
import CGAT.FastaIterator as FastaIterator
import PipelineGeneset as PGeneset
import CGAT.Stats as Stats
import pysam

###################################################################
###################################################################
###################################################################
## read global options from configuration file
P.getParameters( 
    defaults = { 'polyphen_modes' : "" } )

P.PARAMS.update( 
    { "transcripts" :"transcripts.gtf.gz",
      "genes" : 'genes.gtf.gz',
      "annotation" : 'geneset_regions.gff.gz',
      "peptides": 'peptides.fasta',
      "cdna": 'cdna.fasta',
      "cds": 'cds.fasta' } )

PARAMS = P.PARAMS

if os.path.exists("conf.py"):
    execfile("conf.py")

SEPARATOR = "|"

###################################################################
###################################################################
###################################################################
## gene set section
############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], PARAMS['annotation'] )
def buildGeneRegions( infile, outfile ):
    '''annotate genomic regions with reference gene set.

    Only considers protein coding genes. In case of overlapping
    genes, only take the longest (in genomic coordinates).
    Genes not on UCSC contigs are removed.
    '''
    PGeneset.buildGeneRegions( infile, outfile )

############################################################
############################################################
############################################################
@follows( buildGeneRegions )
@files( PARAMS["ensembl_filename_gtf"], PARAMS['genes'] )
def buildGenes( infile, outfile ):
    '''build a collection of exons from the protein-coding
    section of the ENSEMBL gene set. The exons include both CDS
    and UTR.

    The set is filtered in the same way as in :meth:`buildGeneRegions`.
    '''
    PGeneset.buildProteinCodingGenes( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], "gene_info.load" )
def loadGeneInformation( infile, outfile ):
    '''load the transcript set.'''
    PGeneset.loadGeneInformation( infile, outfile )

############################################################
############################################################
############################################################
@files( buildGenes, "gene_stats.load" )
def loadGeneStats( infile, outfile ):
    '''load the transcript set.'''

    PGeneset.loadGeneStats( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], PARAMS["transcripts"] )
def buildTranscripts( infile, outfile ):
    '''build a collection of transcripts from the protein-coding
    section of the ENSEMBL gene set.

    Only CDS are used.
    '''
    PGeneset.buildProteinCodingTranscripts( infile, outfile )

############################################################
############################################################
############################################################
@transform( buildTranscripts, suffix(".gtf.gz"), "_gtf.load" )
def loadTranscripts( infile, outfile ):
    '''load the transcript set.'''
    PGeneset.loadTranscripts( infile, outfile )

############################################################
############################################################
############################################################
@files( buildTranscripts, "transcript_stats.load" )
def loadTranscriptStats( infile, outfile ):
    '''load the transcript set.'''

    PGeneset.loadTranscriptStats( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], "transcript_info.load" )
def loadTranscriptInformation( infile, outfile ):
    '''load the transcript set.'''
    PGeneset.loadTranscriptInformation( infile, 
                                          outfile,
                                          only_proteincoding = PARAMS["ensembl_only_proteincoding"] )

###################################################################
###################################################################
###################################################################
@files( ( (PARAMS["ensembl_filename_pep"], PARAMS["peptides"] ), ) )
def buildPeptideFasta( infile, outfile ):
    '''load ENSEMBL peptide file
    
    *infile* is an ENSEMBL .pep.all.fa.gz file.
    '''
    PGeneset.buildPeptideFasta( infile, outfile )

###################################################################
###################################################################
###################################################################
@files( ( (PARAMS["ensembl_filename_cdna"], PARAMS["cdna"] ), ) )
def buildCDNAFasta( infile, outfile ):
    '''load ENSEMBL peptide file
    
    *infile* is an ENSEMBL .cdna.all.fa.gz file.
    '''
    PGeneset.buildCDNAFasta( infile, outfile )

###################################################################
###################################################################
###################################################################
@follows( loadTranscriptInformation )
@files( [(PARAMS["transcripts"], PARAMS["cds"] ),] )
def buildCDSFasta( infile, outfile ):
    '''build cds sequences from peptide and cds file.
    
    *infile* is an ENSEMBL .cdna.all.fa.gz file.
    '''

    PGeneset.buildCDSFasta( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_pep"], "protein_stats.load" )
def loadProteinStats( infile, outfile ):
    '''load the transcript set.'''

    PGeneset.loadProteinStats( infile, outfile )

############################################################
############################################################
############################################################
@files( ((None, "benchmark.ids"), ) )
def buildBenchmarkSet( infile, outfile ):
    '''build a benchmark set of protein ids.'''
    pass

############################################################
############################################################
############################################################
@files( ((buildBenchmarkSet, "benchmark.input") ,) )
def buildBenchmarkInput( infile, outfile ):

    tmpfile = P.getTempFile()

    dbhandle = sqlite3.connect( PARAMS["database"] )
    cc = dbhandle.cursor()
    statement = '''
    SELECT DISTINCT transcript_id, protein_id FROM peptide_info
    '''
    cc.execute( statement )
    tmpfile.write("transcript_id\tprotein_id\n")
    tmpfile.write("\n".join( [ "\t".join(x) for x in cc ] ) )
    tmpfile.write("\n")
    tmpfilename = tmpfile.name

    statement = '''
    perl %(scriptsdir)s/extract_fasta.pl %(infile)s
    < cds.fasta 
    python %(scripstdir)s/fasta2variants.py --is-cds  
    | python %(scriptsdir)s/substitute_tokens.py 
             --apply=%(tmpfilename)s
    > %(outfile)s
    '''
    P.run()

    os.unlink( tmpfilename )

###################################################################
###################################################################
###################################################################
@transform( "*.input", suffix(".input"), ".features")
def buildPolyphenFeatures( infile, outfile ):
    '''run polyphen on the cluster.

    To do this, first send uniref to all nodes:

    python ~/cgat/cluster_distribute.py 
           --collection=andreas 
           /net/cpp-group/tools/polyphen-2.0.18/nrdb/uniref100*.{pin,psd,psi,phr,psq,pal}
    '''
    
    nsnps = len([ x for x in open(infile)])

    to_cluster = True
    stepsize = max( int(nsnps / 200000.0), 1000 )
    job_array=(0, nsnps, stepsize)
    E.info("running array jobs on %i snps" % nsnps )

    scratchdir = os.path.join(os.path.abspath("."), "scratch")
    try:
        os.mkdir( scratchdir )
    except OSError:
        pass

    resultsdir = outfile + ".dir"
    try:
        os.mkdir( resultsdir )
    except OSError:
        pass

    statement = '''
    /net/cpp-group/tools/polyphen-2.0.18/bin/run_pph_cpp.pl
       -s %(peptides)s
       -b %(polyphen_blastdb)s
       -d %(scratchdir)s
       %(infile)s > %(resultsdir)s/%(outfile)s.$SGE_TASK_ID 2> %(resultsdir)s/%(outfile)s.err.$SGE_TASK_ID
    '''
    P.run()

    to_cluster = False
    job_array=None

    statement = '''find %(resultsdir)s -name "*.err.*" -exec cat {} \; 
                   | gzip 
                   > %(outfile)s.log.gz'''
    P.run()

    statement = '''find %(resultsdir)s -not -name "*.err.*" -exec cat {} \; 
                | gzip 
                > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################
@files( [ (x, "%s_%s.output.gz" % (x[:-len(".features.gz")],y), y ) \
              for x,y in itertools.product( 
            glob.glob( "*.features.gz"), P.asList( PARAMS["polyphen_models"] ) ) ] ) 
def runPolyphen( infile, outfile, model ):
    '''run POLYPHEN on feature tables to classify SNPs.
    '''

    to_cluster = False
    
    # need to run in chunks for large feature files
    statement = """gunzip 
        < %(infile)s
        | %(cmd-farm)s
            --split-at-lines=10000
            --output-header
        "perl %(polyphen_home)s/bin/run_weka_cpp.pl 
           -l %(polyphen_home)s/models/%(model)s.UniRef100.NBd.f11.model
           -p 
           %%STDIN%%"
        | gzip > %(outfile)s 
    """

    P.run()

    return

###################################################################
###################################################################
###################################################################
@transform( buildBenchmarkInput, suffix(".input"), ".load" )
def loadPolyphenWeights( infile, outfile ):
    '''load polyphen input data.'''

    table = "weights"

    statement = '''
    cat < %(infile)s 
    | python %(scriptsdir)s/csv_cut.py snpid counts weight 
    |python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --index=snp_id 
              --table=%(table)s 
    > %(outfile)s
    '''
    P.run()

###################################################################
###################################################################
###################################################################
@transform( runPolyphen, suffix(".output.gz"), ".load")
def loadPolyphen( infile, outfile ):
    '''load polyphen results.

    The comment column is ignored.
    '''
    
    table = P.toTable( outfile )

    statement = '''gunzip
    < %(infile)s
    | perl -p -e "s/o_acc/protein_id/; s/ +//g"
    | cut -f 1-55
    |python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --index=snp_id 
              --index=protein_id
              --table=%(table)s 
              --map=effect:str
    > %(outfile)s
    '''
    P.run()


@follows( loadTranscripts,
          loadTranscriptInformation,
          loadGeneStats,
          loadGeneInformation,
          buildPeptideFasta,
          buildCDSFasta)
def prepare():
    pass

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )

