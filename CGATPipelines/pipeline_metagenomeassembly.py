"""
=============================
Metagenome assembly pipeline
=============================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

The metagenome assembly pipeline takes reads from one or more NGS experiments and
assembles into contigs / scaffolds. Genes present on contigs are predicted using ORF
perdiction software.

Overview
========

The pipeline assumes the data derive from multiple tissues/conditions (:term:`experiment`) 
with one or more biological and/or technical replicates (:term:`replicate`). A :term:`replicate`
within each :term:`experiment` is a :term:`track`.

Assembly stategy
----------------

While there exist many tools for assembling reads from single genomes, only recently has 
software been specifically developed (or extended) to allow for the assembly of metagenomes.
The major factor affecting the ability to assemble a metagenome is the diversity of the sample.
Accurate assembly of long contigs is difficult in the presence of many species at differential
abundances. This is in constrast to single genome assembly where reads are (or should be) uniformly
sampled across the genome.

This pieline therefore uses a range of available software for the assembly of metagenomes.

Considerations
--------------

Metagenomics is a young and rapidly developing field. There is, as yet, no gold standard for 
assembly. It is likely that the presence of multiple, highly related species will lead to 
the assembly of schimeric contigs i.e. contigs derived from more than one species. It is 
generally considered that longer K-mer lengths used in the construction of the de-bruijn
graph (for de-bruijn graph assemblers) will result in fewer chimeras. Nevertheless, longer
k-mers may also result in more, short contigs being produced as a result of a neccessity for
a greater overlap between nodes in the graph. The length of k-mer chosen is also dependent on
the length of reads that you are trying to assemble - longer reads means you can use longer
k-mers. Which k-mer to use in the assembly process is therefore dependent on the data used
and the expected complexity of the sample. We make no effort here to advise on k-mer length.

  
Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.

Input
-----

Reads
+++++

Reads are imported by placing files are linking to files in the :term:`working directory`.

The default file format assumes the following convention:

   <sample>-<condition>-<replicate>.<suffix>

``sample`` and ``condition`` make up an :term:`experiment`, while ``replicate`` denotes
the :term:`replicate` within an :term:`experiment`. The ``suffix`` determines the file type.
The following suffixes/file types are possible:

fastq.gz
   Single-end reads in fastq format.

fastq.1.gz, fastq2.2.gz
   Paired-end reads in fastq format. The two fastq files must be sorted by read-pair.

.. note::

   Quality scores need to be of the same scale for all input files. Thus it might be
   difficult to mix different formats.

Optional inputs
+++++++++++++++

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|ray meta            |>=2.2.0            |Metagenome assembler                            |
+--------------------+-------------------+------------------------------------------------+
|meta-velvet         |>=1.2.02           |Metagenome assembler                            |
+--------------------+-------------------+------------------------------------------------+
|idba_ud             |>=1.1.0            |Metagenome assembler                            |
+--------------------+-------------------+------------------------------------------------+
|metaphlan           |>=1.7.7            |Taxon abundance estimator                       |
+--------------------+-------------------+------------------------------------------------+
|MetaGeneMark        |>=1.0.0            |ORF prediction                                  |
+--------------------+-------------------+------------------------------------------------+
|bedtools            |>=2.17.0           |BED interval analysis suite                     |
+--------------------+-------------------+------------------------------------------------+
|bwa                 |>=0.5.9            |Short read alignment algorithm                  |
+--------------------+-------------------+------------------------------------------------+
|bowtie              |>=0.12.7           |Short read alignment algorithm                  |
+--------------------+-------------------+------------------------------------------------+
|bowtie2             |>=2.0.0            |Short read alignment algorithm                  |
+--------------------+-------------------+------------------------------------------------+
|blastn              |>=2.2.25           |Simlilarity searching algorithm (nucleotides)   |
+--------------------+-------------------+------------------------------------------------+
|blastp              |>=2.2.25           |Simlilarity searching algorithm (proteins)      |
+--------------------+-------------------+------------------------------------------------+
|rpsblast            |>=2.2.25           |Simlilarity searching algorithm (profiles)      |
+--------------------+-------------------+------------------------------------------------+
|hmmer               |>=3                |gene annotation based on hmm models             |
+--------------------+-------------------+------------------------------------------------+



Pipeline output
===============

The main output is the genome assembly - output as a fasta formatted file.
Additional outputs include taxon abundance estimation (metaphlan) and ORF
predictions (MetaGeneMark).

Additional outputs are stored in the database file :file:`csvdb`.

Glossary
========

.. glossary::

Code
====

"""

# load modules
from ruffus import *

import CGAT.Experiment as E
import logging as L
import CGAT.Database as Database
import CGAT.CSV as CSV

import sys, os, re, shutil, itertools, math, glob, time, gzip, collections, random

import numpy, sqlite3
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMetagenomeAssembly as PipelineMetagenomeAssembly
import CGAT.FastaIterator as FastaIterator
import CGAT.Metaphlan as Metaphlan
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMappingQC as PipelineMappingQC
import pysam
import CGAT.Fastq as Fastq

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import CGAT.Pipeline as P
P.getParameters( 
    "pipeline.ini" )


PARAMS = P.PARAMS

###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
import PipelineTracks

# collect fastq.gz tracks
TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
        glob.glob( "*.fastq.gz" ), "(\S+).fastq.gz" ) +\
        PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
            glob.glob( "*.fastq.1.gz" ), "(\S+).fastq.1.gz" )
            
ALL = PipelineTracks.Sample3()
EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", "tissue" ) )
CONDITIONS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", ) )
TISSUES = PipelineTracks.Aggregate( TRACKS, labels = ("tissue", ) )

###################################################################
## Global flags
###################################################################
ASSEMBLERS = P.asList( PARAMS["assemblers"] )
MAPPER = PARAMS["coverage_mapper"]
BOWTIE = MAPPER == "bowtie"
BOWTIE2 = MAPPER == "bowtie2"
BWA = MAPPER == "bwa"

###################################################################
###################################################################
###################################################################
def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''
    dbh = sqlite3.connect( PARAMS["database"] )
    return dbh

###################################################################
###################################################################
###################################################################
# Should reads be pooled
###################################################################
###################################################################
###################################################################
SEQUENCEFILES = ("*.fasta", "*.fasta.gz", "*.fasta.1.gz"
                 , "*.fastq","*.fastq.gz", "*.fastq.1.gz")
SEQUENCEFILES_REGEX = regex(r"(\S+).(fasta$|fasta.gz|fasta.1.gz|fastq$|fastq.gz|fastq.1.gz)")

def pool_out(infiles):
    '''
    return outfile name dependent on
    input pairedness
    '''
    out = {"separate": "1",
                    False: ""}
    inf = infiles[0]
    paired = PipelineMetagenomeAssembly.PairedData().checkPairs(inf)
    if paired:
        paired = paired[0]
    format = PipelineMetagenomeAssembly.PairedData().getFormat(inf)
    outname = "pooled_reads.dir/agg-agg-agg.%s" % format
    return outname

############################################################
@active_if(PARAMS["pool_reads"])
@follows(mkdir("pooled_reads.dir"))
@merge(SEQUENCEFILES
       , pool_out([x for x in glob.glob("*R*.fast*") if not x.endswith(".2.gz") and not x.endswith(".2")])) # bit of a hack
def poolReadsAcrossConditions(infiles, outfile):
    '''
    pool reads across conditions
    '''
    statement = PipelineMetagenomeAssembly.pool_reads(infiles, outfile)
    P.run()

###################################################################
###################################################################
###################################################################
# Taxonomic profiling
###################################################################
###################################################################                                                                                                                                                                          
## load number of reads                                                                                                                                                                                                                      
###################################################################                                                                                                                                                                          
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX,
            r"\1.nreads" )
def countReads( infile, outfile ):
    '''count number of reads in input files.'''
    to_cluster = True
    m = PipelineMapping.Counter()
    statement = m.build( (infile,), outfile )
    P.run()

@merge(countReads, "reads_summary.load" )
def loadReadCounts( infiles, outfile ):
    '''load read counts into database.'''

    to_cluster = False
    outf = P.getTempFile()
    outf.write( "track\ttotal_reads\n")
    for infile in infiles:
        track = P.snip(infile, ".nreads")
        lines = IOTools.openFile( infile ).readlines()
        nreads = int( lines[0][:-1].split("\t")[1])
        outf.write( "%s\t%i\n" % (track,nreads))
    outf.close()
    inname = outf.name

    tablename = P.toTable(outfile)
    statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s --log=%(outfile)s.log 
                  < %(inname)s > %(outfile)s'''
    P.run()
    os.unlink(outf.name)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
## preprocess reads for metaphlan and IDBA
###################################################################                                                                                                                                                                          
@active_if("idba" in ASSEMBLERS and PARAMS["pool_reads"])
@transform(poolReadsAcrossConditions, regex("(\S+).fastq.*gz"), r"\1.fa")
def preprocessIdba(infile, outfile):
    '''
    preprocess pooled reads for IDBA
    '''
    # check for second read in the pair
    if infile.endswith(".fastq.gz"):
        E.info("converting fastq file to fasta file")
        outf = open(outfile, "w")
        for fastq in Fastq.iterate(IOTools.openFile(infile)):
            outf.write("%s\n%s\n" % (">" + fastq.identifier, fastq.seq))
        outf.close()
    elif infile.endswith(".1.gz"):
        read2 = P.snip(infile, ".1.gz") + ".2.gz"
        assert os.path.exists(read2), "file does not exist %s" % read2
    
        statement = '''python %(scriptsdir)s/fastqs2fasta.py 
                   -a %(infile)s 
                   -b %(read2)s 
                   --log=%(infile)s.log 
                   > %(outfile)s'''
        P.run()


###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(SEQUENCEFILES, SEQUENCEFILES_REGEX, r"\1.fa")
def preprocessReads(infile, outfile):
    '''
    create merged fasta file for use with metaphlan 
    and IDBA
    '''
    # check for second read in the pair
    if infile.endswith(".fastq.gz"):
        E.info("converting fastq file to fasta file")
        outf = open(outfile, "w")
        for fastq in Fastq.iterate(IOTools.openFile(infile)):
            outf.write("%s\n%s\n" % (">" + fastq.identifier, fastq.seq))
        outf.close()
    elif infile.endswith(".1.gz"):
        read2 = P.snip(infile, ".1.gz") + ".2.gz"
        assert os.path.exists(read2), "file does not exist %s" % read2
        
        log = infile.replace("fastq.","")
        statement = '''python %(scriptsdir)s/fastqs2fasta.py 
                   -a %(infile)s 
                   -b %(read2)s 
                   --log=%(log)s.log 
                   > %(outfile)s'''
        P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
## annotate metagenomic reads with metaphlan
###################################################################                                                                                                                                                                          
@follows(mkdir("metaphlan.dir"))
@transform(preprocessReads
           , regex("(\S+).fa")
            , r"metaphlan.dir/\1.readmap")
def buildMetaphlanReadmap(infile, outfile):
    '''
    metaphlan is a program used in metagenomics. It assigns
    reads to clades based on specific genetic markers via 
    blastn searching
    '''
    to_cluster = True

    # at present the pipeline will take a set of files
    # and compute the abundances of different taxonomic groups
    # based on ALL reads i.e. paired data are combined into
    # a single file for analysis
    if PARAMS["metaphlan_executable"] == "bowtie2":
        assert os.path.exists(PARAMS["metaphlan_db"] + ".1.bt2"), "missing file %s: Are you sure you have the correct database for bowtie2?" % PARAMS["metaphlan_db"] + ".1.bt2"
        method = "--bowtie2db"
    elif PARAMS["metaphlan_executable"] == "blast":
        assert os.path.exists(PARAMS["metaphlan_db"] + "nin"), "missing file %s: Are you sure you have the correct database for blast?" % PARAMS["metaphlan_db"] + "nin"
        method = "--blastdb"
    statement = PipelineMetagenomeAssembly.Metaphlan().build(infile, method="read_map")
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildMetaphlanReadmap, suffix(".readmap"), ".readmap.load")
def loadMetaphlanReadmaps(infile, outfile):
    '''
    load the metaphlan read maps
    '''
    P.load(infile,outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@merge(loadMetaphlanReadmaps, "metaphlan.dir/taxonomic.counts")
def countMetaphlanTaxonomicGroups(infiles, outfile):
    '''
    count the total number of species that
    were found by metaphlan
    '''
    outf = open(outfile, "w")
    outf.write("track\ttaxon_level\tcount\n")
    taxons = ["_order", "class", "family", "genus", "kingdom", "phylum", "species"]
    dbh = connect()
    cc = dbh.cursor()
    for infile in infiles:
        table = P.toTable(infile)
        track = P.snip(table, "_readmap")
        for taxon in taxons:
            count = cc.execute("""SELECT COUNT(DISTINCT %s) FROM %s""" % (taxon, table)).fetchone()[0]
            outf.write("\t".join([track, taxon, str(count)]) + "\n")
    outf.close()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@follows(mkdir("metaphlan.dir"))
@transform(preprocessReads
           , regex("(\S+).fa")
           , r"metaphlan.dir/\1.relab")
def buildMetaphlanRelativeAbundance(infile, outfile):
    '''
    metaphlan is a program used in metagenomics. It assigns
    reads to clades based on specific genetic markers via 
    blastn searching
    '''
    to_cluster = True
    # at present the pipeline will take a set of files
    # and compute the abundances of different taxonomic groups
    # based on ALL reads i.e. paired data are combined into
    # a single file for analysis
    if PARAMS["metaphlan_executable"] == "bowtie2":
        assert os.path.exists(PARAMS["metaphlan_db"] + ".1.bt2"), "missing file %s: Are you sure you have the correct database for bowtie2?" % PARAMS["metaphlan_db"] + ".1.bt2"
        method = "--bowtie2db"
    elif PARAMS["metaphlan_executable"] == "blast":
        assert os.path.exists(PARAMS["metaphlan_db"] + "nin"), "missing file %s: Are you sure you have the correct database for bowtie2?" % PARAMS["metaphlan_db"] + "nin"
        method = "--blastdb"

    statement = PipelineMetagenomeAssembly.Metaphlan().build(infile, method="rel_ab")
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildMetaphlanRelativeAbundance, suffix(".relab"), ".relab.load")
def loadMetaphlanRelativeAbundances(infile, outfile):
    '''
    load the metaphlan relative abundances
    '''
    P.load(infile,outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@merge(loadMetaphlanRelativeAbundances, "metaphlan.dir/taxonomic.abundances")
def buildMetaphlanTaxonomicAbundances(infiles, outfile):
    '''
    build a file that combines taxonomic abundances 
    from each sample
    '''
    dbh = connect()
    cc = dbh.cursor()
    outf = open(outfile, "w")
    outf.write("track\ttaxon_level\ttaxon\tabundance\tidx\n")
    for infile in infiles:
        table = P.toTable(infile)
        track = P.snip(table, "_relab")
        for data in cc.execute("""SELECT taxon_level, taxon, rel_abundance FROM %s""" % table).fetchall():
            idx = track.split("_")[1]
            outf.write("\t".join([track, data[0], data[1], str(data[2]), idx]) + "\n")
    outf.close()

#########################################
# taxonomic classification targets
#########################################
@follows(loadMetaphlanRelativeAbundances
         , buildMetaphlanTaxonomicAbundances
         , countMetaphlanTaxonomicGroups
         , loadMetaphlanReadmaps)
def metaphlan():
    pass

###################################################################
###################################################################
###################################################################
# Once we have an idea of the species that are in our samples we
# would like to align all of our reads to those genomes to assess
# the expected number of novel species i.e. reads that don't map
# to anything that is predicted to be present
###################################################################
###################################################################
###################################################################
@follows(mkdir("known_genomes.dir"))
@transform(buildMetaphlanRelativeAbundance
           , regex("(\S+)/(\S+).relab")
           , r"known_genomes.dir/\2.species")
def buildPresentSpeciesList(infile, outfile):
    '''
    build files with the list of species that are present
    in each sample
    '''
    to_cluster = True
    statement = '''cat %(infile)s | python %(scriptsdir)s/metaphlan2species.py --log=%(outfile)s.log > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################
@merge(buildPresentSpeciesList, "known_genomes.dir/species_present.tsv")
def buildUnionOfPresentSpecies(infiles, outfile):
    '''
    build a union set of present species for mapping reads to 
    genomes
    '''
    species_set = set()
    for infile in infiles:
        for line in open(infile).readlines():
            species = line[:-1]
            if species not in species_set:
                species_set.add(species)
    outf = open(outfile, "w")
    for species in species_set:
        outf.write("%s\n" % species)


###################################################################
###################################################################
###################################################################
@transform(buildUnionOfPresentSpecies
          , suffix(".tsv")
          , ".fa.gz")
def buildPresentSpeciesMultiFasta(infile, outfile):
    '''
    build a multi fasta file for alignment purposes
    '''
    genomes=PARAMS["known_species_genomesdir"]
    statement = '''python %(scriptsdir)s/species2multifasta.py
                   --level=genus
                   --genomes=%(genomes)s
                   --species=%(infile)s
                   --log=%(outfile)s.log
                   | gzip
                   > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################
@transform(buildPresentSpeciesMultiFasta, suffix(".fa.gz"), ".fa.load")
def loadGenomesAnalysed(infile, outfile):
    '''
    load the species that go into the mapping analysis
    '''
    to_cluster = False
    tmp = P.getTempFilename()
    statement = '''zcat %(infile)s | grep ">" | sed 's/>//g' | sed 's/|/ /g' > %(tmp)s'''
    P.run()

    tablename = P.toTable(outfile)
    statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s --log=%(outfile)s.log
                   < %(tmp)s > %(outfile)s; rm -rf %(tmp)s'''
    
    P.run()

###################################################################
###################################################################
###################################################################
index_suffix = {"bowtie":".ebwt"
          , "bowtie2":".bt2"
          , "bwa":".bwt"}

###################################################################
###################################################################
###################################################################
@transform(buildPresentSpeciesMultiFasta
           , suffix(".gz")
           , ".gz%s" % index_suffix[PARAMS["known_species_aligner"]])
def buildIndexOnPresentSpecies(infile, outfile):
    '''
    build index for the known species genomes
    '''
    if outfile.endswith(".bwt"):
        statement = '''bwa index %(infile)s >& %(outfile)s.log'''
    elif outfile.endswith(".bt2"):
        statement = '''bowtie build -f %(infile)s >& %(outfile)s.log'''
    elif outfile.endswith(".ebwt"):
        statement = '''bowtie build -f %(infile)s >& %(outfile)s.log'''
    P.run()

###################################################################
###################################################################
###################################################################
@follows(buildIndexOnPresentSpecies)
@transform(SEQUENCEFILES
           , SEQUENCEFILES_REGEX
           , add_inputs(buildPresentSpeciesMultiFasta)
           , r"known_genomes.dir/\1.bam")
def mapReadsAgainstKnownSpecies(infiles, outfile):
    '''
    map raw reads against known species
    '''
    genome = os.path.basename(infiles[1])    
    if index_suffix[PARAMS["known_species_aligner"]].endswith(".bwt"):
        infile = infiles[0]
        bwa_index_dir = os.path.dirname(outfile)
        bwa_aln_options = PARAMS["known_species_bwa_aln_options"]
        bwa_sampe_options=PARAMS["known_species_bwa_sampe_options"]
        bwa_threads=PARAMS["known_species_bwa_threads"]
        m = PipelineMapping.BWA()

    elif index_suffix[PARAMS["known_species_aligner"]].endswith(".ebwt"):
        bowtie_index_dir = os.path.dirname(outfile)
        bowtie_options=PARAMS["known_species_bowtie_options"]
        infile, reffile = infiles[0],  os.path.join(bowtie_index_dir, genome) + ".fa.gz"
        m = PipelineMapping.Bowtie( executable = P.substituteParameters( **locals() )["known_species_executable"] )

    elif index_suffix[PARAMS["known_species_aligner"]].endswith(".bt2"):
        bowtie2_index_dir = os.path.dirname(outfile)
        bowtie2_options=PARAMS["known_species_bowtie2_options"]
        infile, reffile = infiles[0],  os.path.join(bowtie2_index_dir, genome) + ".fa.gz"
        m = PipelineMapping.Bowtie2( executable = P.substituteParameters( **locals() )["known_species_executable"] )

    statement = m.build( (infile,), outfile ) 
    P.run()

###################################################################
###################################################################
###################################################################
@transform(mapReadsAgainstKnownSpecies
           , regex("(\S+).dir/(\S+).bam")
           , add_inputs(buildPresentSpeciesMultiFasta)
           , r"\1-\1.dir/\2.picard_stats")
def buildPicardStatsOnKnownSpeciesAlignments(infiles, outfile):
    '''
    build statistics for the alignment of reads against known
    species
    '''
    reffile = infiles[1]
    infile = infiles[0]
    PipelineMappingQC.buildPicardAlignmentStats( infile, 
                                                 outfile,
                                                 reffile )


###################################################################
###################################################################
###################################################################
@jobs_limit(1, "db")
@merge( buildPicardStatsOnKnownSpeciesAlignments, "known_genomes.dir/picard_stats.load" )
def loadPicardStatsOnKnownSpeciesAlignments( infiles, outfile ):
    '''merge alignment stats into single tables.'''

    PipelineMappingQC.loadPicardAlignmentStats( infiles, outfile )

###################################################################
###################################################################
###################################################################
@follows(loadPicardStatsOnKnownSpeciesAlignments)
def presentSpeciesAlignment():
    pass

###################################################################
###################################################################
###################################################################
# Functional assignments using rpsblast and COG categories
###################################################################
###################################################################
###################################################################
@follows(mkdir("function.dir"))
@transform(SEQUENCEFILES
           , SEQUENCEFILES_REGEX
           , r"function.dir/\1.blast.gz")
def runBlastOnRawSequences(infile, outfile):
    '''
    run blast on raw reads for downstream analysis using
    MEGAN assignment to KEGG categories. Outputs blastx
    formatted data after a translated blast against the
    specified database. The chunksize is fixed at 1000 - 
    this may become parametised in the .ini file in the
    future. 
    TODO: Parameterise more blast options
    '''
    db = PARAMS["megan_db"]
    evalue = PARAMS["megan_evalue"]
    temp = P.getTempFilename(".")
    statement = '''fastqToFa %(infile)s %(temp)s; checkpoint 
                  ; cat %(temp)s | python %(scriptsdir)s/farm.py --split-at-regex="^>(\S+)" 
                    --log=%(outfile)s.log 
                    --chunksize=5 "blastx -db %(db)s -evalue %(evalue)s" 
                  | gzip > %(outfile)s; checkpoint
                  ; rm -rf %(temp)s'''
    P.run()

###################################################################
###################################################################
###################################################################
@transform(runBlastOnRawSequences, suffix(".blast.gz"), ".kegg.gz")
def assignKeggFunctions(infile, outfile):
    '''
    assign kegg functions to the aligned reads using
    the lcamapper script that forms part of MEGAN
    '''
    track = P.snip(outfile, ".kegg.gz")
    outf_tax = P.snip(outfile.replace("kegg", "taxonomy"), ".gz")
    statement = '''lcamapper.sh 
                   -k
                   -i %(infile)s
                   -o %(outf_tax)s > %(outfile)s.log; checkpoint
                   ; gzip %(track)s.blast-kegg.txt
                   ; gzip %(outf_tax)s
                   ; mv %(track)s.blast-kegg.txt.gz %(outfile)s
                '''
    P.run()


###################################################
# DEPRECATED but may become useful in the future so
# has been kept in for the time being
###################################################
#     '''
#     run a translated blast on raw sequenced reads
#     '''
#     to_cluster = True    
#     job_options = job_options = " -l mem_free=30G"

#     db = PARAMS["rpsblast_db"]
#     evalue = PARAMS["rpsblast_evalue"]
    
#     # at the moment this only considers
#     # one read in a pair - needs to be adapted for
#     # paired data

#     # converts to fasta on the fly using sed
#     statement = '''zcat %(infile)s
#                   | sed -n '1~4s/^@/>/p;2~4p'
#                   | python %(scriptsdir)s/farm.py 
#                   --split-at-regex="^>(\S+)" 
#                   --chunksize=100000
#                   "rpsblast -db %(db)s 
#                   -evalue %(evalue)s
#                   -soft_masking True
#                   -outfmt 6"
#                   | gzip > %(outfile)s'''
#     P.run()


# ###################################################################                                                                                                                                                                          
# ###################################################################                                                                                                                                                                          
# ###################################################################                                                                                                                                                                          
# @transform(runBlastOnRawSequences
#            , suffix(".gz")
#            , ".cog.gz")
# def assignCOGsToAlignments(infile, outfile):
#     '''
#     assign clusters of orthologous groups ids to the
#     protein alignment
#     '''
#     job_options = " -l mem_free=30G"
#     mapfile = PARAMS["COG_map"]
    
#     statement = '''zcat %(infile)s 
#                   | python %(scriptsdir)s/rpsblast_cdd2cog.py 
#                   --cog-map=%(mapfile)s --log=%(outfile)s.log
#                   | gzip > %(outfile)s'''
#     P.run()

# ###################################################################                                                                                                                                                                          
# ###################################################################                                                                                                                                                                          
# ###################################################################                                                                                                                                                                          
# @transform(assignCOGsToAlignments
#            , suffix(".gz")
#            , add_inputs(countReads)
#            , ".counts.gz")
# def countCOGAssignments(infiles, outfile):
#     '''
#     calculate the proportion of reads that are assigned to each
#     COG and add annotation
#     '''
#     to_cluster = True
#     job_options = " -l mem_free=30G"
#     countsfile = [x for x in infiles[1:len(infiles)] if infiles[0].find(P.snip(x, ".nreads")) != -1]
#     countsfile = countsfile[0]
    
#     # TODO: sort this out for single ended data
#     total = open(countsfile).readline().split("\t")[1][:-1]
# #    total = int(float(total)/2)
#     total = int(total)

#     description_file = PARAMS["COG_description"]
#     inf = infiles[0]

#     statement = '''zcat %(inf)s 
#                    | python %(scriptsdir)s/rpsblast_cog2counts.py
#                    --cog-description=%(description_file)s
#                    --nreads=%(total)s
#                    --log=%(outfile)s.log
#                    | gzip > %(outfile)s'''

#     P.run()

# ###################################################################                                                                                                                                                                         
# ###################################################################                                                                                                                                                                          
# ###################################################################                                                                                                                                                                          
# @transform(countCOGAssignments, suffix(".gz"), ".load")
# def loadCOGCounts(infile, outfile):
#     '''
#     load COG funtion counts
#     '''
#     to_cluster = False
#     # need to preprocess the data - this is dependent
#     # on the map file that was used in the previous step
#     temp = P.getTempFile()
#     function2counts = collections.defaultdict(float)
#     inf = IOTools.openFile(infile)
#     header = inf.readline()
#     for line in inf.readlines():
#         data = line[:-1].split("\t")
#         function, prop = data[1], data[2]
#         function2counts[function] += float(prop)

#     temp.write("funtion\tproportion\n")
#     for function, prop in function2counts.iteritems():
#         temp.write("%s\t%f\n" % (function, prop))
#     temp.close()

#     temp_name = temp.name

#     tablename = P.toTable(outfile)
#     statement = '''python %(scriptsdir)s/csv2db.py
#                    -t %(tablename)s
#                    --log=%(outfile)s.log
#                    < %(temp_name)s > %(outfile)s
#                    '''
#     P.run()


# ###################################################################                                                                                                                                                                         
# ###################################################################                                                                                                                                                                          
# ###################################################################                                                                                                                                                                          
# @transform(loadCOGCounts
#            , suffix(".load")
#            , ".stats")
# def buildRpsblastAlignmentStats(infile, outfile):
#     '''
#     count the proportion of reads that can be aligned 
#     to protein CDD sequences - this is done by subtraction
#     from COG counts = Unassigned COG + Unaligned read
#     '''
#     dbh = connect()
#     cc = dbh.cursor()

#     tablename = P.toTable(infile)
#     total = 0 
#     for data in cc.execute("""SELECT proportion FROM %s""" % tablename).fetchall():
#         prop = data[0]
#         total += prop
    
#     unmapped_in_some_way = 1 - total
#     outf = open(outfile, "w")
#     track = P.snip(infile, ".rpblast.result.cog.counts.load")
#     outf.write("%s\t%f\n" % (track, unmapped_in_some_way))

# ###################################################################                                                                                                                                                                          
# ###################################################################                                                                                                                                                                          
# ###################################################################                                                                                                                                                                          
# @merge(buildRpsblastAlignmentStats, "function.dir/rpsblast_alignment_stats.load")
# def loadRpsblastAlignmentStats(infiles, outfile):
#     '''
#     load alignment statistics
#     '''
#     to_cluster = False
#     tmp = P.getTempFilename()
#     infs = " ".join(infiles)
#     statement = '''cat %(infs)s > %(tmp)s'''
#     P.run()
#     tablename = P.toTable(outfile)
#     statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s --header=track,proportion
#                    --log=%(outfile)s.log < %(tmp)s > %(outfile)s'''
#     P.run()

##################################################################
@follows(runBlastOnRawSequences)
def functional_profile():
    pass

###################################################################                                                                                                                                                                         
# Have reads been pooled
###################################################################                                                                                                                                                                          
SEQUENCE_TARGETS = {1: (poolReadsAcrossConditions, regex("(\S+)/(\S+).(fasta$|fasta.gz|fasta.1.gz|fastq$|fastq.gz|fastq.1.gz)"), "2.contigs.fa")
                    , 0: (SEQUENCEFILES, SEQUENCEFILES_REGEX, "1.contigs.fa")
                    , "": (SEQUENCEFILES, SEQUENCEFILES_REGEX, "1.contigs.fa")}

###################################################################                                                                                                                                                                         
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
## assemble reads with meta-velvet
###################################################################                                                                                                                                                                          
@active_if("metavelvet" in ASSEMBLERS)
@follows(mkdir("metavelvet.dir"))
@transform( SEQUENCE_TARGETS[PARAMS["pool_reads"]][0]
            , SEQUENCE_TARGETS[PARAMS["pool_reads"]][1]
            , r"metavelvet.dir/\%s" % SEQUENCE_TARGETS[PARAMS["pool_reads"]][2])
def runMetavelvet(infile, outfile):
    '''
    run meta-velvet on each track
    '''
    job_options = " -l mem_free=30G"
    statement = PipelineMetagenomeAssembly.Metavelvet().build(infile)
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@jobs_limit(1, "R")
@transform(runMetavelvet
           , suffix(".contigs.fa")
           , ".stats.pdf")
def plotCoverageHistogram(infile, outfile):
    '''
    plot the coverage over kmers
    '''
    inf = P.snip(infile, ".contigs.fa") + ".stats.txt"
    outf = P.snip(inf, ".txt") + ".pdf"
    R('''library(plotrix)''')
    R('''data = read.table("%s", header=TRUE)''' % inf)
    R('''pdf("%s", height = 7, width = 7 )''' % outf)
    R('''weighted.hist(data$short1_cov, data$lgth, breaks=seq(0, 200, by=1))''')
    R["dev.off"]()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(runMetavelvet
           , suffix(".contigs.fa")
           , ".stats.load")
def loadMetavelvetRawStats(infile, outfile):
    '''
    load the assembly stats for meta-velvet
    '''
    inf = P.snip(infile, ".contigs.fa") + ".stats.txt"
    P.load(inf, outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(runMetavelvet, suffix(".contigs.fa"), ".summary.tsv")
def buildMetavelvetStats(infile, outfile):
    '''
    build metavelvet stats:
    N50
    Number of scaffolds 
    Total scaffold length
    '''
    PipelineMetagenomeAssembly.contig_to_stats(infile, outfile, PARAMS)
    
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildMetavelvetStats, regex("(\S+).dir/(\S+).tsv"), r"\1.dir/\1-\2.load")
def loadMetavelvetStats(infile, outfile):
    '''
    load the metavelvet stats
    '''
    P.load(infile, outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
## assemble reads with idba
###################################################################                                                                                                                                                                          
IDBA_TARGETS = {1: (preprocessIdba, regex("(\S+)/(\S+).fa"), "2.contigs.fa")
                , 0: (preprocessReads, regex("(\S+).fa"), "1.contigs.fa")
                , "": (preprocessReads, regex("(\S+).fa"), "1.contigs.fa")}

@active_if("idba" in ASSEMBLERS)
@follows(mkdir("idba.dir"))
@transform(IDBA_TARGETS[PARAMS["pool_reads"]][0]
           , IDBA_TARGETS[PARAMS["pool_reads"]][1]
           , r"idba.dir/\%s" % IDBA_TARGETS[PARAMS["pool_reads"]][2])
def runIdba(infile, outfile):
    '''
    run idba on each track
    '''
    to_cluster = True
    job_options = " -l mem_free=30G"
    statement = PipelineMetagenomeAssembly.Idba().build(infile)
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(runIdba, suffix(".contigs.fa"), ".summary.tsv")
def buildIdbaStats(infile, outfile):
    '''
    build idba stats:
    N50
    Number of scaffolds 
    Total scaffold length
    '''
    PipelineMetagenomeAssembly.contig_to_stats(infile, outfile, PARAMS)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildIdbaStats, regex("(\S+).dir/(\S+).tsv"), r"\1.dir/\1-\2.load")
def loadIdbaStats(infile, outfile):
    '''
    load the idba stats
    '''
    P.load(infile, outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if("ray" in ASSEMBLERS)
@follows(mkdir("ray.dir"))
@transform( SEQUENCE_TARGETS[PARAMS["pool_reads"]][0]
            , SEQUENCE_TARGETS[PARAMS["pool_reads"]][1]
            , r"ray.dir/\%s" % SEQUENCE_TARGETS[PARAMS["pool_reads"]][2])
def runRay(infile, outfile):
    '''
    run Ray on each track
    '''
    to_cluster = True
    job_options=" -pe mpi 1 -q mpi.q -l mem_free=30G "
    statement = PipelineMetagenomeAssembly.Ray().build(infile)
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if("sga" in ASSEMBLERS)
@follows(mkdir("sga.dir"))
@transform( SEQUENCE_TARGETS[PARAMS["pool_reads"]][0]
            , SEQUENCE_TARGETS[PARAMS["pool_reads"]][1]
            , r"sga.dir/\%s" % SEQUENCE_TARGETS[PARAMS["pool_reads"]][2])
def runSGA(infile, outfile):
    '''
    run SGA on each track
    '''
    to_cluster = True
    job_options = " -l mem_free=30G "
    statement = PipelineMetagenomeAssembly.SGA().build(infile) 
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if("soapdenovo" in ASSEMBLERS)
@follows(mkdir("soapdenovo.dir"))
@transform( SEQUENCE_TARGETS[PARAMS["pool_reads"]][0]
            , SEQUENCE_TARGETS[PARAMS["pool_reads"]][1]
            , r"soapdenovo.dir/\%s.cfg" % SEQUENCE_TARGETS[PARAMS["pool_reads"]][2])
def buildSoapdenovoConfig(infile, outfile):
    '''
    run SGA on each track
    '''
    PipelineMetagenomeAssembly.SoapDenovo2().config(infile, outfile, PARAMS)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildSoapdenovoConfig, suffix(".contigs.fa.cfg"), ".contigs.fa")
def runSoapdenovo(infile, outfile):
    '''
    run soapdenovo
    '''
    job_options="-l mem_free=30G"
    statement = PipelineMetagenomeAssembly.SoapDenovo2().build(infile)
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
ASSEMBLY_TARGETS = []
assembly_targets = {"metavelvet": runMetavelvet
                    , "idba": runIdba
                    , "ray": runRay
                    , "sga": runSGA
                    , "soapdenovo": runSoapdenovo
                    }
for x in ASSEMBLERS:
    ASSEMBLY_TARGETS.append(assembly_targets[x])

print ASSEMBLY_TARGETS
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(ASSEMBLY_TARGETS, suffix(".contigs.fa"), ".filtered.contigs.fa")
def filterContigs(infile, outfile):
    '''
    filter contigs if specified in .ini file. If not specified
    then the pipeline will not remove any but will produce a new 
    outfile - this is not space efficient and SHOULD BE CHANGED
    '''
    if not PARAMS["filter"]:
        length = 0
    else:
        length = PARAMS["filter"]
    
    PipelineMetagenomeAssembly.filterContigs(infile, outfile, length)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(filterContigs, suffix(".fa"), ".summary.tsv")
def buildContigStats(infile, outfile):
    '''
    build contig stats:
    N50
    Number of scaffolds 
    Total scaffold length
    max length
    '''
    PipelineMetagenomeAssembly.contig_to_stats(infile, outfile, PARAMS)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildContigStats, regex("(\S+).dir/(\S+).tsv"), r"\1.dir/\1-\2.load")
def loadContigStats(infile, outfile):
    '''
    load the contig stats
    '''
    P.load(infile, outfile)

@follows(loadContigStats)
def contig_stats():
    pass

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@split(loadContigStats, "*/contig.summary.tsv")
def buildContigSummary(infiles, outfile):
    '''
    merge the contig summary statistics
    '''
    stats = collections.defaultdict(list)
    for filepath in infiles:
        dirname = os.path.dirname(filepath)
        stats[dirname].append(os.path.basename(filepath))

    N = PARAMS["scaffold_n"]

    # connect to database
    dbh = connect()
    cc = dbh.cursor()
    for dirname in stats.keys():
        outfname = os.path.join(dirname, "contig.summary.tsv")
        outf = open(outfname, "w")
        outf.write("track\tnscaffolds\tscaffold_length\tN%i\tmean_length\tmedian_length\tmax_length\n" % N)
        for infile in stats[dirname]:
            track = P.snip(infile.split(dirname.split(".dir")[0])[1][1:], ".summary.load")
            table = P.toTable(infile)
            data = cc.execute("""SELECT nscaffolds
                                 , scaffold_length
                                 , N50
                                 , mean_length
                                 , median_length
                                 , max_length FROM %s""" % table).fetchone()
            outf.write("\t".join(map(str, [track, data[0], data[1], data[2], data[3], data[4], data[5]])) + "\n")
        outf.close()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildContigSummary, suffix(".tsv"), ".load")
def loadContigSummary(infile, outfile):
    '''
    load contig summary stats for each assembler
    '''
    outname = P.snip(os.path.dirname(infile), ".dir") + "_" + os.path.basename(infile) + ".load"
    P.load(infile, outname)
    P.touch(outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(filterContigs, suffix(".fa"), ".lengths.tsv")
def buildContigLengths(infile, outfile):
    '''
    output lengths for each contig in each of the assemblies
    '''
    PipelineMetagenomeAssembly.build_scaffold_lengths(infile, outfile, PARAMS)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildContigLengths, suffix(".lengths.tsv"), ".lengths.load")
def loadContigLengths(infile, outfile):
    '''
    load contig lengths
    '''
    outname = P.snip(os.path.dirname(infile), ".dir") + "_" + P.snip(os.path.basename(infile), ".tsv") + ".load"
    P.load(infile, outname, "--index=scaffold_name")
    P.touch(outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(filterContigs, suffix(".fa"), ".gc.tsv")
def buildContigGCContent(infile, outfile):
    '''
    build the GC content for each contig
    '''
    statement = '''cat %(infile)s 
                   | python %(scriptsdir)s/fasta2table.py 
                   --section=cpg
                   --log=%(outfile)s.log
                   > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildContigGCContent, suffix(".gc.tsv"), ".gc.load")
def loadContigGCContent(infile, outfile):
    '''
    load contig GC content
    '''
    outname = P.snip(os.path.dirname(infile), ".dir") + "_" + P.snip(os.path.basename(infile), ".tsv") + ".load"
    P.load(infile, outname, "--index=id")
    P.touch(outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(filterContigs, suffix(".fa"), ".blast")
def runBlastOnContigs(infile, outfile):
    '''
    run blast on the contigs for downstream taxonomic assignment
    runs a translated blast (x) and outputs blastx format
    or imput into MEGAN
    '''
    db = PARAMS["megan_db"]
    evalue = PARAMS["megan_evalue"]
    statement = '''cat %(infile)s 
                  | python %(scriptsdir)s/farm.py --split-at-regex="^>(\S+)" 
                    --log=%(outfile)s.log 
                    --chunksize=100 "blastx -db %(db)s -evalue %(evalue)s" > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(runBlastOnContigs, suffix(".blast"), ".lca")
def runLCA(infile, outfile):
    '''
    run the lowest common ancestor algorithm
    on the blast output to assign contigs to 
    taxa - from mtools. Runs with defaults at
    the moment
    '''
    statement = '''lcamapper.sh 
                   -i %(infile)s
                   -o %(outfile)s > %(outfile)s.log'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(runLCA, suffix(".lca"), ".taxa.gz")
def parseLCA(infile, outfile):
    '''
    tabulate LCA output into nice format
    '''
    statement = '''cat %(infile)s 
                   | python %(scriptsdir)s/lca2table.py
                     --summarise=individual
                     --log=%(outfile)s.log
                   | sed -e 's/order/_order/g'
                   | gzip > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@jobs_limit(1, "db")
@transform(parseLCA, suffix(".gz"), ".load")
def loadLCA(infile, outfile):
    '''
    load LCA results
    '''
    tablename = P.snip(os.path.dirname(infile), ".dir")+"_"+os.path.basename(P.snip(infile, ".gz"))
    tablename = P.toTable(tablename + ".load")
    statement = '''zcat %(infile)s | python %(scriptsdir)s/csv2db.py
                  -t %(tablename)s
                  --index=id
                  --log=%(outfile)s.log
                  > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(filterContigs, suffix(".fa"), ".tetra")
def buildTetranucleotideFreq(infile, outfile):
    '''
    calculate the tetranucleotide frequency for
    contigs
    '''
    statement = '''cat %(infile)s | python %(scriptsdir)s/fasta2kmercontent.py
                   -k 4 --log=%(outfile)s.log > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildTetranucleotideFreq, suffix(".tetra"), ".tetra.load")
def loadTetranucleotideFreq(infile, outfile):
    '''
    load tetranucleotide frequency matrix
    '''
    P.load(infile, outfile, "--index=contig")

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
## gene finding using MetaGeneMark
###################################################################                                                                                                                                                                          
@transform(filterContigs, suffix(".fa"), ".genes.tsv")
def findGenesUsingMetaGeneMark(infile, outfile):
    '''
    Use the MetaGeneMark hmm software to predict genes.
    Output is tsv - similar to gff format but with
    sequences
    '''
    to_cluster=True
    mparams = PARAMS["metagenemark_model_params"]
    statement = '''gmhmmp -a -d -f G -m %(mparams)s -o %(outfile)s %(infile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(findGenesUsingMetaGeneMark, regex("(\S+).tsv"), r"\1.gff.gz")
def parseGenesGff(infile, outfile):
    '''
    parse the genes file
    '''
    to_cluster = True
    statement = '''cat %(infile)s | python %(scriptsdir)s/formatMetagenemark.py 
                                    --format gff 
                                    --log=%(outfile)s.log 
                                     | gzip > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(findGenesUsingMetaGeneMark, regex("(\S+).tsv"), r"\1.fasta.gz")
def parseGenesFasta(infile, outfile):
    '''
    parse the genes file
    '''
    to_cluster = True
    statement = '''cat %(infile)s | python %(scriptsdir)s/formatMetagenemark.py 
                                    --format fasta 
                                    --log=%(outfile)s.log 
                                    | sed 's/DNA /DNA_/g' 
                                    | gzip > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(findGenesUsingMetaGeneMark, regex("(\S+).tsv"), r"\1.aa.gz")
def parseGenesAa(infile, outfile):
    '''
    parse the genes file
    '''
    to_cluster = True
    statement = '''cat %(infile)s | python %(scriptsdir)s/formatMetagenemark.py 
                                    --format aa 
                                    --log=%(outfile)s.log 
                                    | sed 's/Protein /Protein_/g' 
                                    | gzip > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(parseGenesAa, suffix(".aa.gz"), ".essential.hmm.gz")
def assignEssentialGenesToContigs(infile, outfile):
    '''
    assign essential genes to contigs
    '''
    dirname = os.path.dirname(infile)
    essential = PARAMS["hmmer_hmm"]
    tempdir = P.getTempDir(".")

    statement = '''zcat %(infile)s > %(tempdir)s/orfs.fa;
                   hmmsearch --tblout %(tempdir)s/hmm.out --cut_tc
                   --notextw  %(essential)s %(tempdir)s/orfs.fa;
                    tail -n+4 %(tempdir)s/hmm.out | sed 's/ * / /g' | cut -f1,4 -d " " 
                   | gzip > %(outfile)s'''
    P.run()
    statement = '''rm -rf %(tempdir)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(assignEssentialGenesToContigs, suffix(".gz")
           , add_inputs(parseGenesGff)
           , ".contigs.gz")
def postprocessEssentialGeneAssignments(infiles, outfile):
    '''
    need to add the contig that each orf is associates with to the 
    file
    '''
    track = P.snip(os.path.basename(infiles[0]), ".essential.hmm.gz")
    genes, gff = infiles[0], [inf for inf in infiles[1:] if inf.find(track) != -1][0]
    protein2contig = {}
    for gff in GTF.iterator(IOTools.openFile(gff)):
        protein2contig["Protein_" + str(gff.gene_id)] = gff.contig

    # output contig associated with protein id
    outf = IOTools.openFile(outfile, "w")
    outf.write("contig\torf\thmm_profile\n")
    for line in IOTools.openFile(genes).readlines():
        data = line[:-1].split(" ")
        protein, profile = data[0], data[1]
        outf.write("\t".join([protein2contig[protein], protein, profile]) + "\n")
    outf.close()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(postprocessEssentialGeneAssignments, suffix(".gz"), ".load")
def loadEssentialGeneAssignments(infile, outfile):
    '''
    load assignments of essential genes 
    '''
    P.load(infile, outfile, "--index=contig")

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@follows(mkdir("genes.dir"))
@transform(parseGenesAa, regex("(\S+).dir/(\S+).aa.gz"), r"genes.dir/\1_\2.blast.result.gz")
def runBlastOnAminoAcidSequences(infile, outfile):
    '''
    look for homology with known genes
    '''
    to_cluster = True
    db = PARAMS["blastp_db"]
    evalue = PARAMS["blastp_evalue"]
    if PARAMS["blastp_ungapped"]:
        ungapped = "-ungapped"
    else:
        ungapped = ""

    statement = '''zcat %(infile)s 
                  | python %(scriptsdir)s/farm.py --split-at-regex="^>(\S+)" 
                    --log=%(outfile)s.log 
                    --chunksize=1000 "blastp -db %(db)s -evalue %(evalue)s
                    -outfmt '6 qseqid qlen sseqid sgi sacc slen qstart qend sstart send evalue bitscore score length pident mismatch gaps frames staxids sscinames'"
                    | gzip > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(runBlastOnAminoAcidSequences, suffix(".result"), r".result.load")  
def loadBlastOnAminoAcidSequences(infile, outfile):
    '''
    load blastp results
    '''
    P.load(infile, outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform([parseGenesGff,parseGenesFasta,parseGenesAa]
           ,regex("(\S+).dir/(\S+).genes.(\S+).gz")
           , r"\1.dir/\1_\2.genes.\3.tsv.gz")
def buildGeneTables(infile, outfile):
    '''
    build gene tables
    '''
    to_cluster = True
    if infile.endswith(".gff.gz"): 
        outf = gzip.open(outfile, "w")
        outf.write("chr\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattributes\n")
        for line in gzip.open(infile).readlines():
            outf.write(line)
        outf.close()
    else:
        statement = '''zcat %(infile)s | python %(scriptsdir)s/fasta2table.py 
                                         -s sequence 
                                         --log=%(outfile)s.log | gzip > %(outfile)s'''
        P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
jobs_limit(1, "db")
@transform(buildGeneTables, regex("(\S+)/(\S+).genes.(\S+).tsv.gz"), r"\1/\2.genes.\3.load")
def loadGeneTables(infile, outfile):
     '''
     load genes from metagenemaek analysis
     '''
     if infile.find("gff") != -1:
         P.load(infile, outfile)
     else:
         P.load(infile, outfile)





###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@follows(loadGeneTables)
def metagenemark():
    pass

@follows(loadBlastOnAminoAcidSequences)
def geneSimilarity():
    pass

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
## build indices for mapping - this is for coverage analysis
###################################################################                                                                                                                                                                          
@active_if(BOWTIE)
@transform(filterContigs, suffix(".fa"), ".ebwt")
def buildAssemblyBowtieIndices(infile, outfile):
    '''
    build bowtie indices
    '''
    to_cluster = True
    outbase = P.snip(infile, ".fa")
    directory = os.path.dirname(infile)
    statement = '''bowtie-build -f %(infile)s %(outbase)s'''
    P.run()
    P.touch(outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if(BOWTIE2)
@transform(filterContigs, suffix(".fa"), ".bt2")
def buildAssemblyBowtie2Indices(infile, outfile):
    '''
    build bowtie indices
    '''
    to_cluster = True
    outbase = P.snip(infile, ".fa")
    statement = '''bowtie2-build -f %(infile)s %(outbase)s'''
    P.run()
    P.touch(outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if(BWA)
@transform(filterContigs, suffix(".fa"), ".fa.bwt")
def buildAssemblyBWAIndices(infile, outfile):
    '''
    build bwa indices
    '''
    to_cluster = True
    statement = '''bwa index %(infile)s'''
    P.run()
    P.touch(outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
## map 
###################################################################                                                                                                                                                                          
index = {"bowtie":buildAssemblyBowtieIndices
         , "bowtie2":buildAssemblyBowtie2Indices
         , "bwa":buildAssemblyBWAIndices}
INDEX = index[MAPPER]

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if("metavelvet" in ASSEMBLERS)
@transform(SEQUENCEFILES
           , SEQUENCEFILES_REGEX
           , add_inputs(INDEX, runMetavelvet)
           , r"metavelvet.dir/\1.filtered.contigs.bam")
def mapReadsAgainstMetavelvetContigs(infiles, outfile):
    '''
    map reads against metavelvet contigs
    '''
    inf = infiles[0]
    to_cluster = True
    index_dir = os.path.dirname(outfile)

    if "agg" not in infiles[1]:
        genome = re.search(".*R[0-9]*", infiles[0]).group(0) + ".filtered.contigs.fa"
    else:
        genome = "agg-agg-agg.filtered.contigs"

    if infiles[1].endswith(".bt2") or infiles[1].endswith(".ebwt"):
        infile, reffile = infiles[0],  os.path.join(index_dir, genome) + ".fa"
        m = PipelineMapping.Bowtie( executable = P.substituteParameters( **locals() )["bowtie_executable"] )

    elif infiles[1].endswith("bwt"):
        genome = genome + ".fa"
        job_options= " -l mem_free=%s" % (PARAMS["bwa_memory"])
        bwa_index_dir = index_dir
        bwa_aln_options = PARAMS["bwa_aln_options"]
        bwa_sampe_options=PARAMS["bwa_sampe_options"]
        bwa_threads=PARAMS["bwa_threads"]
        m = PipelineMapping.BWA(remove_non_unique = True)
    statement = m.build( (inf,), outfile ) 
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if("idba" in ASSEMBLERS)
@transform(SEQUENCEFILES
           , SEQUENCEFILES_REGEX
           , add_inputs(INDEX, runIdba)
           , r"idba.dir/\1.filtered.contigs.bam")
def mapReadsAgainstIdbaContigs(infiles, outfile):
    '''
    map reads against idba contigs
    '''
    inf = infiles[0]
    to_cluster = True
    index_dir = os.path.dirname(outfile)

    if "agg" not in infiles[1]:
        genome = re.search(".*R[0-9]*", infiles[0]).group(0) + ".filtered.contigs.fa"
    else:
        genome = "agg-agg-agg.filtered.contigs.fa"

    if infiles[1].endswith(".bt2") or infiles[1].endswith(".ebwt"):
        infile, reffile = infiles[0],  os.path.join(index_dir, genome)
        m = PipelineMapping.Bowtie( executable = P.substituteParameters( **locals() )["bowtie_executable"] )

    elif infiles[1].endswith("bwt"):
        genome = genome
        job_options= " -l mem_free=%s" % (PARAMS["bwa_memory"])
        bwa_index_dir = index_dir
        bwa_aln_options = PARAMS["bwa_aln_options"]
        bwa_sampe_options=PARAMS["bwa_sampe_options"]
        bwa_threads=PARAMS["bwa_threads"]
        m = PipelineMapping.BWA(remove_non_unique = True)
    statement = m.build( (inf,), outfile ) 
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if("ray" in ASSEMBLERS)
@transform(SEQUENCEFILES
           , SEQUENCEFILES_REGEX
           , add_inputs(INDEX, runRay)
           , r"ray.dir/\1.filtered.contigs.bam")
def mapReadsAgainstRayContigs(infiles, outfile):
    '''
    map reads against Ray contigs
    '''
    inf = infiles[0]
    to_cluster = True
    index_dir = os.path.dirname(outfile)

    if "agg" not in infiles[1]:
        genome = re.search(".*R[0-9]*", infiles[0]).group(0) + ".filtered.contigs.fa"
    else:
        genome = "agg-agg-agg.filtered.contigs.fa"

    if infiles[1].endswith(".bt2") or infiles[1].endswith(".ebwt"):
        infile, reffile = infiles[0],  os.path.join(index_dir, genome) + ".fa"
        m = PipelineMapping.Bowtie( executable = P.substituteParameters( **locals() )["bowtie_executable"] )

    elif infiles[1].endswith("bwt"):
        genome = genome
        job_options= " -l mem_free=%s" % (PARAMS["bwa_memory"])
        bwa_index_dir = index_dir
        bwa_aln_options = PARAMS["bwa_aln_options"]
        bwa_sampe_options=PARAMS["bwa_sampe_options"]
        bwa_threads=PARAMS["bwa_threads"]
        m = PipelineMapping.BWA(remove_non_unique = True)
    statement = m.build( (inf,), outfile ) 
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if("sga" in ASSEMBLERS)
@transform(SEQUENCEFILES
           , SEQUENCEFILES_REGEX
           , add_inputs(INDEX, runSGA)
           , r"sga.dir/\1.filtered.contigs.bam")
def mapReadsAgainstSGAContigs(infiles, outfile):
    '''
    map reads against Ray contigs
    '''
    inf = infiles[0]
    to_cluster = True
    index_dir = os.path.dirname(outfile)

    if "agg" not in infiles[1]:
        genome = re.search(".*R[0-9]*", infiles[0]).group(0) + ".filtered.contigs.fa"
    else:
        genome = "agg-agg-agg.filtered.contigs.fa"

    if infiles[1].endswith(".bt2") or infiles[1].endswith(".ebwt"):
        infile, reffile = infiles[0],  os.path.join(index_dir, genome) + ".fa"
        m = PipelineMapping.Bowtie( executable = P.substituteParameters( **locals() )["bowtie_executable"] )

    elif infiles[1].endswith("bwt"):
        genome = genome
        job_options= " -l mem_free=%s" % (PARAMS["bwa_memory"])
        bwa_index_dir = index_dir
        bwa_aln_options = PARAMS["bwa_aln_options"]
        bwa_sampe_options=PARAMS["bwa_sampe_options"]
        bwa_threads=PARAMS["bwa_threads"]
        m = PipelineMapping.BWA(remove_non_unique = True)
    statement = m.build( (inf,), outfile ) 
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if("soapdenovo" in ASSEMBLERS)
@transform(SEQUENCEFILES
           , SEQUENCEFILES_REGEX
           , add_inputs(INDEX, runSoapdenovo)
           , r"soapdenovo.dir/\1.filtered.contigs.bam")
def mapReadsAgainstSoapdenovoContigs(infiles, outfile):
    '''
    map reads against Ray contigs
    '''
    inf = infiles[0]
    to_cluster = True
    index_dir = os.path.dirname(outfile)

    if "agg" not in infiles[1]:
        genome = re.search(".*R[0-9]*", infiles[0]).group(0) + ".filtered.contigs.fa"
    else:
        genome = "agg-agg-agg.filtered.contigs.fa"

    if infiles[1].endswith(".bt2") or infiles[1].endswith(".ebwt"):
        infile, reffile = infiles[0],  os.path.join(index_dir, genome) + ".fa"
        m = PipelineMapping.Bowtie( executable = P.substituteParameters( **locals() )["bowtie_executable"] )

    elif infiles[1].endswith("bwt"):
        genome = genome
        job_options= " -l mem_free=%s" % (PARAMS["bwa_memory"])
        bwa_index_dir = index_dir
        bwa_aln_options = PARAMS["bwa_aln_options"]
        bwa_sampe_options=PARAMS["bwa_sampe_options"]
        bwa_threads=PARAMS["bwa_threads"]
        m = PipelineMapping.BWA(remove_non_unique = True)
    statement = m.build( (inf,), outfile ) 
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
ALIGNMENT_TARGETS = []
alignment_targets = {"metavelvet":mapReadsAgainstMetavelvetContigs
                     , "idba":mapReadsAgainstIdbaContigs
                     , "ray":mapReadsAgainstRayContigs
                     , "sga":mapReadsAgainstSGAContigs
                     , "soapdenovo":mapReadsAgainstSoapdenovoContigs}
for x in ASSEMBLERS:
    ALIGNMENT_TARGETS.append(alignment_targets[x])

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(ALIGNMENT_TARGETS
           , regex("(\S+).dir/(\S+).bam")
           , r"\1.dir/\1_\2.alignment_stats")
def buildAlignmentStats(infile, outfile):
    '''
    use bam2stats to get alignment statistics
    '''
    statement = '''cat %(infile)s | python %(scriptsdir)s/bam2stats.py 
                   --log=%(outfile)s.log - 
                   > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildAlignmentStats, suffix("_stats"), "_stats.load")
def loadAlignmentStats(infile, outfile):
    '''
    load bam2stats results
    '''
    P.load(infile, outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(ALIGNMENT_TARGETS
           , regex("(\S+).dir/(\S+).bam")
           , r"\1.dir/\1_\2.picard_stats")
def buildPicardStats(infile, outfile):
    '''build alignment stats using picard.
    Note that picards counts reads but they are in fact alignments.
    '''
    if PARAMS["pool_reads"]:
        reffile = os.path.join(os.path.dirname(infile),"agg-agg-agg.filtered.contigs.fa")
    else:
        reffile = P.snip(infile, ".bam") + ".fa"
    PipelineMappingQC.buildPicardAlignmentStats( infile, 
                                                 outfile,
                                                 reffile )

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
#@jobs_limit( 1, "db" )
@merge( buildPicardStats, "picard_stats.load" )
def loadPicardStats( infiles, outfile ):
    '''merge alignment stats into single tables.'''

    PipelineMappingQC.loadPicardAlignmentStats( infiles, outfile )

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@follows(*ALIGNMENT_TARGETS)
@transform(ALIGNMENT_TARGETS,
           suffix(".bam"),
           add_inputs(loadPicardStats),
           ".coverage.gz")
def buildCoverageOverContigs(infiles, outfile):
    '''
    build histograms of the coverage over each of the contigs
    '''
    to_cluster = True
    
    bam = infiles[0]
    track = os.path.dirname(bam)[:-len(".dir")]+"_"+ P.snip(os.path.basename(bam), ".bam")

    # nnect to database
    dbh = connect()
    cc = dbh.cursor()

    # get number of passed filter aligned reads from picard stats
    scale_factor = cc.execute("""SELECT PF_READS_ALIGNED FROM picard_stats_alignment_summary_metrics
                              WHERE track == '%s'""" % track).fetchone()[0]

    scale_factor = 1/(float(scale_factor)/1000000)

    statement = '''genomeCoverageBed -ibam %(bam)s -scale %(scale_factor)f -d | gzip > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildCoverageOverContigs, suffix(".gz"), ".stats.gz")
def buildCoverageStats(infile, outfile):
    '''
    build coverage statistics - mean and standard deviation
    '''
    job_options = " -l mem_free=30G"
    to_cluster = True
    statement = '''zcat %(infile)s | python %(scriptsdir)s/coverage2stats.py --log=%(outfile)s.log | gzip > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildCoverageStats
           , suffix(".gz")
           , add_inputs(buildContigLengths)
           , ".postprocess.gz")
def postprocessCoverageStats(infiles, outfile):
    '''
    genomeCoverageBed outputs only non-zero depth. Add a "0" to 
    contigs that have zero coverage
    '''
    stats_file =  infiles[0]
    inf = IOTools.openFile(stats_file)
    header = inf.readline()

    if PARAMS["pool_reads"]:
        contigs = [x for x in infiles[1:len(infiles)] if x.find(os.path.dirname(stats_file)) != -1][0]
    else:
        contigs = stats_file.replace(".coverage.stats.gz", ".lengths.tsv")

    contig2stats = {}
    for line in inf.readlines():
        data = line[:-1].split("\t")
        contig, mean, sd = data[0], data[1], data[2]
        contig2stats[contig] = (mean, sd)

    inf2 = open(contigs)
    header2 = inf2.readline()
    outf = gzip.open(outfile, "w")
    outf.write(header)
    for line in inf2.readlines():
        data = line[:-1].split("\t")
        contig, length = data[0], data[1]
        if contig in contig2stats.keys():
            outf.write("%s\t%s\t%s\n" % (contig, contig2stats[contig][0], contig2stats[contig][1]))
        else:
            outf.write("%s\t0\t0\n" % contig)
    outf.close()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@jobs_limit(1, "db")
@transform(postprocessCoverageStats, suffix(".postprocess.gz"), ".load")
def loadCoverageStats(infile, outfile):
    '''
    load coverage stats
    '''
    tablename = P.toTable(P.snip(os.path.dirname(infile), ".dir") + "_%s" % os.path.basename(outfile))
    statement = '''zcat %(infile)s | python %(scriptsdir)s/csv2db.py 
                -t %(tablename)s 
                --index=contig 
                --log=%(outfile)s.log > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@follows(loadCoverageStats)
def coverage():
    pass

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@follows(loadEssentialGeneAssignments
         , loadLCA
         , loadContigGCContent
         , loadContigLengths
         , loadCoverageStats
         , loadEssentialGeneAssignments)
def loadContigAttributes():
    pass

####################
# full targets
####################
@follows( loadReadCounts
          , contig_stats
          , metaphlan
          , loadContigAttributes
          , coverage
          , loadGeneTables)
def full():
    pass

####################
# report building
####################
@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''
    E.info( "starting documentation build process from scratch" )
    P.run_report( clean = True )

@follows( mkdir( "report" ) )
def update_report():
    '''update report.'''
    E.info( "updating documentation" )
    P.run_report( clean = False )


if __name__== "__main__":
    if sys.argv[1] == "plot":
        pipeline_printout_graph("test.pdf", "pdf",[full], no_key_legend=True,
                                size = (4,4),
                                user_colour_scheme = {"colour_scheme_index" :1})
    else:
        sys.exit( P.main(sys.argv) )
    



