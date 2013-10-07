###############################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 
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
==========================
metagenome simulation
bechmark pipeline
==========================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.

Description
============

The pipeline is used to assess the results of a simulated metagenomic analysis using the
pipeline_metagenomeassembly.py pipeline. Configuration is therefore required to
point to the data that has been analysed. 

Results of this pipeline are used for simulated data where we are aware of the ground
truth. Three main features of the assembly are analysed:

* Taxonomic groups analysis - includes abundance estimations and false positive and
  negative rates

* Genome coverage - using the genomes that were used for simulations we align contigs
  back to the references and assess the proportion of eah genome that is covered
  by contigs with >99% similarity

* Assessing chimericity of assembled contigs

TODO
----

* Identifying alternative useful criteria for the benchmarking of simulated metagenome
  assemblies.

* gene annotations - given that we have a set of known protein coding genes from teh 
  reference annotations we can assess the quality of de novo gene predictions made
  through the pipeline_metagenomeassembly.py pipleline.

Code
====

"""

# load modules
from ruffus import *

import time
import CGAT.Experiment as E
import logging as L
import CGAT.Database as Database
import CGAT.CSV as CSV

import sys
import os
import re
import shutil
import itertools
import math
import glob
import time
import gzip
import collections
import random

import numpy
import sqlite3
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError
import CGATPipelines.PipelineMapping as PipelineMapping
import CGAT.FastaIterator as FastaIterator
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMappingQC as PipelineMappingQC
import CGAT.Bed as Bed
import CGAT.Nucmer as Nucmer
import pysam
import CGAT.Fastq as Fastq
import sqlite3
import CGATPipelines.PipelineMetagenomeBenchmark as PipelineMetagenomeBenchmark
import CGATPipelines.PipelineMetagenomeAssembly as PipelineMetagenomeAssembly

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

SEQUENCE_FILES = glob.glob("*.fastq.1.gz*") + glob.glob("*.fastq.gz")
GENOMES = glob.glob(os.path.join(PARAMS["genomes_genomesdir"], "*.fna"))
CONTIGS = glob.glob("*.fa")
ALIGNMENTS = glob.glob("*.bam")

######################################################
######################################################
######################################################
def dbList(xset):
    '''
    return a list for checking inclusion from a db query
    '''
    return "(" + ",".join(["'" + x + "'" for x in xset]) + ")"

######################################################
######################################################
######################################################
# The first benchmarking is to do with the estimated
# taxonomic abundances - metaphlan predicts these
# but we know the ground truth based on our simulations
######################################################
######################################################
######################################################
@follows(mkdir("taxonomy.dir"))
@merge(GENOMES
       , r"taxonomy.dir/gi_accessions.tsv")
def buildGiAccessionNumbers(infiles, outfile):
    '''
    build a list of gi accession numbers from the genomes
    in order to get assignments of taxa from them
    '''
    # first line in each file is the contig name
    outf = open(outfile, "w")
    for inf in infiles:
        outf.write(open(inf).readline().split("|")[1] + "\n")
    outf.close()
        
######################################################
######################################################
######################################################
@merge([buildGiAccessionNumbers] + [os.path.join(PARAMS["taxonomy_taxdir"], x) for x in ["ncbi.map", "gi_taxid_nucl.dmp.gz", "ncbi.lvl", "nodes.dmp"]]
       , "taxonomy.dir/gi2taxa.tsv")
def buildGi2Taxa(infiles, outfile):
    '''
    associate each input genome gi identifier to
    different taxonomic levels
    '''
    to_cluster = True
    gi_accessions = infiles[0]
    ncbi_map = infiles[1]
    gi2taxid_map = infiles[2]
    codes = infiles[3]
    nodes = infiles[4]
    statement = '''python %(scriptsdir)s/gi2parents.py 
                   -g %(gi_accessions)s
                   -m %(ncbi_map)s
                   -n %(gi2taxid_map)s
                   -c %(codes)s
                   -t %(nodes)s 
                   --log=%(outfile)s.log > %(outfile)s'''
    P.run()

######################################################
######################################################
######################################################
@transform(buildGiAccessionNumbers, suffix(".tsv"), ".load")
def loadGiAccessionNumbers(infile, outfile):
    '''
    load the gi acesssion numbers
    '''
    P.load(infile, outfile)

###################################################
###################################################
###################################################
@transform(SEQUENCE_FILES
           , regex("(\S+).(fastq.gz|fastq.1.gz)")
           , add_inputs(buildGi2Taxa)
           , r"taxonomy.dir/\1.taxonomy.relab")
def buildTrueTaxonomicRelativeAbundances(infiles, outfile):
    '''
    get species level relative abundances for the simulateds
    data. This involes creating maps between different identifiers
    from the NCBI taxonomy. This is so that the results are comparable
    to species level analysis from metaphlan
    The gi_taxid_nucl is a huge table and therefore this function
    takes an age to run - can think of optimising this somehow
    '''
    to_cluster = True
    PipelineMetagenomeBenchmark.buildTrueTaxonomicRelativeAbundances(infiles, outfile)

###################################################
###################################################
###################################################
@jobs_limit(1, "db")
@transform(buildTrueTaxonomicRelativeAbundances, suffix(".relab"), ".relab.load")
def loadTrueTaxonomicAbundances(infile, outfile):
    '''
    load the true taxonomic relative abundances
    '''
    P.load(infile, outfile)
    
###################################################
###################################################
###################################################
@follows(mkdir("taxonomy.dir"))
@transform(glob.glob(os.path.join(os.path.join(PARAMS["results_resultsdir"], "metaphlan.dir"), "*.relab"))
       , regex("(\S+)/(\S+).relab")
       , r"taxonomy.dir/metaphlan_\2.taxonomy.relab.load")
def loadEstimatedTaxonomicRelativeAbundances(infile, outfile):
    '''
    load metaphlan taxonomic abundance estimations
    '''
    P.load(infile, outfile)

###################################################
###################################################
###################################################
@transform(loadTrueTaxonomicAbundances
           , suffix(".load")
           , add_inputs(loadEstimatedTaxonomicRelativeAbundances)
           , ".png")
def plotRelativeAbundanceCorrelations(infiles, outfile):
    '''
    plot the correlation between the estimated 
    relative abundance of species and the true
    relative abundances - done on the shared set
    '''
    PipelineMetagenomeBenchmark.plotRelativeAbundanceCorrelations(infiles, outfile)

###################################################
###################################################
###################################################
@transform(loadTrueTaxonomicAbundances
           , suffix(".load")
           , add_inputs(loadEstimatedTaxonomicRelativeAbundances)
           , ".fp")
def calculateFalsePositiveRate(infiles, outfile):
    '''
    calculate the false positive rate in taxonomic
    abundances
    '''
    PipelineMetagenomeBenchmark.calculateFalsePositiveRate(infiles, outfile)

###################################################
###################################################
###################################################
@transform(calculateFalsePositiveRate, suffix(".fp"), ".fp.load")
def loadFalsePositiveRate(infile, outfile):
    '''
    load false positive rate calculations
    '''
    P.load(infile, outfile)

###################################################
###################################################
###################################################
@follows(loadFalsePositiveRate
         , plotRelativeAbundanceCorrelations)
def taxonomy():
    pass

###################################################
###################################################
###################################################
# look at assembly statistics after filtering 
# for contigs that are above a certain coverage
###################################################
###################################################
###################################################
COVERAGE_FILES = glob.glob(os.path.join(PARAMS["results_resultsdir"], "*/*.coverage.load"))
@follows(mkdir("contig_stats.dir"))
@transform(CONTIGS, regex("(\S+).fa")
           , add_inputs(*COVERAGE_FILES)
           , r"contig_stats.dir/\1.filtered.fa")
def filterContigsByCoverage(infiles, outfile):
    '''
    filter contigs by their average base coverage
    '''
    P.submit("PipelineMetagenomeBenchmark", "filterByCoverage", infiles = infiles, outfiles = outfile)

###################################################
###################################################
###################################################
@transform(filterContigsByCoverage, suffix(".fa"), ".stats")
def buildFilteredContigStats(infile, outfile):
    '''
    build contig stats for the filtered set
    '''
    PARAMS["filter"] = None
    PARAMS["scaffold_n"] = 50
    PipelineGenomeAssembly.contig_to_stats(infile, outfile, PARAMS)

###################################################
###################################################
###################################################
@transform(filterContigsByCoverage, suffix(".fa"), ".lengths.tsv")
def buildFilteredContigLengths(infile, outfile):
    '''
    output lengths for each contig in each of the assemblies
    '''
    PARAMS["filter"] = None
    PipelineGenomeAssembly.build_scaffold_lengths(infile, outfile, PARAMS)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildFilteredContigLengths, suffix(".lengths.tsv"), ".lengths.load")
def loadFilteredContigLengths(infile, outfile):
    '''
    load contig lengths
    '''
    outname = P.snip(os.path.dirname(infile), ".dir") + "_" + P.snip(os.path.basename(infile), ".tsv") + ".load"
    P.load(infile, outname)
    P.touch(outfile)

###################################################
###################################################
###################################################
@transform(buildFilteredContigStats, suffix(".filtered.stats"), ".filtered.load")
def loadFilteredContigStats(infile, outfile):
    '''
    load the filtered contig stats
    '''
    P.load(infile, outfile)

###################################################
###################################################
###################################################
def alignmentTargets(genome_files, contig_files):
    '''
    generator object to produce filenames for 
    aligning contigs to known ncbi genomes
    '''
    parameters = []
    for genome, contig in itertools.product(genome_files, contig_files):
        outfile = os.path.join("alignment.dir", P.snip(contig, ".contigs.fa") + "_vs_"  + P.snip(os.path.basename(genome), ".fna")) + ".delta"
        parameters.append([genome, outfile, contig])
    return parameters

###################################################
###################################################
###################################################
@follows(mkdir("alignment.dir"))
@files(alignmentTargets(GENOMES, CONTIGS))
def alignContigsToReference(infile, outfile, param):
    '''
    align the contigs to the reference genomes
    using nucmer
    '''
    print infile, param

    to_cluster = True

    reffile, contigfile = infile, param
    pattern = P.snip(os.path.basename(outfile), ".delta")
    statement = '''nucmer -p %(pattern)s %(reffile)s %(contigfile)s'''
    P.run()
    outf = os.path.basename(outfile)
    statement = '''mv %(outf)s alignment.dir'''
    P.run()

###################################################
###################################################
###################################################
@transform(alignContigsToReference, regex("alignment.dir/(\S+).delta"), r"alignment.dir/\1.filtered")
def filterAlignments(infile, outfile):
    '''
    filter alignments to retain only those that
    have > 99% identity to the reference
    '''
    to_cluster = True
    statement = '''delta-filter -q -i 99 %(infile)s > %(outfile)s'''
    P.run()

###################################################
###################################################
###################################################
@transform(filterAlignments, suffix(".filtered"), ".coords")
def buildAlignmentCoordinates(infile, outfile):
    '''
    build coordinates file from alignment delta
    file
    '''
    to_cluster = True
    statement = '''show-coords -T -r %(infile)s > %(outfile)s'''
    P.run()

###################################################
###################################################
###################################################
@transform(buildAlignmentCoordinates, suffix(".coords"), ".bed.gz")
def createAlignmentBedFiles(infile, outfile):
    '''
    create bed files - the intervals are with respect to the 
    reference genome
    intervals are merged to form a non redundant alignment set
    '''
    # has to be output from show coords in tab format
    # also have to be sorted for mergeBed
    to_cluster = True
    statement = '''cat %(infile)s
                   | python %(scriptsdir)s/nucmer2bed.py -t bed4 --log=%(outfile)s.log 
                   | mergeBed -i - 
                   | gzip > %(outfile)s'''
    P.run()

###################################################
###################################################
###################################################
@merge(createAlignmentBedFiles, "alignment.dir/alignments.size")
def buildAlignmentSizes(infiles, outfile):
    '''
    use bed files to sum the total number of bases
    that are aligned to the genomes
    '''
    outf = open(outfile, "w")
    outf.write("genome\tsize\n")
    for infile in infiles:
        genome = P.snip(os.path.basename(infile), ".bed.gz")
        c = 0
        inf = IOTools.openFile(infile)
        for bed in Bed.iterator(inf):
            c += bed.end - bed.start
        outf.write("%s\t%s\n" % (genome, str(c)))
    outf.close()

###################################################
###################################################
###################################################
@follows(mkdir("genome_stats.dir"))
@transform(glob.glob(os.path.join(PARAMS["genomes_genomesdir"], "*.fna"))
           , regex("(\S+)/(\S+).fna")
           , r"genome_stats.dir/\2.size")
def collectGenomeSizes(infile, outfile):
    '''
    output the genome sizes for each genome
    '''
    to_cluster = True
    outf = open(outfile, "w")
    outf.write("genome\tlength\n")
    # assume single fasta entry
    for fasta in FastaIterator.iterate(IOTools.openFile(infile)):
        name = P.snip(os.path.basename(infile), ".fna")
        length = len(list(fasta.sequence))
        outf.write("%s\t%s\n" % (name, str(length)))
    outf.close()

###################################################
###################################################
###################################################
@merge(collectGenomeSizes, "genome_stats.dir/genome.sizes")
def mergeGenomeSizes(infiles, outfile):
    '''
    merge the genome sizes into a summary file
    '''
    all_data = []
    outf = open(outfile, "w")
    header = open(infiles[0]).readline()
    outf.write(header)
    for infile in infiles:
        inf = open(infile)
        header = inf.readline()
        outf.write(inf.readline())
    outf.close()

###################################################
###################################################
###################################################
@follows(mkdir("expected_genome_coverage.dir"))
@transform(SEQUENCE_FILES
           , regex("(\S+).(fastq.gz|fastq.1.gz)")
           , add_inputs(mergeGenomeSizes)
           , r"expected_genome_coverage.dir/\1.expected.coverage")
def buildExpectedGenomeCoverage(infiles, outfile):
    '''
    build the expected coverage over the genomes
    in the sample based on read depth and length
    '''
    P.submit("PipelineMetagenomeBenchmark", "buildExpectedCoverageOverGenomes", infiles = infiles, outfiles = outfile)


###################################################
###################################################
###################################################
@follows(mkdir("genome_coverage.dir"))
@merge([buildAlignmentSizes, mergeGenomeSizes], "genome_coverage.dir/genome.coverage")
def buildCoverageOverGenomes(infiles, outfile):
    '''
    create file with the coverage over each of the 
    simulated genomes
    '''
    P.submit("PipelineMetagenomeBenchmark", "buildCoverageOverGenomes", infiles = infiles, outfiles = outfile)

###################################################
###################################################
###################################################
@transform(buildExpectedGenomeCoverage
           , regex("(\S+)/(\S+).coverage")
           , add_inputs(buildCoverageOverGenomes)
       , r"genome_coverage.dir/\2_observed.tsv")
def mergeExpectedAndObservedGenomeCoverage(infiles, outfile):
    '''
    merge the expected and actual estimates
    of genome coverage
    '''
 
    expected = open(infiles[0])
    expected_header = expected.readline()
    observed = open(infiles[1])
    observed_header = observed.readline()

    expected_data = {}
    E.info("reading expected coverage over genomes")
    for line in expected.readlines():
        data = line[:-1].split("\t")
        gi, coverage = data[0], data[1]
        expected_data[gi] = coverage

    outf = open(outfile, "w")
    E.info("writing results")
    outf.write("track\tgi\tspecies\tobserved\texpected\n")
    for line in observed.readlines():
        data = line[:-1].split("\t")
        track, gi, species, coverage = data[0], data[1], "_".join(data[2].split("_")[5:7]), data[3]
        outf.write("%s\t%s\t%s\t%s\t%s\n" % (track, gi, species, coverage, expected_data[gi]))
    outf.close()

###################################################
###################################################
###################################################
@transform(mergeExpectedAndObservedGenomeCoverage, suffix(".tsv"), ".load")
def loadExpectedAndObservedGenomeCoverage(infile, outfile):
    '''
    load the combined table for observed and expected
    genome coverage
    '''
    P.load(infile, outfile)

###################################################
###################################################
###################################################
@split(buildCoverageOverGenomes, "genome_coverage.dir/*.coverage.png")
def plotCoverageOverGenomes(infile, outfile):
    '''
    plot the percent coverage over each genome
    '''
    PipelineMetagenomeBenchmark.plotCoverageOverGenomes(infile, outfile)

###################################################
###################################################
###################################################
# assess chimeric contigs - Chimeric contigs are
###################################################
###################################################
###################################################
def chimeraTargets(alignment_files, contig_files):
    '''
    generator object to produce filenames for 
    scoring chimericity
    '''
    parameters = []
    for alignment, contig in itertools.product(genome_files, contig_files):
        outfile = os.path.join("chimeras.dir", P.snip(alignment, ".bam") + ".chimeras")
        parameters.append( [outfile, alignment, contig] )
    return parameters

###################################################
###################################################
###################################################
@follows(mkdir("expected_contigs.dir"))
@transform(ALIGNMENTS
           , regex("(\S+).bam")
           , add_inputs(CONTIGS)
           , r"expected_contigs.dir/\1.species_map.tsv")
def buildSpeciesMap(infiles, outfile):
    '''
    build species map file for input into
    contigs2random_samples.py
    '''
    to_cluster = True
    bam = infiles[0]
    contig = [x for x in infiles[1] if P.snip(x, ".contigs.fa") == P.snip(bam, ".bam")][0]
    statement = ''' cat %(contig)s | python %(scriptsdir)s/bam2species_map.py -b %(bam)s --log=%(outfile)s.log > %(outfile)s'''
    P.run()

###################################################
###################################################
###################################################
@transform(buildSpeciesMap, suffix(".tsv"), ".renamed.tsv")
def renameGenomesInSpeciesMap(infile, outfile):
    '''
    read names need to be converted to correct species genome
    names for input into buildExpectedContigs
    '''
    to_cluster = True
    genomes_dir = PARAMS["genomes_genomesdir"]
    statement = '''cat %(infile)s | python %(scriptsdir)s/species_map2species_map.py 
                   -g %(genomes_dir)s 
                   --log=%(outfile)s.log > %(outfile)s'''
    P.run()

###################################################
###################################################
###################################################
@transform(CONTIGS, regex("(\S+).fa")
           , add_inputs(renameGenomesInSpeciesMap)
           , r"expected_contigs.dir/\1.expected.fa")
def buildExpectedContigs(infiles, outfile):
    '''
    build an expectation set of contigs - sample from
    original genomes. Provides an opportunity 
    to assess the expected chimericity based on the 
    species that are present in the sample
    '''
    to_cluster = False
    contig = infiles[0]

    inputs = [infile for infile in infiles if infile != contig]
    species_map = [infile for infile in inputs if P.snip(os.path.basename(infile), ".species_map.renamed.tsv") == P.snip(contig, ".contigs.fa")][0]

    genomes_dir = PARAMS["genomes_genomesdir"]
    statement = '''cat %(contig)s | python %(scriptsdir)s/contigs2random_sample.py 
                   -m %(species_map)s 
                   -g %(genomes_dir)s 
                   --log=%(outfile)s.log > %(outfile)s'''
    P.run()

###################################################
###################################################
###################################################
@transform(buildExpectedContigs, suffix(".fa"), ".bt2")
def buildBowtie2Indices(infile, outfile):
    '''
    build bowtie indices
    '''
    to_cluster = True
    outbase = P.snip(infile, ".fa")
    statement = '''bowtie2-build -f %(infile)s %(outbase)s'''
    P.run()
    P.touch(outfile)

###################################################
###################################################
###################################################
@follows(buildBowtie2Indices)
@transform(buildBowtie2Indices
           , regex("(\S+)/(\S+).bt2")
           , add_inputs(SEQUENCE_FILES)
           , r"\1/\2.bam")
def mapReadsWithBowtieAgainstExpectedContigs(infiles, outfile):
    '''
    map reads against contigs with bowtie
    '''
    to_cluster = True
    
    bowtie_index_dir = "expected_contigs.dir"
    genome = os.path.basename(re.search(".*R[0-9]*", infiles[0]).group(0) + ".contigs.expected")
    for seq in infiles[1]:
        to_cluster = True
        infile, reffile = seq,  genome + ".fa"
        m = PipelineMapping.Bowtie( executable = P.substituteParameters( **locals() )["bowtie_executable"] )
        statement = m.build( (infile,), outfile ) 
        P.run()

###################################################
###################################################
###################################################
@follows(mkdir("chimeras.dir"))
@transform(ALIGNMENTS, regex("(\S+).bam"), r"chimeras.dir/\1.bam")
def linkAlignmentFiles(infile, outfile):
    '''
    not ideal but makes the next task easier to implement
    '''
    statement = '''cd chimeras.dir; checkpoint; ln -s ../%(infile)s .; checkpoint ln -s ../%(infile)s.bai .'''
    P.run()

###################################################
###################################################
###################################################
@transform([linkAlignmentFiles, mapReadsWithBowtieAgainstExpectedContigs]
           , regex("(\S+)/(\S+).bam")
           , r"chimeras.dir/\2.chimera")
def buildChimerasBasedOnReads(infile, outfile):
    '''
    this function is an alternative to counting a contig as a chimera
    if it aligns to more than one genome. A contig is likely to align
    to multiple genomes with high idenitity if there contains very similar
    genomes in the sample. This is true of our simulation that contains
    subspecies of the same species e.g. B.fragilis subspecies

    A more appropriate method for assessing chimericity is to score
    each contig with a chimericity score. The chimericity score is the 
    ratio of "good" alignments / "bad" alignments. An alignment is considered
    "good" if it is from the species from which the majority of alignments
    from that contig are derived'''

    P.submit("PipelineMetagenomeBenchmark", "buildChimerasBasedOnReads"
             , infiles = infile, outfiles = outfile)
        
###################################################
###################################################
###################################################
@merge(buildChimerasBasedOnReads, "chimeras.dir/observed_expected.chimera")
def mergeChimeras(infiles, outfile):
    '''
    merge chimeras to have all expected and observed in
    a single file
    '''
    pass


###################################################
###################################################
###################################################
@transform(buildChimerasBasedOnReads, suffix(".chimera"), ".chimera.load")
def loadChimericityScores(infile, outfile):
    '''
    load the chimericity scores
    '''
    P.load(infile, outfile)



###################################################
###################################################
###################################################
@follows(taxonomy
         , loadExpectedAndObservedGenomeCoverage
         , loadChimericityScores)
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
    sys.exit( P.main(sys.argv) )










