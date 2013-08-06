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

The pipeline is used to assess the results of a metagenomic analysis using the
pipeline_metagenomeassembly.py pipeline. Configuration is therefore required to
point to the data that has been analysed. 

Results of this pipeline are used for simulated data where we are aware of the ground
truth. Three main features of the assembly are analysed:

* Taxonomic groups analysis - includes abundance estimations and false positive and
  negative rates

* Genome coverage - using the genomes that were used for simulations we align contigs
  back to the references and assess the proportion of eah genome that is covered
  by contigs with >99% similarity

TODO
----

* gene annotations - given that we have a set of known protein coding genes from teh 
  reference annotations we can assess the quality of de novo gene predictions made
  through the pipeline_metagenomeassembly.py pipleline.

Code
====

"""

# load modules
from ruffus import *

import time
import Experiment as E
import logging as L
import Database, CSV

import sys, os, re, shutil, itertools, math, glob, time, gzip, collections, random

import numpy, sqlite3
import GFF, GTF, IOTools, IndexedFasta
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError
import PipelineMapping
import FastaIterator
import PipelineMapping
import PipelineMappingQC
import Bed
import Nucmer
import pysam
import Fastq
import sqlite3

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import Pipeline as P
P.getParameters( 
    "pipeline.ini" )

PARAMS = P.PARAMS

SEQUENCE_FILES = glob.glob("*.fastq.1.gz*")
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
@transform(os.path.join(PARAMS["taxonomy_taxdir"], "names.dmp")
           , regex("(\S+)/(\S+).dmp")
           , r"taxonomy.dir/\2.dmp.gz")
def parseTaxonomyNamesFile(infile, outfile):
    '''
    parse the ncbi taxonomy file
    '''
    inf = open(infile)
    header = inf.readline()
    outf = gzip.open(outfile, "w")
    outf.write("taxid\ttaxname\tdescription\n")
    for line in inf.readlines():
        data = "".join(line[:-1].split("\t")).split("|")
        taxid, name, description = data[0], data[1], data[3]
        outf.write("%s\t%s\t%s\n" % (taxid, name, description))
    outf.close()

######################################################
######################################################
######################################################
@transform(parseTaxonomyNamesFile, suffix(".dmp.gz"), ".load")
def loadTaxonomyNames(infile, outfile):
    '''
    load taxonomic names
    '''
    P.load(infile, outfile, "--index=taxid")

######################################################
######################################################
######################################################
@follows(mkdir("taxonomy.dir"))
@transform(os.path.join(PARAMS["taxonomy_taxdir"], "categories.dmp")
           , regex("(\S+)/(\S+).dmp")
           , r"taxonomy.dir/\2.dmp.gz")
def parseCategoriesFile(infile, outfile):
    '''
    parse categories file
    '''
    inf = open(infile)
    outf = gzip.open(outfile, "w")
    outf.write("kingdom\tspecies_id\ttaxid\n")
    for line in inf.readlines():
        data = line[:-1].split("\t")
        outf.write("%s\t%s\t%s\n" % (data[0], data[1], data[2]))
    outf.close()

######################################################
######################################################
######################################################
@follows(mkdir("taxonomy.dir"))
@transform(parseCategoriesFile, suffix(".dmp.gz"), ".load")
def loadCategoriesFile(infile, outfile):
    '''
    load the taxonomic categories file
    '''
    P.load(infile, outfile, "--index=taxid")

###################################################
###################################################
###################################################
@follows(mkdir("taxonomy.dir"))
@transform(os.path.join(PARAMS["taxonomy_taxdir"], "gi_taxid_nucl.dmp.gz")
           , regex("(\S+)/(\S+).dmp.gz")
           , r"taxonomy.dir/\2.load")
def loadGi2Taxid(infile, outfile):
    '''
    load gi to taxid mapping
    '''
    P.load(infile, outfile, "--header=gi,taxid --index=taxid")

###################################################
###################################################
###################################################
@transform(SEQUENCE_FILES
           , regex("(\S+).fastq.1.gz")
           , r"taxonomy.dir/\1.taxonomy.relab")
def buildTrueTaxonomicRelativeAbundances(infile, outfile):
    '''
    get species level relative abundances for the simulateds
    data. This involes creating maps between different identifiers
    from the NCBI taxonomy. This is so that the results are comparable
    to species level analysis from metaphlan
    The gi_taxid_nucl is a huge table and therefore this function
    takes an age to run - can think of optimising this somehow
    '''
    to_cluster = True

    total = 0
    rel_abundance = collections.defaultdict(int)
    for fastq in Fastq.iterate(IOTools.openFile(infile)):
        total += 1
        gi = fastq.identifier.split("|")[1]
        rel_abundance[gi] += 1
    for gi, ab in rel_abundance.iteritems():
        rel_abundance[gi] = float(ab)/total

    dbh = sqlite3.connect(PARAMS["database"])
    cc = dbh.cursor()
    result = collections.defaultdict(float)
    for gi in rel_abundance.keys():
        E.info("processing gi %s" % gi)
        taxid = cc.execute("""SELECT taxid FROM gi_taxid_nucl WHERE gi == '%s'""" % gi).fetchone()[0]
        species_id = cc.execute("""SELECT species_id FROM categories WHERE taxid == '%s'""" % taxid).fetchone()[0]
        species_name = cc.execute("""SELECT taxname FROM names WHERE taxid == '%s' AND description == 'scientific name'""" % species_id).fetchone()[0]
        abundance = rel_abundance[gi]
        E.info("mapped gi %s to taxid: %s, species_id: %s, species_name: %s" % (str(gi), str(taxid), str(species_id), species_name))
        result[species_name] += abundance

    outf = open(outfile, "w")
    outf.write("species_name\trelab\n")
    for species_name, abundance in result.iteritems():
        # create names consistent with metaphlan
        species_name = species_name.replace(" ", "_")
        outf.write("%s\t%f\n" % (species_name, abundance))
    outf.close()

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
    # connect to database
    dbh = sqlite3.connect(PARAMS["database"])
    cc = dbh.cursor()

    true_file = infiles[0]
    temp = P.getTempFile()
    temp.write("true\testimate\n")
    for estimate_file in infiles[1:]:
        if os.path.basename(estimate_file)[len("metaphlan_"):] == os.path.basename(true_file):
            tablenames = [P.toTable(os.path.basename(true_file)), P.toTable(os.path.basename(estimate_file))]
            # get data
            statement = """SELECT a.relab, b.rel_abundance
                           FROM %s as a, %s as b
                           WHERE b.taxon_level == "species"
                           AND a.species_name == b.taxon""" % (tablenames[0], tablenames[1])
            for data in cc.execute(statement).fetchall():
                true, estimate = data[0], data[1]
                temp.write("%f\t%f\n" % (true, estimate))
    temp.close()
    print temp.name

    inf = temp.name
    R('''data <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % inf)
    R('''png("%s")''' % outfile)
    main_name = P.snip(outfile, ".png")
    R('''data$estimate <- data$estimate/100''')
    R('''plot(data$estimate, data$true, pch = 16, main = "%s", xlab = "estimated relative abundance", ylab = "observed relative abundance")''' % main_name)
    R('''text(0.05, y = 0.35, labels = paste("r = ", round(cor(data$estimate, data$true),2)), cex = 2)''')
    R["dev.off"]()
    os.unlink(inf)

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

    # connect to database
    dbh = sqlite3.connect(PARAMS["database"])
    cc = dbh.cursor()

    true_file = infiles[0]
    true_set = set()
    estimate_set = set()
    for estimate_file in infiles[1:]:
        if os.path.basename(estimate_file)[len("metaphlan_"):] == os.path.basename(true_file):
            tablenames = [P.toTable(os.path.basename(true_file)), P.toTable(os.path.basename(estimate_file))]

            for species in cc.execute("""SELECT species_name FROM %s""" % tablenames[0]).fetchall():
                true_set.add(species[0])
            for species in cc.execute("""SELECT taxon FROM %s WHERE taxon_level == 'species'""" % tablenames[1]).fetchall():
                if species[0].find("_unclassified") != -1: continue
                estimate_set.add(species[0])
    
    total_estimate = len(estimate_set)
    total_true = len(true_set)

    E.info("counting false positives and false negatives")
    print estimate_set.difference(true_set)
    nfp = len(estimate_set.difference(true_set))
    nfn = len(true_set.difference(estimate_set))
    ntp = len(estimate_set.intersection(true_set))

    E.info("writing results")
    track = P.snip(os.path.basename(true_file), ".load")
    outf = open(outfile, "w")
    outf.write("track\ttp_rate\tfp_rate\tfn_rate\n")
    outf.write("\t".join(map(str, [track, float(ntp)/total_estimate, float(nfp)/total_estimate, float(nfn)/total_true])) + "\n")
    outf.close()

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
@follows(loadTaxonomyNames
         ,loadCategoriesFile
         ,loadGi2Taxid
         , loadFalsePositiveRate
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
COVERAGE_FILES = glob.glob(os.path.join(PARAMS["results_resultsdir"], "*/*.coverage"))
@follows(mkdir("contig_stats.dir"))
@transform(COVERAGE_FILES
           , regex("(\S+)/(\S+).dir/(\S+).coverage")
           , r"contig_stats.dir/\2_\3.coverage.load")
def loadCoverageData(infile, outfile):
    '''
    load coverage data into database
    '''
    to_cluster = True
    tablename = P.toTable(outfile)
    database = os.path.join(PARAMS["results_resultsdir"], PARAMS["database"])
    dbh = sqlite3.connect(database)
    cc = dbh.cursor()
    temp = P.getTempFile()
    temp.write("contig_id\tacoverage\n")
    for data in cc.execute("""SELECT contig_id, AVG(coverage) FROM %s GROUP BY contig_id""" % tablename).fetchall():
        temp.write("\t".join(list(data)) + "\n")
    temp.close()
    P.load(temp.name, outfile)
    os.unlink(temp.name)

###################################################
###################################################
###################################################
@transform(CONTIGS, regex("(\S+).fa")
           , add_inputs(loadCoverageData)
           , "contig_stats.dir/\1.filtered")
def filterContigsByCoverage(infiles, outfile):
    '''
    filter contigs by their average base coverage
    '''
    fcoverage = PARAMS["coverage_filter"]
    contig_file = infiles[0]
    dbh = sqlite3.connect(PARAMS["database"])
    cc = dbh.cursor()
    for infile in infiles[1:]:
        print contig_file, P.snip(os.path.basename(infile), ".load")




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
        additional_input = add_inputs(contig)
        parameters.append( [outfile, genome, contig] )
    return parameters

###################################################
###################################################
###################################################
@follows(mkdir("alignment.dir"))
@parallel(alignmentTargets(GENOMES, CONTIGS))
def alignContigsToReference(outfile, param1, param2):
    '''
    align the contigs to the reference genomes
    using nucmer
    '''
    to_cluster = True

    reffile, contigfile = param1, param2
    pattern = P.snip(os.path.basename(outfile), ".delta")
    statement = '''nucmer -p %(pattern)s %(reffile)s %(contigfile)s'''
    P.run()
    outf = os.path.basename(outfile)
    statement = '''mv %(outf)s alignment.dir'''
    P.run()

###################################################
###################################################
###################################################
@transform(alignContigsToReference, suffix(".delta"), ".filtered")
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
@follows(mkdir("genome_coverage.dir"))
@merge([buildAlignmentSizes, mergeGenomeSizes], "genome_coverage.dir/genome.coverage")
def buildCoverageOverGenomes(infiles, outfile):
    '''
    create file with the coverage over each of the 
    simulated genomes
    '''
    outf = open(outfile, "w")
    outf.write("genome\tpcoverage\n")
    alignments = {}
    genomes = {}
    for infile in infiles:
        inf = open(infile)
        header = inf.readline()
        for line in inf.readlines():
            data = line[:-1].split("\t")
            name, value = data[0], data[1]
            if infile.find("alignments") != -1:
                alignments[name] = value
            else:
                genomes[name] = value
    for name in alignments.keys():
        outf.write("%s\t%f\n" % (name, float(alignments[name])/float(genomes[name])))
    outf.close()

###################################################
###################################################
###################################################
@transform(buildCoverageOverGenomes, suffix(".coverage"), ".coverage.pdf")
def plotCoverageOverGenomes(infile, outfile):
    '''
    plot the percent coverage over each genome
    '''
    R('''data <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)''' % infile)
    # rename 
    R('''names <- unlist(strsplit(rownames(data), "NC_*_", perl = T))''')
    R('''names <- unlist(strsplit(names, "_complete"))''')
    R('''names <- names[seq(2, length(names), 3)]''')

    R('''rownames(data) <- names''')
    # plot
    R('''pdf("%s", height = 8, width = 15)''' % outfile)
    R('''par(mai = c(5,1,1,1))''')
    R('''barplot(data[,1][order(data[,1])]
         , names.arg = rownames(data)[order(data[,1])]
         ,las = 2
         , cex.names = 0.75
         , ylim = c(0,1)
         , col = "blue")''')
    R('''abline(h = 0.95, cex = 2)''')
    R('''abline(h = 0.90, cex = 2, lty = 2)''')
    R["dev.off"]()


############
# DEPRECATED
############
###################################################
###################################################
###################################################
# assess chimeric contigs - Chimeric contigs are
# defined as those contigs that have multiple
# exact matches at > 99% id with > 1 reference 
###################################################
###################################################
###################################################
# @follows(mkdir("chimeras.dir"))
# @merge([buildAlignmentCoordinates, PARAMS["contigs_contigs"]]
#        , "chimeras.dir/chimeras.tsv")
# def buildChimeraContigs(infiles, outfile):
#     '''
#     build a set of chimeric contigs. Defined as aligning
#     to > 1 genome with > 99% identity
#     '''
#     contigs = [infile for infile in infiles if not infile.endswith(".coords")][0]
#     coordfiles = [infile for infile in infiles if infile.endswith(".coords")]

#     E.info("reading contigs into memory")
#     contigs = [fasta.title for fasta in FastaIterator.iterate(IOTools.openFile(contigs))]
#     total_contigs = len(contigs)

#     E.info("reading coordinate files")
#     alignments = collections.defaultdict(set)
#     for coordfile in coordfiles:
#         name = P.snip(os.path.basename(coordfile), ".coords")
#         inf = open(coordfile)
#         for data in Nucmer.Coords()(inf):
#             alignments[name].add(data[-1])

#     E.info("checking chimeric contigs")
#     not_aligned = 0
#     result = collections.defaultdict(set)
#     for contig in contigs:
#         for name in alignments.keys():
#             if contig in alignments[name]:
#                 result[contig].add(name)

#     E.info("writing results")
#     outf = open(outfile, "w")
#     outf.write("contig\tgenome_count\n")
#     for contig, genome in result.iteritems():
#         outf.write("%s\t%s\n" % (contig, ", ".join(list(genome))))
#     outf.close()
    
# ###################################################
# ###################################################
# ###################################################
# @transform(buildChimeraContigs, suffix(".tsv"), ".count.load")
# def loadChimeraContigs(infile, outfile):
#     '''
#     load chimera counts
#     '''
#     P.load(infile, outfile)

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
@follows(mkdir("chimeras.dir"))
@transform(ALIGNMENTS, regex("(\S+).bam"), r"chimeras.dir/\1.chimera")
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

    bamfile = infile
    
    E.info("reading bam file")
    samfile = pysam.Samfile(bamfile)
    
    result = collections.defaultdict(list)
    E.info("iterating over reads")
    for read in samfile.fetch():
        species = read.qname.split("ln")[0].split("from")[1]
        species = species[1:-1]
        result[read.tid].append(species)

    E.info("calculating proportion of 'good' reads")
    outf = open(outfile, "w")
    outf.write("contig\tchimericity\n")
    for contig in result.keys():
        ngood = collections.defaultdict(int)
        for species in result[contig]:
            ngood[species] += 1
        chimericity = float(max(ngood.values()))/sum(ngood.values())
        outf.write("%s\t%s\n" % contig, (str(chimericity)))
    outf.close()
        

@follows(taxonomy
         , plotCoverageOverGenomes
         , buildChimerasBasedOnReads)
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










