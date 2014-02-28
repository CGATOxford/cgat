import sys
import re
import os
import tempfile
import collections
import shutil
import gzip
import sqlite3

import CGAT.IOTools as IOTools
import CGAT.Pipeline as P
import CGAT.Experiment as E
import CGATPipelines.PipelineMapping as PipelineMapping
import CGAT.FastaIterator as FastaIterator
import CGAT.Fastq as Fastq
import glob
import collections
import CGATPipelines.PipelineTracks as PipelineTracks
import CGAT.Metaphlan as Metaphlan
import numpy as np
import shutil
import pysam
import sqlite3
from rpy2.robjects import r as R

P.getParameters(
    "pipeline.ini")

PARAMS = P.PARAMS


################################################
################################################
################################################
def buildTrueTaxonomicRelativeAbundances(infiles, outfile):
    '''
    get species level relative abundances for the simulateds
    data. This involes creating maps between different identifiers
    from the NCBI taxonomy. This is so that the results are comparable
    to species level analysis from metaphlan
    '''
    levels = ["species", "genus", "family", "order", "class", "phylum"]
    taxa = open(infiles[1])
    header = taxa.readline()
    gi2taxa = collections.defaultdict(list)
    for line in taxa.readlines():
        data = line[:-1].split("\t")
        gi, strain, species, genus, family, order, _class, phylum = data[
            0], data[1], data[2], data[3], data[4], data[5], data[6], data[7]
        gi2taxa[gi] = (species, genus, family, order, _class, phylum)

    outf = open(outfile, "w")
    outf.write("level\ttaxa\trelab\n")
    for i in range(len(levels)):
        total = 0
        result = collections.defaultdict(int)
        for fastq in Fastq.iterate(IOTools.openFile(infiles[0])):
            total += 1
            gi = fastq.identifier.split("|")[1]
            result[gi2taxa[gi][i]] += 1
        for taxa, value in result.iteritems():
            outf.write("%s\t%s\t%s\n" %
                       (levels[i], taxa, float(value) / total))
    outf.close()

################################################
################################################
################################################


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
    temp.write("true\testimate\tlevel\n")
    for estimate_file in infiles[1:]:
        if os.path.basename(estimate_file)[len("metaphlan_"):] == os.path.basename(true_file):
            tablenames = [P.toTable(os.path.basename(true_file)), P.toTable(
                os.path.basename(estimate_file))]
            # get data for each taxonomic level
            for taxa in ["phylum", "class", "order", "family", "genus", "species"]:
                statement = """SELECT a.relab, b.rel_abundance, a.level
                           FROM %s as a, %s as b
                           WHERE b.taxon_level == "%s"
                           AND a.taxa == b.taxon""" % (tablenames[0], tablenames[1], taxa)
                for data in cc.execute(statement).fetchall():
                    true, estimate, level = data[0], data[1], data[2]

                    temp.write("%f\t%f\t%s\n" % (true, estimate, level))
    temp.close()

    inf = temp.name
    R('''data <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' %
      inf)
    R('''library(ggplot2)''')
    R('''data$estimate <- data$estimate/100''')
    R('''ggplot(data, aes(true, estimate, colour = level)) + geom_point() + geom_smooth(method = "lm")''')
    R('''ggsave("%s")''' % outfile)

    out_cors = P.snip(outfile, ".pdf") + ".cors"
    R('''cors <- data.frame("level" = c("phylum", "class", "order", "family", "genus", "species"), "cor" = rep(0, 6))''')
    R('''cors[1,2] <- cor(data$true[data$level == "phylum"], data$estimate[data$level == "phylum"])''')
    R('''cors[2,2] <- cor(data$true[data$level == "class"], data$estimate[data$level == "class"])''')
    R('''cors[3,2] <- cor(data$true[data$level == "order"], data$estimate[data$level == "order"])''')
    R('''cors[4,2] <- cor(data$true[data$level == "family"], data$estimate[data$level == "family"])''')
    R('''cors[5,2] <- cor(data$true[data$level == "genus"], data$estimate[data$level == "genus"])''')
    R('''cors[6,2] <- cor(data$true[data$level == "species"], data$estimate[data$level == "species"])''')
    R('''write.table(cors, file = "%s", row.names = F, sep = "\t")''' %
      out_cors)

    # do the same at the low end - not for species
    R('''data$estimate <- data$estimate/100''')
    R('''ggplot(data[data$true < 0.75 & data$level != "species",], aes(true, estimate, colour = level)) + geom_point() + geom_smooth(method = "lm", se = F)''')
    outf = P.snip(outfile, ".pdf") + ".lowest.pdf"
    R('''ggsave("%s")''' % outf)

    out_cors = P.snip(outfile, ".pdf") + ".lowest.cors"
    R('''cors <- data.frame("level" = c("phylum", "class", "order", "family", "genus", "species"), "cor" = rep(0, 6))''')
    R('''cors[1,2] <- cor(data$true[data$level == "phylum" & data$true < 0.75], data$estimate[data$level == "phylum" & data$true < 0.75])''')
    R('''cors[2,2] <- cor(data$true[data$level == "class" & data$true < 0.75], data$estimate[data$level == "class" & data$true < 0.75])''')
    R('''cors[3,2] <- cor(data$true[data$level == "order" & data$true < 0.75], data$estimate[data$level == "order" & data$true < 0.75])''')
    R('''cors[4,2] <- cor(data$true[data$level == "family" & data$true < 0.75], data$estimate[data$level == "family" & data$true < 0.75])''')
    R('''cors[5,2] <- cor(data$true[data$level == "genus" & data$true < 0.75], data$estimate[data$level == "genus" & data$true < 0.75])''')
    R('''write.table(cors, file = "%s", row.names = F, sep = "\t")''' %
      out_cors)

    os.unlink(inf)

################################################
################################################
################################################


def calculateFalsePositiveRate(infiles, outfile):
    '''
    taxonomy false positives and negatives etc
    '''
    # connect to database
    dbh = sqlite3.connect(PARAMS["database"])
    cc = dbh.cursor()

    levels = ["phylum", "class", "order", "family", "genus", "species"]
    tablename_true = P.toTable(infiles[0])

    # get corresponding estimate file
    tablename_estimate = P.toTable(os.path.basename([inf for inf in infiles[
                                   1:] if os.path.basename(inf)[len("metaphlan_"):] == os.path.basename(infiles[0])][0]))

    outf = open(outfile, "w")
    track = P.snip(os.path.basename(infiles[0]), ".taxonomy.relab.load")
    for level in levels:
        for cutoff in [0, 1]:
            true_set = set()
            estimate_set = set()
            for taxa in cc.execute("""SELECT taxa FROM %s WHERE level == '%s' AND relab > %f""" % (tablename_true, level, float(cutoff) / 100)):
                true_set.add(taxa[0])
            for taxa in cc.execute("""SELECT taxon FROM %s WHERE taxon_level == '%s' AND rel_abundance > %f""" % (tablename_estimate, level, float(cutoff))):
                estimate_set.add(taxa[0])
            total_true = len(true_set)
            total_estimate = len(estimate_set)
            tp = true_set.intersection(estimate_set)
            fp = estimate_set.difference(true_set)

            fp_rate = float(len(fp)) / total_estimate
            tp_rate = float(len(tp)) / total_true
            outf.write("%s\t%f\t%f\t%s\t%s\n" %
                       (level, fp_rate, tp_rate, track, str(cutoff)))
    outf.close()

################################################
################################################
################################################


def filterByCoverage(infiles, outfile):

    fcoverage = PARAMS["coverage_filter"]
    contig_file = infiles[0]
    dbh = sqlite3.connect(
        os.path.join(PARAMS["results_resultsdir"], PARAMS["database"]))
    cc = dbh.cursor()
    contigs = set()
    for infile in infiles[1:]:
        dirsplit = infile.split("/")
        infile = os.path.join(
            PARAMS["results_resultsdir"], dirsplit[-2].split(".dir")[0] + "-" + dirsplit[-1])
        tablename = P.toTable(os.path.basename(infile))
        if P.snip(contig_file, ".fa") == P.snip(os.path.basename(infile), ".coverage.load"):
            statement = """SELECT contig_id ave FROM
                           (SELECT contig_id, AVG(coverage) as ave FROM %s GROUP BY contig_id)
                           WHERE ave > %i""" % (tablename, PARAMS["coverage_filter"])
            for data in cc.execute(statement).fetchall():
                contigs.add(data[0])
    outf = open(outfile, "w")
    print contigs
    for fasta in FastaIterator.iterate(IOTools.openFile(contig_file)):
        identifier = fasta.title.split(" ")[0]
        if identifier in contigs:
            outf.write(">%s\n%s\n" % (identifier, fasta.sequence))
    outf.close()

################################################
################################################
################################################


def buildChimerasBasedOnReads(infile, outfile):

    bamfile = infile

    E.info("reading bam file")
    samfile = pysam.Samfile(bamfile)

    result = collections.defaultdict(list)
    E.info("iterating over reads")
    for read in samfile.fetch():
        # get gi as genome identifier
        species = read.qname.split("|")[1]
        result[samfile.getrname(read.tid)].append(species)

    E.info("calculating proportion of 'good' reads")
    outf = open(outfile, "w")
    outf.write("contig\tchimericity\n")
    for contig in result.keys():
        ngood = collections.defaultdict(int)
        for species in result[contig]:
            ngood[species] += 1
        chimericity = float(
            1 - float(max(ngood.values())) / sum(ngood.values()))
        outf.write("%s\t%s\n" % (contig, str(chimericity)))
    outf.close()

################################################
################################################
################################################


def buildExpectedCoverageOverGenomes(infiles, outfile):
    '''
    take sequence files and estimate the theoretical
    coverage we would get over genomes in the 
    sample i.e. at 1X coverage
    '''

    # if paired end then will have to multiply
    # by two
    multiply = False
    if infiles[0].endswith(".fastq.1.gz"):
        multiply = True

    # the theoretical coverage is defined as
    # (read length (L) * no. reads (N)) / genome size (G) (bp)

    # get genome sizes into memory
    genomes = open(infiles[1])
    header = genomes.readline()
    genome_sizes = {}
    for line in genomes.readlines():
        data = line[:-1].split("\t")
        gi = data[0].split("_")[1]
        size = data[1]
        genome_sizes[gi] = size

    # get the expected genome size
    expected_genome_sizes = collections.defaultdict(int)
    E.info("iterating over fastq file")
    for fastq in Fastq.iterate(IOTools.openFile(infiles[0])):
        gi = fastq.identifier.split("|")[1]
        expected_genome_sizes[gi] += 1
    E.info("iterating over fastq file: DONE")

    # get the proportion of each genome covered
    outf = open(outfile, "w")
    outf.write("gi\texpected_coverage\n")
    for gi, size in expected_genome_sizes.iteritems():
        if multiply:
            size = size * 2
        if gi not in genome_sizes:
            E.warn("could not find gi no. %s in dictionary" % gi)
            continue
        proportion_coverage = float(size) / float(genome_sizes[gi])
        if proportion_coverage > 1:
            proportion_coverage = 1
        outf.write("%s\t%f\n" % (gi, proportion_coverage))
    outf.close()


################################################
################################################
################################################
def buildCoverageOverGenomes(infiles, outfile):
    '''
    create file with the coverage over each of the 
    simulated genomes
    '''
    outf = open(outfile, "w")
    outf.write("track\tgi\tgenome\tpcoverage\n")
    alignments = {}
    genomes = {}

    inf1 = open(infiles[0])
    inf2 = open(infiles[1])
    header2 = inf2.readline()
    header1 = inf1.readline()
    for line in inf2.readlines():
        data = line[:-1].split("\t")
        name, value = data[0], data[1]
        genomes[name] = value

    for line in inf1.readlines():
        data = line[:-1].split("\t")
        name, value = data[0], data[1]
        gname = name.split("_")
        gi = gname[3]
        gname = "_".join(gname[2:])
        name = name.split("_vs")[0]
        alignments[(name, gname, gi)] = (
            gname, float(value) / float(genomes[gname]))

    for track, proportion in alignments.iteritems():
        outf.write("%s\t%s\t%s\t%f\n" %
                   (track[0], track[2], proportion[0], proportion[1]))
    outf.close()

################################################
################################################
################################################


def plotCoverageOverGenomes(infile, outfile):
    '''
    plot the percent coverage over each genome
    '''
    names_list = set()
    inf = open(infile)
    header = inf.readline()
    for line in inf.readlines():
        name = line[:-1].split("\t")[0]
        names_list.add(name)

    R('''data <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' %
      infile)
    for name in names_list:
        outname = os.path.join("genome_coverage.dir", name) + ".coverage.png"
        R('''data2 <- data[data$track == "%s",]''' % name)
        R('''png("%s", height = 500, width = 1000)''' % outname)
        R('''par(mai = c(2,1,1,1))''')
        R('''barplot(data2$pcoverage[order(data2$pcoverage)]
             , names.arg = data2$gi[order(data2$pcoverage)]
             , las = 2
             , cex.names = 0.75
             , ylim = c(0,1)
             , col = "blue")''')
        R('''abline(h = 0.95, cex = 2)''')
        R('''abline(h = 0.90, cex = 2, lty = 2)''')
        R["dev.off"]()
