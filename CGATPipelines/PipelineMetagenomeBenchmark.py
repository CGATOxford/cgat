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
    "pipeline.ini" )

PARAMS = P.PARAMS


################################################
################################################
################################################
def buildTrueTaxonomicRelativeAbundances(infile, outfile):
    '''
    get species level relative abundances for the simulateds
    data. This involes creating maps between different identifiers
    from the NCBI taxonomy. This is so that the results are comparable
    to species level analysis from metaphlan
    The gi_taxid_nucl is a huge table and therefore this function
    takes an age to run - can think of optimising this somehow
    '''
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

    inf = temp.name
    R('''data <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % inf)
    R('''png("%s")''' % outfile)
    main_name = P.snip(outfile, ".png")
    R('''data$estimate <- data$estimate/100''')
    R('''plot(data$estimate, data$true, pch = 16, main = "%s", xlab = "estimated relative abundance", ylab = "observed relative abundance")''' % main_name)
    R('''text(0.05, y = 0.35, labels = paste("r = ", round(cor(data$estimate, data$true),2)), cex = 2)''')
    R('''abline(lm(data$true~data$estimate))''')
    R["dev.off"]()
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

    true_file = infiles[0]
    full_true_set = set()
    full_estimate_set = set()
    pointone_true_set = set()
    pointone_estimate_set = set()
    one_true_set = set()
    one_estimate_set = set()
    five_true_set = set()
    five_estimate_set = set()
    ten_true_set = set()
    ten_estimate_set = set()

    for estimate_file in infiles[1:]:
        if os.path.basename(estimate_file)[len("metaphlan_"):] == os.path.basename(true_file):
            tablenames = [P.toTable(os.path.basename(true_file)), P.toTable(os.path.basename(estimate_file))]

            for species in cc.execute("""SELECT species_name, relab FROM %s""" % tablenames[0]).fetchall():
                full_true_set.add(species[0])
                if species[1] > 0.001:
                    pointone_true_set.add(species[0])
                if species[1] > 0.01:
                    one_true_set.add(species[0])
                if species[1] > 0.05:
                    five_true_set.add(species[0])
                if species[1] > 0.1:
                    ten_true_set.add(species[0])
    
            for species in cc.execute("""SELECT taxon, rel_abundance 
                                         FROM %s WHERE taxon_level == 'species' 
                                         AND rel_abundance""" % tablenames[1]).fetchall():
                if species[0].find("_unclassified") != -1: continue
                full_estimate_set.add(species[0])
                if species[1] > 0.1:
                    pointone_estimate_set.add(species[0])
                if species[1] > 1:
                    one_estimate_set.add(species[0])
                if species[1] > 5:
                    five_estimate_set.add(species[0])
                if species[1] > 10:
                    ten_estimate_set.add(species[0])

    total_full_true = len(full_true_set)
    total_full_estimate = len(full_estimate_set)

    total_pointone_true = len(pointone_true_set)
    total_pointone_estimate = len(pointone_estimate_set)

    total_one_true = len(one_true_set)
    total_one_estimate = len(one_estimate_set)

    total_five_true = len(five_true_set)
    total_five_estimate = len(five_estimate_set)

    total_ten_true = len(ten_true_set)
    total_ten_estimate = len(ten_estimate_set)

    E.info("counting false positives and false negatives")
    nfp_full = len(full_estimate_set.difference(full_true_set))
    nfn_full = len(full_true_set.difference(full_estimate_set))
    ntp_full = len(full_estimate_set.intersection(full_true_set))

    nfp_pointone = len(pointone_estimate_set.difference(pointone_true_set))
    nfn_pointone = len(pointone_true_set.difference(pointone_estimate_set))
    ntp_pointone = len(pointone_estimate_set.intersection(pointone_true_set))

    nfp_one = len(one_estimate_set.difference(one_true_set))
    nfn_one = len(one_true_set.difference(one_estimate_set))
    ntp_one = len(one_estimate_set.intersection(one_true_set))

    nfp_five = len(five_estimate_set.difference(five_true_set))
    nfn_five = len(five_true_set.difference(five_estimate_set))
    ntp_five = len(five_estimate_set.intersection(five_true_set))

    nfp_ten = len(ten_estimate_set.difference(ten_true_set))
    nfn_ten = len(ten_true_set.difference(ten_estimate_set))
    ntp_ten = len(ten_estimate_set.intersection(ten_true_set))

    out = ["full", "pointone", "one", "five", "ten"]
    outf = open(outfile, "w")
    outf.write("track\ttp_rate\tfp_rate\tfn_rate\n")
    for x in out:
        track = P.snip(os.path.basename(true_file), ".load") + "_%s" % x
        total_true = eval("total_%s_true" % x)
        total_estimate = eval("total_%s_estimate" % x)

        ntp = eval("ntp_%s" % x)
        nfp = eval("nfp_%s" % x)
        nfn = eval("nfn_%s" % x)
        
        tpr = float(ntp) / total_estimate
        fpr = float(nfp) / total_estimate
        fnr = float(nfn) / total_true
        
        outf.write("\t".join(map(str, [track, ntp, nfp, nfn])) + "\n")
    outf.close()

################################################
################################################
################################################
def filterByCoverage(infiles, outfile):

    fcoverage = PARAMS["coverage_filter"]
    contig_file = infiles[0]
    dbh = sqlite3.connect(os.path.join(PARAMS["results_resultsdir"],PARAMS["database"]))
    cc = dbh.cursor()
    contigs = set()
    for infile in infiles[1:]:
        dirsplit = infile.split("/")
        infile = os.path.join(PARAMS["results_resultsdir"], dirsplit[-2].split(".dir")[0] + "-" + dirsplit[-1])
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
        chimericity = float(1 - float(max(ngood.values()))/sum(ngood.values()))
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
            size = size*2
        if gi not in genome_sizes:
            E.warn("could not find gi no. %s in dictionary" % gi)
            continue
        proportion_coverage = float(size)/float(genome_sizes[gi])
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
        alignments[(name, gname, gi)] = (gname, float(value)/float(genomes[gname]))

    for track, proportion in alignments.iteritems():
        outf.write("%s\t%s\t%s\t%f\n" % (track[0], track[2], proportion[0], proportion[1]))
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

    R('''data <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
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
