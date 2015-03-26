'''
classes and utility functions for pipeline_metagenomecommunities.py

Different assembly tools will use different inputs. Some can take
fasta files whereas others will take fastq and in either case can be
paired-end (in the same or different files) or single end

'''

import os
import collections
import shutil

import CGAT.IOTools as IOTools
import CGAT.Pipeline as P
import CGAT.Experiment as E
import CGAT.FastaIterator as FastaIterator
import CGAT.Fastq as Fastq
import glob
import CGAT.Metaphlan as Metaphlan
import numpy as np
from rpy2.robjects import r as R
import CGAT.GTF as GTF


def normaliseKraken(infile, outfile):
    '''
    normalise kraken counts by nreads/million mapped
    '''
    inf = IOTools.openFile(infile)
    header = inf.readline().replace("rel_abundance", "rpm")
    mapped = 0

    # will have to iterate over the file twice 
    for line in inf.readlines():
        data = line[:-1].split("\t")
        count = int(data[-1])
        mapped += count
    inf.close()
    
    inf = IOTools.openFile(infile)
    inf.readline()
    outf = IOTools.openFile(outfile, "w")
    outf.write(header)
    for line in inf.readlines():
        data = line[:-1].split("\t")
        count = int(data[-1])/(float(mapped)/1000000)
        outf.write("\t".join(map(str, data[:-1] + [count])) + "\n")
    outf.close()

###################################################################
###################################################################
###################################################################
    
def countContributingReads(infile, outfile):
    '''
    count number of reads with a taxnomic assignment
    '''
    levels = [
        "phylum", "class", "order", "family", "genus", "species"]
    result = collections.OrderedDict()
    for level in levels:
        result[level] = 0
    inf = IOTools.openFile(infile)
    header = inf.readline().split("\t")

    # column indices
    indices = [3,5,7,9,11,13] 
    total = 0
    for line in inf.readlines():
        total += 1
        data = line[:-1].split("\t")
        phylum, _class, order, family, genus, species = [data[i] for i in indices]
        if phylum != "NA":
            result["phylum"] += 1
        if _class != "NA":
            result["class"] += 1
        if order != "NA":
            result["order"] += 1
        if family != "NA":
            result["family"] += 1
        if genus != "NA":
            result["genus"] += 1
        if species != "NA":
            result["species"] += 1
    outf = open(outfile, "w")
    outf.write("level\tn_reads\tpct_reads\n")
    for level, count in result.iteritems():
        nreads, prop = count, float(count) / total
        outf.write("\t".join([level, str(nreads), str(prop*100)]) + "\n")
    outf.close()

###################################################################
###################################################################
###################################################################

def plotMDS(infile, outfile):
    '''
    perform multidimensional scaling of normalised
    counts
    '''
    outname_matrix = P.snip(outfile, ".pdf") + ".tsv"
    R('''library(gtools)''')
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''rownames(dat) <- dat$taxa
         dat <- dat[,1:ncol(dat)-1]
         dat <- dat[, mixedsort(colnames(dat))]
         conds <- unlist(strsplit(colnames(dat), ".R[0-9]"))[seq(1, ncol(dat)*2, 2)]
         conds <- unlist(strsplit(conds, ".", fixed = T))[seq(2, length(conds)*2, 2)]
         dat <- as.matrix(t(dat))
         dist <- dist(dat)
         ord1 <- cmdscale(dist)
         ord2 <- as.data.frame(ord1)
         ord2$cond <- conds
         ggplot(ord2, aes(x = V1, y = V2, colour = cond)) + geom_point(size = 3) + scale_colour_manual(values = c("brown", "darkGreen", "slateGrey", "darkBlue"))
         ggsave("%s")''' % outfile)

###################################################################
###################################################################
###################################################################

def testDistSignificance(infile, outfile):
    '''
    use permanova to test significance of differences
    in distance metrics
    '''
    R('''library(gtools)''')
    R('''library(ggplot2)''')
    R('''library(vegan)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''rownames(dat) <- dat$taxa
         dat <- dat[,1:ncol(dat)-1]
         dat <- dat[, mixedsort(colnames(dat))]
         conds <- unlist(strsplit(colnames(dat), ".R[0-9]"))[seq(1, ncol(dat)*2, 2)]
         conds <- unlist(strsplit(conds, ".", fixed = T))[seq(2, length(conds)*2, 2)]
         d <- dist(t(dat))
         paov <- adonis(d~conds)
         write.table(paov$aov.tab, file = "%s")''' % outfile)
         
###################################################################
###################################################################
###################################################################

def barplotAbundance(infile, outfile):
    '''
    barplot normalised counts abundance
    '''
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''rownames(dat) <- dat$taxa''')
    R('''dat <- dat[,2:ncol(dat)]''')
    R('''library(ggplot2)''')
    R('''library(reshape)''')
    R('''library(gtools)''')

    R('''sums <- apply(dat, 2, sum)''')
    R('''dat2 <- data.frame(t(apply(dat, 1, function (x) x/sums)))''')
    R('''dat2$taxa <- rownames(dat)''')
    R('''dat2 <- dat2[, mixedsort(colnames(dat2))]''')
    R('''dat2 <- melt(dat2)''')
    R('''ggplot(dat2, aes(x = variable, y = value, fill = taxa, stat = "identity")) + geom_bar() + opts(axis.text.x=theme_text(angle=90))''')
    R('''ggsave("%s")''' % outfile)

###################################################################
###################################################################
###################################################################

def MAPlot(infile,
           threshold_stat,
           p_threshold, 
           fc_threshold,
           outfile):
    '''
    MA plot the results
    '''
    
    if threshold_stat== "p":
        p="P.Value"
    elif threshold_stat == "padj":
        p="adj.P.Val"
    else:
        p="adj.P.Val"

    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''dat$sig <- ifelse(dat$%s < %f & abs(dat$logFC) > %f, 1, 0)''' % (p, p_threshold, fc_threshold))
    R('''ggplot(dat, aes(x = AveExpr, y = logFC, colour = factor(sig))) + geom_point(alpha = 0.5) + scale_colour_manual(values = c("black", "blue"))''')
    R('''ggsave("%s")''' % outfile)
   
###################################################################
###################################################################
###################################################################

def plotHeatmap(results, 
                norm_matrix, 
                threshold_stat,
                p_threshold, 
                fc_threshold,
                outfile):
    '''
    plot heatmap of differentially abundant genes
    '''
    if threshold_stat== "p":
        p="P.Value"
    elif threshold_stat == "padj":
        p="adj.P.Val"
    else:
        p="adj.P.Val"

    temp = P.getTempFilename(".")
    R('''library(gplots)''')
    R('''library(gtools)''')
    E.info("reading data")
    R('''mat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % norm_matrix)
    R('''rownames(mat) <- mat$taxa
         mat <- as.matrix(mat[,1:ncol(mat)-1])''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % results)
    E.info("data loaded")

    R('''diff.genes <- unique(dat$taxa[dat$%s < %f & abs(dat$logFC) > %f])''' % (p, p_threshold, fc_threshold))

    ##############################
    # this is a hack
    # to avoid errors when 
    # a single differential
    # abundant feature is found
    ##############################
    R('''write.table(diff.genes,file = "%s", row.names = F, sep = "\t")''' % temp)
    
    tmp = open(temp)
    tmp.readline()
    if len(tmp.readlines()) == 1: 
        P.touch(outfile)
    else:
        R('''mat <- mat[as.character(diff.genes), ]
             samples <- colnames(mat)
             mat <- as.data.frame(t(apply(mat, 1, scale)))
             colnames(mat) <- samples
         mat <- mat[, mixedsort(colnames(mat))]
         colours = colorRampPalette(c("blue", "white", "red"))(75)
         png("%s", height = 1000, width = 500)
         heatmap.2(as.matrix(mat), 
                   trace = "none", 
                   scale = "none", 
                   col = colours,
                   Colv = F,
                   dendrogram = "row",
                   margins = c(15, 15))
             dev.off()''' % outfile)

    os.unlink(temp)

###################################################################
###################################################################
###################################################################

def annotate(infile, annotation_file, outfile):
    '''
    annotate infile with annotations from 
    annotation gtf file
    '''
    inf = open(infile)
    header = inf.readline()
    include = set()

    E.info("reading genes to keep")
    for line in inf.readlines():
        data = line[:-1].split("\t")
        gene_id = data[8].strip('"')
        include.add(gene_id)

    E.info("reading annotations file")
    annotations = {}
    for gtf in GTF.iterator(IOTools.openFile(annotation_file)):
        if gtf.gene_id in include:
            annotations[gtf.gene_id] = [gtf.gene_name, gtf.species, gtf.description]

    inf = open(infile)
    header = inf.readline()
    
    E.info("writing results with annotations")
    outf = open(outfile, "w")
    outf.write(header.strip("\n") + "\tgene_name\tspecies_centroid\tdescription\n")
    for line in inf.readlines():
        data = line[:-1].split("\t")
        gene_id = data[8].strip('"')
        try:
            outf.write("\t".join(data + annotations[gene_id]) + "\n")
        except KeyError:
            outf.write("\t".join(data + ["NA", "NA", "NA"]) + "\n")
    outf.close()


###################################################################
###################################################################
###################################################################

def rarefactionCurve(infile, 
                     outfile,
                     rdir,
                     f = 1000000,
                     step = 1000000):
    '''
    perform and plot rarefaction curves 
    on species richness - the max rarefaction sample
    will be the minimum count across datasets
    '''
    R('''source("%s/metagenomic_diversity.R")''' % rdir)
    R('''library(gtools)''')
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''rownames(dat) <- dat$taxa''')
    R('''dat <- dat[,2:ncol(dat)]''')
    R('''tdat <- data.frame(t(dat))''')
    R('''sums <- apply(tdat, 1, sum)''')
    R('''min.count <- min(sums)''')
    R('''tdat <- tdat[mixedsort(rownames(tdat)),]''')
    R('''conds <- unlist(strsplit(rownames(tdat), ".R[0-9]"))[seq(1, nrow(tdat)*2, 2)]''')
    R('''conds <- unlist(strsplit(conds, ".", fixed = T))[seq(2, length(conds)*2, 2)]''')
    R('''rf <- rarefaction(tdat, from = %i, to = min.count, step = %i, groups = conds)''' % (f, step))
    R('''plotRarefaction(rf)''')
    R('''ggsave("%s")''' % outfile)


###################################################################
###################################################################
###################################################################

def testRichness(infile, 
                 outfile,
                 rdir,
                 sample):
    '''
    test significance of richness between
    groups using kruskal wallis test
    '''
    R('''source("%s/metagenomic_diversity.R")''' % rdir)
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)''' % infile)
    R('''tdat <- data.frame(t(dat))''')
    R('''conds <- unlist(strsplit(rownames(tdat), ".R[0-9]"))[seq(1, nrow(tdat)*2, 2)]''')
    R('''conds <- unlist(strsplit(conds, ".", fixed = T))[seq(2, length(conds)*2, 2)]''')
    R('''k <- richness.test(tdat, %i, groups = conds)''' % sample)
    R('''write.table(data.frame(p = k$p.value, stat = k$statistic), file = "%s", sep = "\t", row.names = F)''' % outfile)

###################################################################
###################################################################
###################################################################

def barplotDiversity(infile,
                     outfile,
                     rdir,
                     ind = "shannon"):
    '''
    barplot diversity
    '''
    R('''source("%s/metagenomic_diversity.R")''' % rdir)
    R('''library(gtools)''')
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''rownames(dat) <- dat$taxa''')
    R('''dat <- dat[,2:ncol(dat)]''')
    R('''tdat <- data.frame(t(dat))''')
    R('''tdat <- tdat[mixedsort(rownames(tdat)),]''')
    R('''conds <- unlist(strsplit(rownames(tdat), ".R[0-9]"))[seq(1, nrow(tdat)*2, 2)]''')
    R('''conds <- unlist(strsplit(conds, ".", fixed = T))[seq(2, length(conds)*2, 2)]''')
    R('''plotDiversity(tdat, index = "%s", groups = conds)''' % ind)
    R('''ggsave("%s")''' % outfile)

###################################################################
###################################################################
###################################################################

def testDiversity(infile,
                     outfile,
                     rdir,
                     ind = "shannon"):
    '''
    barplot diversity
    '''
    R('''source("%s/metagenomic_diversity.R")''' % rdir)
    R('''library(gtools)''')
    R('''library(ggplot2)''')
    R('''dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")''' % infile)
    R('''rownames(dat) <- dat$taxa''')
    R('''dat <- dat[,2:ncol(dat)]''')
    R('''tdat <- data.frame(t(dat))''')
    R('''tdat <- tdat[mixedsort(rownames(tdat)),]''')
    R('''conds <- unlist(strsplit(rownames(tdat), ".R[0-9]"))[seq(1, nrow(tdat)*2, 2)]''')
    R('''conds <- unlist(strsplit(conds, ".", fixed = T))[seq(2, length(conds)*2, 2)]''')
    R('''k <- div.test(tdat, groups = factor(conds))''')
    R('''write.table(data.frame(p = k$p.value, stat = k$statistic), file = "%s", sep = "\t", row.names = F)''' % outfile)


