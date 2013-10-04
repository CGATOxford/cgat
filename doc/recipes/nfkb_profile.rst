
What is the binding profile of NFKB across gene models?
========================================================

After processing RNA-seq data through alignment, gene/transcript abundance estimation and differential
expression analysis, we are left with an unannotated list of differentially expressed genes. We may want
to know whether there is evidence to suggest that these genes are regulated by a transcription factor
of interest. We can answer this using ChIP-seq data that we ourselves have generated or by using
public resources such as ENCODE.

For example, we have carried out an RNA-seq experiment in lymphoblastoid cell lines (LCLs) looking at
the effect of TNF-a stimulation on gene expression. Using one of the many tools for conducting
differential expression analysis we have arrived at a set of 133 genes that are up-regulated when
LCLs are stimulated with TNF-a. 

We know that the main transcription factor that drives expression of inflammatory genes using other
immune stimulators is NFKB. We would therefore like to answer the question:

Is there evidence to support a role for NFKB in the regulation of genes regulated by TNF-a in LCLs?

ENCODE have produced many ChIP-seq data sets and by a stroke of luck they have NFKB ChIP-seq data in 
TNF-a stimulated LCLs. In an exploratory phase of the analysis, we would like to see what the profile
of NFKB binding is across genes i.e does it bind predominantly at the TSS, exons or 3' UTR. We 
can do this fairly easily with a few files and a few commands.

The input files that we require are:

* A gtf file of protein coding gene transcripts. 

* A bam file of NFKB ChIP-seq reads (NFKB.bam).

Genesets can be easily downloaded from ENSEMBL or UCSC. For example we can download the ENSEMBL reference geneset
by typing::

    wget ftp://ftp.ensembl.org/pub/release-73/gtf/homo_sapiens/Homo_sapiens.GRCh37.73.gtf.gz -o logfile

We can then take all protein coding genes from thsi geneset with an awk statement::

    zcat Homo_sapiens.GRCh37.73.gtf.gz  | awk '$2=="protein_coding"' | gzip > coding_geneset.gtf.gz

Using the CGAT tool bam2geneprofile we can assess the binding profile of NFKB across gene models::

    cgat bam2geneprofile --bamfile=NFKB.bam 
                         --gtffile=coding_geneset.gtf.gz 
                         --method=geneprofile 
                         --profile_normalization=counts
                         --output-filename-pattern=nfkb_profile_%s


This statement will produce a matrix as an output file named "nfkb_profile.geneprofile.matrix.tsv.gz" with the following format::

   +---+--------+----------+--------------+
   |bin|region  |region_bin|counts        |
   +---+--------+----------+--------------+
   |0  |upstream|0         |0.22691292876 |
   +---+--------+----------+--------------+
   |1  |upstream|1         |0.224274406332|
   +---+--------+----------+--------------+
   |2  |upstream|2         |0.221635883905|
   +---+--------+----------+--------------+
   |3  |upstream|3         |0.192612137203|
   +---+--------+----------+--------------+
   |4  |upstream|4         |0.221635883905|
   +---+--------+----------+--------------+
   |5  |upstream|5         |0.213720316623|
   +---+--------+----------+--------------+
   |6  |upstream|6         |0.213720316623|
   +---+--------+----------+--------------+
   |7  |upstream|7         |0.200527704485|
   +---+--------+----------+--------------+
   |8  |upstream|8         |0.20580474934 |
   +---+--------+----------+--------------+
 

This data is amenable to further manipulation and visualisation. For example, we can use R to produce a profile plot
over the gene model. Start R and type::

   R version 2.15.2 (2012-10-26) -- "Trick or Treat"
   Copyright (C) 2012 The R Foundation for Statistical Computing
   ISBN 3-900051-07-0
   Platform: x86_64-unknown-linux-gnu (64-bit)

   R is free software and comes with ABSOLUTELY NO WARRANTY.
   You are welcome to redistribute it under certain conditions.
   Type 'license()' or 'licence()' for distribution details.

   R is a collaborative project with many contributors.
   Type 'contributors()' for more information and
   'citation()' on how to cite R or R packages in publications.

   Type 'demo()' for some demos, 'help()' for on-line help, or
   'help.start()' for an HTML browser interface to help.
   Type 'q()' to quit R.

   > profile <- read.csv("nfkb_profile.geneprofile.matrix.tsv.gz", header = T, stringsAsFactors = F, sep = "\t")
 
   > plot(profile$bin, profile$counts, cex = 0)   

   > lines(profile$bin, profile$counts, col = "blue")

   > abline(v = c(1000, 2000))

   > mtext("upstream", adj = 0.1)
    
   > mtext("exons", adj = 0.5)

   > mtext("downstream", adj = 0.9)


This set of commands will produce the figure shown.

.. image:: ../plots/nfkb_profile.pdf 


This plot displays the predominance of NFKB binding at transcription start sites of protein coding genes. 

While NFKB binds to the TSSs of protein coding genes, it also binds to many intergenic regions of the genome. In addition
to meta-gene profiles we may also want to know the chromatin state at which NFKB binding occurs. For example, we can
integrate additional histone modification ChIP-seq data from the ENCODE project. H3K4me3 and H3K4me1 mark promoters and 
enhancers, respectively. We would like to visualise the profile of these marks at all the genomic locations of
NFKB binding.

The input files for this analysis are:

* bed intervals describing NFKB peaks (e.g. from a MACS peak calling analysis)

* H3K4me1 and H3K4me3 alignment files in .bam format.







