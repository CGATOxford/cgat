===========================
Using CGAT tools - Recipes
===========================

In this section you will find representative examples for using tools
developed in CGAT. The recipes presented aim to provide intuitive
real-life examples of CGAT script use for the analysis of genomic
datasets. If there is a tool in the CGAT collection for which you
would like a use case then please post a request on the `CGAT users
group`_ website.

The recipes are implemented as ipython_ notebooks.

..
   .. toctree::
      :maxdepth: 2

      recipes/gat
      recipes/metagenome_contigs_kmers

:download:`Recipe02 Plotting read-density in Intervals <recipes/Recipe02-BAMCoverageHistograms.html>`
	  Illustrate how to plot read density in a large number of
	  genomic intervals, for example in the ChIP-seq experiment.

:download:`Recipe03 Calculating CpG content in promotors <recipes/Recipe03-CpGInLncRNA.html>`
	  Illustrate how to caluclate CpG content from
	  a gene list of genes.


:download:`Recipe04 Plotting meta-gene profiles <recipes/Recipe04-MetaGeneProfiles.html>`
	  Illustrate how to compute meta-gene profiles.

:download:`Recipe05 Stripping Bam Files <recipes/Recipe05-StrippingBamFiles.html>`
	  Illustrate how to remove sequence and/or
	  quality information from BAM files and
	  how to put it back.

.. _CGAT users group: https://groups.google.com/forum/?fromgroups#!forum/cgat-user-group

Tools relevant for RNA-seq
==========================

Below is a list of tools that are of interest for RNA-seq
analysis.

:doc:`scripts/gtf2table`

   This script takes a :term:`GTF` formatted file and outputs for each
   gene or transcript one or more annotations that are computed by
   integrating the gene/transcripts with additional data. For example,
   by giving it a :term:`BAM` formatted file, it will output the
   number of reads overlapping with a particular transcript. Adding
   another :term:`GTF` file with a reference gene set will annotate
   each transcript according to the overlap and classify it as
   fragment, novel, etc. The tabular output from one or more
   :term:`gtf2table` runs together with additional annotation data can
   be used to extract transcripts/genes of interest, for example, to
   select genes of high G+C content that have anti-sense expression
   and are novel.

:doc:`scripts/counts2counts`

   This script manipulates a counts table such as one created by
   :doc:`scripts/gtf2table` or :term:`featureCounts`. It can
   normalize, filter, or compute random permutations for post-hoc
   power analysis.

:doc:`scripts/counts2table` 

   (an older version of this is :doc:`scripts/runExpression`). This script
   takes a table with count data and applies statistical tests to
   detect differences between sample groups. It wraps several methods
   of relevance in RNA-seq analysis such as :term:`DEseq` and
   :term:`EdgeR`.

:doc:`scripts/gtf2gtf`

   This script allows manipulation of GTF files such as sorting,
   filtering, renaming but also manipulation of gene models such as
   combining transcripts into genes, etc.

:doc:`scripts/diff_gtf`

   This scripts compares two gene sets and outputs the number of
   shared and unique genes, exons and bases. These are standard
   metrics used in gene-prediction.

:doc:`scripts/gtfs2tsv`

   This scripts compares two gene sets and outputs lists of 
   shared and unique genes.

:doc:`scripts/gtf2tsv`
   
   This script converts a :term:`GTF` formatted file into tabular
   format including the transcript/gene attributes. This is useful for
   uploading the geneset into a database.

:doc:`scripts/expression2distance`

   This script generates a distance matrix for time-series data.

:doc:`scripts/expression2expression`

   This script normalizes and transforms RNA-seq time-series data.

All these scripts work from and output standard genomic file formats
and are thus easily integrated with other tools such as
:term:`bedtools`.

There are some additional scripts for gene-set manipulation that might
be of interest:

:doc:`scripts/gtfs2graph`

:doc:`scripts/diff_transcript_sets`


