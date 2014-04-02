.. _LcaAbundances:


================================================================
Lowest Common Ancestor (LCA) - Relative abundance estimations
================================================================


Relative abundance
====================

The metagenome communities pipeline attempts to estimate the relative abundance of 
taxonomic groups within each sample using the lowest common ancestor approach as 
implemented in MEGAN. Briefly, each sequence is aligned to a database (protein or nucleotide)
using BLAST. A taxonomic assignment is made for each read based on the lowest common ancestor 
of each of the reads alignments i.e. the sequence will be assigned at the species level if
all of its alignments come from the same species. However, where there is high sequence similarity between
species, the sequence will be assigned at a higher node in the NCBI taxonomy - a node that
represents the lowest common ancestor of each of the alignments. In contrast to metaphlan,
LCA attempts to assign every read to a taxonomy. Note that only the relative abundances of those 
groups that are present at > 1% are shown.


.. report:: Lca.LcaRelativeAbundance
   :render: table

   Relative abundances estimated by LCA


.. report:: Lca.LcaRelativeAbundance
   :render: bar-plot
   :groupby: slice

   Relative abundances as estimated by LCA


Contributing reads
====================

The number of reads that were used for estimating taxonomic abundances is described below. It is likely that many reads
will not contribute to the distribution as many will not map sequences with high enough similarity.
The number of reads contributing to abundance estimations using LCA is greater than metaphlan due to the
reduced database that is used to estimate relative abundances using metaphlan.


.. report:: Lca.LcaContributingReads
   :render: table

   Total number and proportion of reads that contribute to assignment of taxonomic level


Total number of taxa
=======================

Below are plots displaying the total number of taxa that were identified in each sample. This provides
some information about the diversity of the sample.

.. report:: Lca.LcaTotalTaxa
   :render: table

   Total taxa


.. report:: Lca.LcaTotalTaxa
   :render: bar-plot

   Total taxa


.. _paper: http://www.ncbi.nlm.nih.gov/pubmed/22688413
