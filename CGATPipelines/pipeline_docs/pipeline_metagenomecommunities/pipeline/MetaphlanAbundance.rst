.. _Abundances:


===========================================
Metaphlan - Relative abundance estimations
===========================================


Species Abundance distribution
===============================

Below we describe the relative abundance distribution of species identified in the samples by Metaphlan.

.. report:: Metaphlan.SpeciesAbundanceDistribution
   :render: bar-plot
   
   Relative abundance distribution of species



Relative abundance
====================

The metagenome assembly pipeline also attempts to estimate the relative abundance of 
taxonomic groups within each sample. At present it supports the use of metaphlan.
See the metaphlan paper_ for more information. You must remember that these are
onlyt the relative abundances based on reads that could be assigned to a specific
clade.

The method uses a blast/bowtie alignment to a known set of clade-specific markers
that allows us to estimate the relative abundance of taxa at different resolutions i.e.
kingdom, phylum, class, order, family, genus, species.

Below is a per sample description of the relative abundance of each taxon at each
resolution (% of total reads that could be mapped to the respective level).

Note that only the relative abundances of those groups that are present at > 1% are shown.


.. report:: Metaphlan.RelativeAbundance
   :render: table

   Relative abundances


.. report:: Metaphlan.RelativeAbundance
   :render: bar-plot
   :groupby: slice

   Relative abundances as estimated by metaphlan


Contributing reads
====================

The number of reads that were used for estimating taxonomic abundances is described below. It is likely that many reads
will not contribute to the distribution as many will not map to clade-specific markers in the database. In addition
this description is of reads that are assigned at the highest taxonomic level i.e. phylum and 
there are likely many fewer reads that are assigned at the level of the species. In addition to this, relative abundance
estimations by metaphlan are not based upon all reads that are assigned i.e. low abundance taxa are not represented in the
relative abundance estimations. The nu,ber of contributing reads is therefore an overestimation.


.. report:: Metaphlan.ContributingReads
   :render: table

   Total number and proportion of reads that contribute to assignment of taxonomic level


Total number of taxa
=======================

Below are plots displaying the total number of taxa that were identified in each sample. This provides
some information about the diversity of the sample.

.. report:: Metaphlan.TotalTaxa
   :render: table

   Total taxa


.. report:: Metaphlan.TotalTaxa
   :render: bar-plot
   :groupby: slice

   Total taxa


.. _paper: http://www.ncbi.nlm.nih.gov/pubmed/22688413
