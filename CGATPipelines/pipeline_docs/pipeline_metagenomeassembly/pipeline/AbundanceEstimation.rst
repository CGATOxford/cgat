.. _Abundances:


===============================
Relative taxonomic abundances
===============================


Relative abundance
===================


The metagenome assembly pipeline also attempts to estimate the relative abundance of 
taxonomic groups within each sample. At present it supports the use of metaphlan.
See the metaphlan paper_ for more information. You must remember that these are
onlyt the relative abundances based on reads that could be assigned to a specific
clade.

The method uses a balst/bowtie alignment to a known set of clade-specific markers
that allows us to estimate the relative abundance of taxa at different resolutions i.e.
kingdom, phylum, class, order, family, genus, species.

Below is a per sample description of the relative abundance of each taxon at each
resolution (% of total reads that could be mapped to the respective level).

Note that only the relative abundances of those groups that are present at > 1% are shown.


.. report:: RelativeAbundance.RelativeAbundance
   :render: table

   Relative abundances


.. report:: RelativeAbundance.RelativeAbundance
   :render: stacked-bar-plot
   :mpl-rc: legend.fontsize=10
   :legend-location: lower-right

   Relative abundances


The number of reads that were used for assigning taxonomic groups is described below. It is likely that many reads
will not contribute to the distribution as many will not map to clade-specific markers in the database. In addition
this description if of reads that are assigned at the highest taxonomic level i.e. phylum and 
there are likely many fewer reads that are assigned at the level of the species.


.. report:: RelativeAbundance.ContributingReads
   :render: table

   Proportion of reads that contribute to assignment of taxonomic level


.. report:: RelativeAbundance.ContributingReads
   :render: interleaved-bar-plot

   Proportion of reads that contribute to assignment of taxonomic level




Total number of taxa
====================

Below are plots displaying the total number of taxa that were identified in each sample.

.. report:: RelativeAbundance.TotalSpecies
   :render: table

   Total taxa


.. report:: RelativeAbundance.TotalSpecies
   :render: interleaved-bar-plot
   :mpl-rc: legend.fontsize=10
   :legend-location: lower-right

   Total taxa


.. _paper: http://www.ncbi.nlm.nih.gov/pubmed/22688413
