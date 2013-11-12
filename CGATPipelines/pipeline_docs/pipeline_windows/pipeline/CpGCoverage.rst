============
CpG analysis
============

.. _CpGEnrichment:

CpG Enrichment
==============

The following plot shows the distribution of percentage CpG in
regions covered by reads in the individual data tracks.

.. report:: CpGCoverage.CpGDistribution
   :render: line-plot
   :transform: histogram
   :xrange: 0,0.2
   :as-lines:
   :split-at: 10

   Distribution of CpG in regions covered by reads

.. _CpGCoverage:

CpG coverage
============

This section displays the coverage of CpG dinucleotides by reads
within the read data.

.. report:: CpGCoverage.CpGCoverage
   :render: line-plot
   :transform: aggregate
   :xrange: 0,20
   :yrange: 0,1
   :tf-aggregate: normalized-total,cumulative
   :xtitle: read depth

   Cumulative plot of coverage of CpG dinucleotides with 
   reads.

The following table shows the number of CpG in total,
within protein coding sequence and within repeats.

.. report:: CpGCoverage.CpGContext
   :render: table
   :tracks: total,protein_coding,repeats
   :force:

   Number of CpG dinucleotides in protein coding regions and repeats.

.. report:: CpGCoverage.CpGContext
   :render: pie-plot
   :tracks: total,protein_coding,repeats
   :pie-first-is-total: other

   Number of CpG dinucleotides in protein coding regions and repeats.

.. report:: CpGCoverage.CpGContext
   :render: table
   :force:
   :transpose:

   Number of CpG dinucleotides in different genomic regions.


