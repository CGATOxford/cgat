============
CpG analysis
============

.. _CpGEnrichment:

CpG Enrichment
==============

The following plot shows the distribution of percentage CpG in
regions covered by reads in the individual data tracks. These plots
can give and indication if there is any enrichment for CpG in
genomic fragments pulled down.

.. report:: CpGCoverage.CpGDistributionInTags
   :render: line-plot
   :xrange: 0,0.2
   :as-lines:
   :split-at: 10
   :split-always: genomic
   :xtitle: CpG content / %
   :groupby: all
   :layout: column-2

   Distribution of CpG in regions covered by reads

.. report:: CpGCoverage.CpGDistributionInTags
   :render: line-plot
   :transform: aggregate
   :tf-aggregate: normalized-total,cumulative
   :xrange: 0,0.2
   :as-lines:
   :split-at: 10
   :split-always: genomic
   :xtitle: CpG content / %
   :groupby: all
   :layout: column-2

   Cumulative distribution of CpG in regions covered by reads

.. _CpGCoverage:

CpG coverage
============

This section displays the coverage of CpG dinucleotides by tags. These
plots help checking if the sequencing depth has been sufficient to
permit differential methylation analysis. Note that this is a
nucleotide-level analysis, while in an MeDIP-seq a sequence fragment
will span multiple CpG sites and usually clusters of CpG are analyzed.

.. report:: CpGCoverage.CpGCoverageByReads
   :render: line-plot
   :as-lines:
   :transform: aggregate
   :xrange: 0,20
   :yrange: 0,1
   :tf-aggregate: normalized-total,reverse-cumulative
   :xtitle: read depth
   :split-at: 10
   :layout: column-2
   :groupby: all

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


