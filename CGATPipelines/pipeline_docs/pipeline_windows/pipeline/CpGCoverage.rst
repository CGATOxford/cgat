===========
CpGCoverage
===========

This section displays the coverage of CpG dinucleotides by reads
within the read data.

.. report:: CpGCoverage.CpGCoverage
   :render: line-plot
   :transform: aggregate
   :logscale: y
   :xrange: 0,50
   :tf-aggregate: cumulative

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
   :render: table
   :tracks: total,protein_coding,repeats
   :pie-first-is-total: other

   Number of CpG dinucleotides in protein coding regions and repeats.

.. report:: CpGCoverage.CpGContext
   :render: table
   :force:

   Number of CpG dinucleotides in different genomic regions.

