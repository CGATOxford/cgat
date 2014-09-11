==========
Complexity
==========

This section looks at library duplication levels and estimated complexity using metrics from the Picard function MarkDuplicates.

Please note that complexity plots are only produced for paired-end data (that contains non-optical duplicates). If in doubt please check for the presence of a histogram in the ".....picard_duplication_metrics" file(s) in the bam directory. 

Picard MarkDuplicates metrics for each BAM file for each :term:`track`. See the 
`Picard metrics <http://picard.sourceforge.net/picard-metric-definitions.shtml#AlignmentSummaryMetrics>`_
for a definition of the field contents.
(note: PF=pass filter, reads that pass the vendor's filter criteria).

.. report:: Mapping.PicardDuplicationSummary
   :render: table

Estimated return from extra sequencing.
=======================================

Data table

.. report:: Mapping.DuplicationMetricsTable
   :render: table

All histograms. In the plot the X axis represent the factor by which the present amount of sequencing is multiplied by: the current depth represents a coverage multiple of 1. The y axis represents the actual additional coverage predicted.

.. report:: Mapping.DuplicationMetrics
   :render: line-plot
   :as-lines:
   :yrange: 0,
   :xrange: 0,20
   :ytitle: additional actual coverage

   Individual plots

.. report:: Mapping.DuplicationMetrics
   :render: line-plot
   :as-lines:
   :xrange: 0, 20
   :yrange: 0,
   :groupby: track
   :ytitle: additional actual coverage
   :layout: column-3
   :width: 200

   Grouped plots

