Bed formatted intervals
=======================

This page summarizes some basic counts in the BED formatted
files in the annotation pipeline. Bed files typically contain
genomic intervals.

.. report:: AnnotationReport.BedSummaryIntervals
   :render: interleaved-bar-plot
   :logscale: y

   Number of intervals in various BED files. Note that
   the Y-axis is log-scaled.


.. report:: AnnotationReport.BedSummaryBases
   :render: interleaved-bar-plot
   :logscale: y

   Number of bases covered by intervals in various
   BED files. Note that this number can be larger
   than the genome size as intervals can overlap.
   Note that the y-axis is log-scaled.


Per contig summaries
---------------------

.. report:: AnnotationReport.BedSummaryIntervalsPerContig
   :render: matrix-plot
   :logscale: z
   :mpl-figure: figsize=(20,20)
   :transform-matrix: sort
   :transform: pivot
   :pivot-index: track
   :pivot-value: nintervals
   :pivot-column: contig

   Number of intervals in various annotation files

.. report:: AnnotationReport.BedSummaryBasesPerContig
   :render: matrix-plot
   :logscale: z
   :mpl-figure: figsize=(20,20)
   :transform-matrix: sort
   :transform: pivot
   :pivot-index: track
   :pivot-value: nbases
   :pivot-column: contig

   Number of nucleotides covered in various annotation
   files
