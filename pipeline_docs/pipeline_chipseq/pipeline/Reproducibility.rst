===============
Reproducibility
===============

Overlap in bases
================

.. report:: Overlaps.OverlapsBasesNormalized
   :render: matrix-plot
   :transform-matrix: square,correspondence-analysis
   :tracks: replicates
   
   Percent overlap of bases between sets. The matrix plot displays for each pair
   of set the quotient of bases that are commen between both sets and bases
   that are present in any of the sets (intersection / union). 

.. report:: Overlaps.OverlapsBasesNormalized
   :render: matrix
   :transform-matrix: square,correspondence-analysis
   :tracks: replicates
   
   Percent overlap of bases between sets. The matrix plot displays for each pair
   of set the quotient of bases that are commen between both sets and bases
   that are present in any of the sets (intersection / union). 

Overlap in intervals
====================

.. report:: Overlaps.OverlapsExonsPercent
   :render: matrix-plot
   :tracks: replicates
   :transform-matrix: square,sort,symmetric-max

   Percent overlap between sets. The matrix plot displays for each pair
   of set the quotient of intervals that are common between both sets.
   The matrix has been symmetrized

.. report:: Overlaps.OverlapsExonsPercent
   :render: matrix
   :tracks: replicates
   :transform-matrix: square,sort

   Percent overlap between sets. The matrix plot displays for each pair
   of set the quotient of intervals that are common between both sets.

Peak heights
============

The following plot shows the maximum (peak) interval value for each set.
The peak value is the maximum number of reads falling into the
bins that constitute an interval. The peak is the position with the maximum
number of reads.

.. report:: ChipseqReport.IntervalPeakValues
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total,reverse-cumulative
   :as-lines:
   :tracks: replicates

   Distribution of the number of reads at the peak within an interval.
   The distribution list the proportion of intervals of a certain peak
   value or more.
