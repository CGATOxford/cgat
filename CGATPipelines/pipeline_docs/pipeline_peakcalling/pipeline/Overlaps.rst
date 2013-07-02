 .. _OverlapBetweenSets:

*****************************
Overlap between sets
*****************************

This page presents the overlap between sets.

Plots of pairwise overlap between sets
======================================

The following table lists the percent of intervals in a track that overlap with exons in another track.

.. report:: Overlaps.OverlapsExonsPercent
   :render: matrix-plot
   :transform-matrix: square,correspondence-analysis

   Percent overlap of intervals between sets. The matrix plot displays for each set in a row 
   the percentage of intervals in that set that overlap with intervals
   in the set in the column.  Note that this plot need not be symmetric if the
   numbers of intervals in the two sets are very different.


.. report:: Overlaps.OverlapsBasesNormalized
   :render: matrix-plot
   :transform-matrix: square,correspondence-analysis

   Percent overlap of bases between sets. The matrix plot displays for each pair
   of set the quotient of bases that are commen between both sets and bases
   that are present in any of the sets (intersection / union).

Tables of overlap between sets
==============================

The following tables list the number of exons in a track 
that overlap with exons in another track.

.. report:: Overlaps.OverlapsExonsPercent
   :render: matrix
   :transform-matrix: square,correspondence-analysis

   Percent overlap sets counted by exons

.. report:: Overlaps.OverlapsExonsCounts
   :render: matrix
   :transform-matrix: square,correspondence-analysis

   Overlap between sets counted by exons

.. report:: Overlaps.OverlapsBasesPercent
   :render: matrix
   :transform-matrix: square,correspondence-analysis

   Percent base overlap between sets


Reproducibility ROC curves
==========================

The following ROC curves correlate various peak measures with 
reproducibility. The ROC curves are computed by sorting all 
intervals by value. True positives are those intervals that 
are found in all replicates.

.. report:: Overlaps.OverlapROC
   :render: line-plot
   :tracks: merged
   :layout: column-3
   :as-lines: 
   :width: 300
   :groupby: track

   ROC curves examining the reproducibility between sets.


