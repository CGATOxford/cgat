=========================
Transcription Start Sites
=========================

This part of the pipeline examines the distances of intervals
to transcription start sites. Transcription start sites are
usually refered to distances to the :term:`reference gene set`.

Overlap with TSS
----------------

Number of TSS start sites that intervals overlap.

.. report:: TSS.TSSOverlap
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total

   Number of TSS start sites that intervals overlap.

Distance to closest TSS
-----------------------

The following plots show the distance of each 
interval to the closest TSS.

.. report:: TSS.TSSClosest
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,100000
   :yrange: 0,1
   :tf-aggregate: normalized-total,cumulative
   :as-lines:

   Histogram of distances of interval to TSS

.. report:: TSS.TSSClosest
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,5000
   :tf-range: 0,100000,100
   :yrange: 0,1
   :tf-aggregate: normalized-total,cumulative
   :as-lines:

   Histogram of distances of interval to TSS

Closest upstream TSS
--------------------

The following plots show the distance of each 
interval to the closest TSS that is upstream
of the interval.

.. report:: TSS.TSSClosestUpstream
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,100000
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:

   Histogram of distances to closest upstream TSS 

.. report:: TSS.TSSClosestUpstream
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,5000
   :tf-range: 0,100000,100
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:

   Histogram of distances to closest upstream TSS

Closest downstream
------------------

The following plots show the distance of each 
interval to the closest TSS that is downstream
of the intervals.

.. report:: TSS.TSSClosestDownstream
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,100000
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:

   Histogram of distances to closest downstream TSS

.. report:: TSS.TSSClosestDownstream
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,5000
   :tf-range: 0,100000,100
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:

   Histogram of distances to closest downstream TSS


