*************************
Transcription start sites
*************************

This part of the pipeline examines the distances of intervals
to transcription start sites. Transcription start sites are
taken from the :term:`reference gene set`.

Overlap with TSS
----------------

Number of intervals that overlap a :term:`tss`.

.. report:: TSS.TSSOverlap
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total

   Number of intervals that overlap a :term:`tss`.

Distance to closest TSS
-----------------------

The following two plots show the distance of each 
interval to its closest :term:`tss`. The 
direction of the :term:`tss` is ignored.
Both plots show the same data at different scales.

An interval that overlaps a :term:`tss` is recorded
as having a distance of 0.

.. report:: TSS.TSSClosest
   :render: line-plot
   :transform: histogram
   :xrange: 0,100000
   :yrange: 0,1
   :tf-aggregate: normalized-total,cumulative
   :logscale: x
   :as-lines:
   :groupby: all

   Histogram of distances of intervals to closest :term:`tss`.

.. report:: TSS.TSSClosest
   :render: line-plot
   :transform: histogram
   :xrange: 0,5000
   :tf-range: 0,100000,100
   :yrange: 0,1
   :tf-aggregate: normalized-total,cumulative
   :logscale: x
   :as-lines:
   :groupby: all

   Histogram of distances of intervals to closest :term:`tss`.

Closest upstream TSS
--------------------

The following two plots show the distance of each 
interval to the closest :term:`tss` in upstream direction
looking from the :term:`tss`.

.. report:: TSS.TSSClosestUpstream
   :render: line-plot
   :transform: histogram
   :xrange: 0,100000
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:
   :groupby: all

   Histogram of distances to closest upstream :term:`tss`.

.. report:: TSS.TSSClosestUpstream
   :render: line-plot
   :transform: histogram
   :xrange: 0,5000
   :tf-range: 0,100000,100
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:
   :groupby: all

   Histogram of distances to closest upstream :term:`tss`.

Closest downstream
------------------

The following plots show the distance of each 
interval to the closest TSS that is downstream
of the intervals.

.. report:: TSS.TSSClosestDownstream
   :render: line-plot
   :transform: histogram
   :xrange: 0,100000
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:
   :groupby: all

   Histogram of distances to closest downstream TSS

.. report:: TSS.TSSClosestDownstream
   :render: line-plot
   :transform: histogram
   :xrange: 0,5000
   :tf-range: 0,100000,100
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:
   :groupby: all

   Histogram of distances to closest downstream TSS

Statistical association
=======================

I computed the statistical significance of intervals with transcription start site.
Briefly, the distance of intervals to the closest :term:`TSS` upstream or downstream
was computed both for the observed intervals and also for a set of simulated intervals.

.. Note:
   Intervals overlapping a TSS were removed. If they are kept in they overwhelm
   the signals as indeed, there is a large proportion of intervals overlapping 
   a TSS.

Distance to TSS
---------------

The ``tss`` workspace includes more or less the full genome, but segmented by
the presence of a transcription start site.

.. report:: TSS.TSSDistances
   :render: table

   Table with significance results. The statistical significance tests if the median
   of the distribution is closer than expected. View the plots for a more detailed	
   analysis.

Intergenic workspace
--------------------

The ``intronic`` workspace includes all genomic segments that are between
protein coding genes.

.. report:: TSS.IntergenicDistances
   :render: table

   Table with significance results. The statistical significance tests if the median
   of the distribution is closer than expected. View the plots for a more detailed	
   analysis.

Intronic workspace
------------------

The ``intronic`` workspace includes all genomic segments that are covered by protein
coding genes introns. Thus one gene might contribute several segments.

.. report:: TSS.IntronicDistances
   :render: table

   Table with significance results. The statistical significance tests if the median
   of the distribution is closer than expected. View the plots for a more detailed	
   analysis.

Genic workspace
---------------

The ``genic`` workspace includes all genomic segments that are covered by protein
coding genes - exons and introns.

.. report:: TSS.GenicDistances
   :render: table

   Table with significance results. The statistical significance tests if the median
   of the distribution is closer than expected. View the plots for a more detailed	
   analysis.
