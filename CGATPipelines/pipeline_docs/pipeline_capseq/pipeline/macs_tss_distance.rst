==========================
Transcription Start Sites
==========================

This part of the pipeline examines the distances of CAPseq intervals
from transcription start sites. Here two sets of transcription start sites are considered: 
one TSS per protein-coding Ensembl transcript and one TSS per protein-coding Ensembl gene.
In the case of single TSS per gene the the first base of the merge gene model was selected as the representative 
TSS in the absence of tissue specific expression information.

Overlap with TSSs
==================

Number of transcription start sites that each intervals overlaps:

.. report:: macs_tss_distance.transcriptTSSOverlap
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total

   Transcript TSS

.. report:: macs_tss_distance.geneTSSOverlap
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total

   Gene TSS

Distance to closest TSS
=======================

The following plots show the genomic distance of each 
interval to the closest TSS. The frequency is normalised to the total number of intervals.

Transcript TSS
--------------

.. report:: macs_tss_distance.TSSClosest
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,100000
   :yrange: 0,1
   :tf-range: 0,1000000,100
   :tf-aggregate: normalized-total,cumulative
   :as-lines:
   :xtitle: Genomic_distance_(bp)

   Histogram of distances of interval to TSS (100bp intervals up to 100kb)

.. report:: macs_tss_distance.TSSClosest
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,5000
   :tf-range: 0,1000000,100
   :yrange: 0,1
   :tf-aggregate: normalized-total,cumulative
   :as-lines:
   :xtitle: Genomic_distance_(bp)

   Histogram of distances of interval to TSS (100bp intervals up to 5kb)

Gene TSS
--------

.. report:: macs_tss_distance.geneTSSClosest
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,100000
   :yrange: 0,1
   :tf-range: 0,1000000,100
   :tf-aggregate: normalized-total,cumulative
   :as-lines:
   :xtitle: Genomic_distance_(bp)

   Histogram of distances of interval to TSS (100bp intervals up to 100kb)

.. report:: macs_tss_distance.geneTSSClosest
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,5000
   :tf-range: 0,1000000,100
   :yrange: 0,1
   :tf-aggregate: normalized-total,cumulative
   :as-lines:
   :xtitle: Genomic_distance_(bp)

   Histogram of distances of interval to TSS (100bp intervals up to 5kb)

Closest Upstream TSS
=====================

The following plots show the distance of each 
interval to the closest TSS that is upstream
of the interval.

Transcript TSS
--------------

.. report:: macs_tss_distance.transcriptTSSClosestUpstream
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,100000
   :tf-range: 0,1000000,100
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:
   :xtitle: Genomic_distance_(bp)

   Histogram of distances to closest upstream TSS (100bp intervals up to 100kb)

.. report:: macs_tss_distance.transcriptTSSClosestUpstream
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,5000
   :tf-range: 0,1000000,100
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:
   :xtitle: Genomic_distance_(bp)

   Histogram of distances to closest upstream TSS (100bp intervals up to 5kb)

Gene TSS
--------------

.. report:: macs_tss_distance.geneTSSClosestUpstream
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,100000
   :tf-range: 0,1000000,100
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:
   :xtitle: Genomic_distance_(bp)

   Histogram of distances to closest upstream TSS (100bp intervals up to 100kb)

.. report:: macs_tss_distance.geneTSSClosestUpstream
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,5000
   :tf-range: 0,1000000,100
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:
   :xtitle: Genomic_distance_(bp)

   Histogram of distances to closest upstream TSS (100bp intervals up to 5kb)

Closest Downstream TSS
==========================

The following plots show the distance of each 
interval to the closest TSS that is downstream
of the intervals.

Transcript TSS
--------------

.. report:: macs_tss_distance.transcriptTSSClosestDownstream
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,500000
   :tf-range: 0,1000000,100
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:
   :xtitle: Genomic_distance_(bp)

   Histogram of distances to closest downstream TSS

.. report:: macs_tss_distance.transcriptTSSClosestDownstream
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,10000
   :tf-range: 0,1000000,100
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:
   :xtitle: Genomic_distance_(bp)

   Histogram of distances to closest downstream TSS (100bp intervals up to 10kb)

Gene TSS
--------------

.. report:: macs_tss_distance.geneTSSClosestDownstream
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,500000
   :tf-range: 0,1000000,100
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:
   :xtitle: Genomic_distance_(bp)

   Histogram of distances to closest downstream TSS

.. report:: macs_tss_distance.geneTSSClosestDownstream
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,10000
   :tf-range: 0,1000000,100
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:
   :xtitle: Genomic_distance_(bp)

   Histogram of distances to closest downstream TSS (100bp intervals up to 10kb)



