=========================
MACS Genomic Annotations
=========================

Intervals are annotated using a reference gene set. Briefly, transcripts 
of all ENSEMBL protein coding genes are merged per gene. If two genes 
overlap, the longer one is taken. With the resultant list of exons, 
genomic segments are annotated as exon, CDS, UTR and flank, intronic, intergenic. Next,
intervals are classified according to their overlap with these segments. As intervals
are usually much longer than the expected binding site, the overlap for most
classes is only required to be partial. The classes are:

+---------------+---------------------------------------------------------------------------------+
|Term           | Definition                                                                      |
+---------------+---------------------------------------------------------------------------------+
|cds            |mostly part of a CDS. These are intervals within a protein-coding exon           |
+---------------+---------------------------------------------------------------------------------+
|utr            |mostly part of UTR. These are intervals within the UTR of a gene                 |
+---------------+---------------------------------------------------------------------------------+
|intergenic     |mostly intergenic. These are intervals within an intergenic region               |
|               |i.e. more than 5kb from the closest exon                                         |
+---------------+---------------------------------------------------------------------------------+
|upstream       |not any of the above and partly upstream of a gene. These are intervals that     |
|               |might overlap part of the UTR or the 5kb segment before to the 5'-terminal       |
|               |exon of a gene.                                                                  |
+---------------+---------------------------------------------------------------------------------+
|downstream     |not any of the above and partly downstream of a gene. These are intervals        |
|               |that might overlap part of the UTR or the 5kb segment after to the 3'-terminal   |
|               |exon of a gene.                                                                  |
+---------------+---------------------------------------------------------------------------------+
|intronic       |not any of the above and partly intronic. Note that these could also include     |
|               |promotors of short alternative transcripts that skip one or more of the first    |
|               |exons.                                                                           |
+---------------+---------------------------------------------------------------------------------+
|ambiguous      |none of the above                                                                |
+---------------+---------------------------------------------------------------------------------+


.. report:: macs_annotations.AllAnnotations
   :render: matrix 

   Summary of all genomic annotations

.. report:: macs_annotations.AllAnnotations
   :render: interleaved-bar-plot

   Chart of all genomic annotations

.. report:: macs_annotations.AllAnnotations
   :render: pie-plot
   :layout: column-2

   Chart of all genomic annotations

.. report:: macs_annotations.AnnotationsBases
   :render: matrix
   :transform-matrix: normalized-row-max

   Bases overlapping genomic annotations

Transcription Start Sites
=========================

This part of the pipeline examines the distances of intervals
to transcription start sites. Transcription start sites are
usually refered to distances to the :term:`reference gene set`.

Overlap with TSS
----------------

Number of TSS start sites that intervals overlap.

.. report:: macs_annotations.TSSOverlap
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total

   Number of TSS start sites that intervals overlap.

Distance to closest TSS
-----------------------

The following plots show the genomic distance of each 
interval to the closest TSS. The frequency is normalised to the total number of intervals.

.. report:: macs_annotations.TSSClosest
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

.. report:: macs_annotations.TSSClosest
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

The following plots show the genomic distance (bp) of each 
interval to the closest TSS.

.. report:: macs_annotations.TSSClosest
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,100000
   :yrange: 0,70000
   :tf-range: 0,1000000,100
   :tf-aggregate: cumulative
   :as-lines:
   :xtitle: Genomic_distance_(bp)
   :ytitle: Intervals

   Histogram of distances of interval to TSS (100bp intervals up to 100kb)

.. report:: macs_annotations.TSSClosest
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,5000
   :tf-range: 0,1000000,100
   :yrange: 0,20000
   :tf-aggregate: cumulative
   :as-lines:
   :xtitle: Genomic_distance_(bp)
   :ytitle: Intervals

   Histogram of distances of interval to TSS (100bp intervals up to 5kb)

Closest upstream TSS
--------------------

The following plots show the distance of each 
interval to the closest TSS that is upstream
of the interval.

.. report:: macs_annotations.TSSClosestUpstream
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,100000
   :tf-range: 0,1000000,100
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:
   :xtitle: Genomic_distance_(bp)

   Histogram of distances to closest upstream TSS 

.. report:: macs_annotations.TSSClosestUpstream
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,5000
   :tf-range: 0,1000000,100
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:
   :xtitle: Genomic_distance_(bp)

   Histogram of distances to closest upstream TSS

Closest downstream
------------------

The following plots show the distance of each 
interval to the closest TSS that is downstream
of the intervals.

.. report:: macs_annotations.TSSClosestDownstream
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,100000
   :tf-range: 0,1000000,100
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:
   :xtitle: Genomic_distance_(bp)

   Histogram of distances to closest downstream TSS

.. report:: macs_annotations.TSSClosestDownstream
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,5000
   :tf-range: 0,1000000,100
   :tf-aggregate: normalized-total,cumulative
   :yrange: 0,1
   :as-lines:
   :xtitle: Genomic_distance_(bp)

   Histogram of distances to closest downstream TSS

TSS Profile
------------

.. report:: macs_annotations.TSSProfile
   :render: line-plot
   :transform: histogram
   :groupby: track
   :xrange: -3000,3000
   :tf-range: -1000000,1000000,200
   :as-lines:
   :layout: column-2

   Density of CAP-seq intervals around TSS

TTS Profile
------------

.. report:: macs_annotations.TTSProfile
   :render: line-plot
   :transform: histogram
   :groupby: track
   :xrange: -3000,3000
   :tf-range: -1000000,1000000,200
   :as-lines:
   :layout: column-2

   Density of CAP-seq intervals around TTS

Genomic Repeats
===============

The following plots show the number of binding intervals that overlap repeats.

.. report:: macs_annotations.RepeatOverlap
   :render: table
   :transform-matrix: normalized-row-total

   Proportion of intervals overlapping repeats

.. _FigureRepeatOverlap:

.. report:: macs_annotations.RepeatOverlap
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total

   Proportion of intervals overlapping repeats

