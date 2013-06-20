.. _Annotations:

=============
 Annotations
=============

Transcript models are annotated using a reference gene set. The principal three 
classes are

known:
   transcripts overlapping the reference gene set.

novel:
   transcripts not overlapping the reference gene set.

ambiguous:
   transcripts that partially overlap the reference gene set.

Overview
========

.. report:: Trackers.AllAnnotations
   :render: table
   :slices: all  

   Table with all annotations

.. report:: Trackers.AnnotationsBases
   :render: matrix
   :slices: all
   :transform-matrix: normalized-row-max

   Bases overlapping annotations.

.. _AnnotationsFirstLevel:

First level classification
==========================

.. report:: Trackers.BasicAnnotations
   :render: stacked-bar-plot
   :slices: all  
   :transform-matrix: normalized-row-total

   First level annotations.

.. report:: Trackers.BasicAnnotations
   :render: table
   :slices: all  

   Table with the first level annotations.

.. report:: Trackers.BasicAnnotations
   :render: matrix
   :transform-matrix: normalized-row-total
   :slices: all  

   Table with the first level annotations.

.. _AnnotationsKnown:

Annotation of known transcripts
===============================

.. report:: Trackers.KnownAnnotations
   :render: stacked-bar-plot
   :slices: known  

   Annotations of known transcripts.

.. report:: Trackers.KnownAnnotations
   :render: stacked-bar-plot
   :slices: known  
   :transform-matrix: normalized-row-total

   Annotations of known transcripts.

.. report:: Trackers.KnownAnnotations
   :render: table
   :slices: known  

   Annotations of known transcripts.

.. report:: Trackers.KnownAnnotations
   :render: matrix
   :transform-matrix: normalized-row-max
   :slices: known  

   Annotations of known transcripts.

.. _AnnotationsNovel:

Annotation of novel transcripts
===============================

.. report:: Trackers.UnknownAnnotations
   :render: stacked-bar-plot
   :slices: unknown  

   Annotations of novel transcripts.

.. report:: Trackers.UnknownAnnotations
   :render: stacked-bar-plot
   :slices: unknown  
   :transform-matrix: normalized-row-total

   Annotations of novel transcripts.

.. report:: Trackers.UnknownAnnotations
   :render: matrix
   :transform-matrix: normalized-row-max
   :slices: unknown  

   Annotations of novel transcripts.

.. report:: Trackers.UnknownAnnotations
   :render: table
   :slices: unknown  

   Annotations of novel transcripts.

