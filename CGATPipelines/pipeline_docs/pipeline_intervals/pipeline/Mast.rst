====
MAST
====

==================
Presence of motifs
==================

MAST is used to look for the presence of canonical motifs in the dataset.
Canonical motifs are known beforehand.

Summary
=======

The following plots show the percent of peaks that are explained by
the presence of a motif. The percentage of peaks explained is
estimated using MAST-curves (:pmid:`19160518`). MAST-curves define
an E-Value cutoff with the largest difference between peaks with motifs
and expected number of peaks with motif by chance.

.. report:: Motifs.MastQuickSummary
   :render: interleaved-bar-plot
   :transform: filter
   :tf-fields: MC-explained / %
   :layout: column-3
   :width: 200

   Percentage of peaks explained by motifs.

.. report:: Motifs.MastQuickSummary
   :render: table
   :transform: filter                      
   :tf-fields: MC-With-Motifs,MC-Evalue,MC-explained / %

   Percentage of peaks explained by motifs.


Specificity of Motifs
=====================

The following plots show the distribution of E-values of matches
within an interval and within a same-sized interval left and right of
the interval. A good motif/interval combination should have lower
evalues within intervals compared to left/right of an interval.

.. report:: Motifs.MastMotifEvalues
   :render: box-plot
   :layout: column-4
   :width: 200

   Distribution of Evalues.

Location of motifs within intervals
===================================

Motifs should be located centrally within predicted binding intervals.

The following plots show distance of motifs within intervals from the
middle of an interval.

A strong motif displays a sigmoidal curve, while a weak/unspecific
motif creates a diagonal.

[slow so disabled by default]

.. .. report:: Motifs.MastMotifLocationMiddle
..   :render: line-plot
..   :transform: histogram
..   :as-lines:
..   :tf-aggregate: normalized-total,cumulative
..   :layout: column-3
..   :width: 200

   Location of motifs within intervals. If several
   motifs are within an interval, the midpoint
   of all motifs is used. The x-axis shows the
   distance of the motif to the middle of the
   interval.

.. Control intervals
.. +++++++++++++++++

.. The following plots show the relative location of motifs within
.. *control* intervals, random genomic locations of the same size.
.. These plots should all show a straight line.

.. .. report:: Motifs.MastControlLocationMiddle
..    :render: line-plot
..    :transform: histogram
..    :as-lines:
..    :tf-aggregate: normalized-total,cumulative
..    :layout: column-3
..    :width: 200

..    Location of motifs within *control* intervals.
..    If several motifs are within an interval, the midpoint
..    of all motifs is used. The x-axis shows the
..    distance of the motif to the peak.

Motifs and peak strength
========================

The following plot shows the proportion of peaks of a certain height
that contain a motif. For each combination of track and motif there
are two lines:
   * proportion with motif: proportion of peaks of height equal or
     greater than X that contain a motif. This usually starts at 100%
     for high peaks (high values of X) but then drops off to the
     number of peaks that contain a motif.

   * recall: proportion of peaks of height equal or greater than
     X. This is a reverse cumulative distribution starting at low
     values for high X and then increasing to 100% as X decreases.

[slow so disabled by default]

.. .. report:: Motifs.MastPeakValWithMotif
..   :render: line-plot
..   :groupby: track
..   :as-lines:
..   :layout: column-3
..   :width: 200

   Proportion of intervals with a certain peakvalue or higher
   that contain a motif.

.. Motifs and interval locations
.. =============================

.. .. report:: Motifs.AnnotationsMotifs
..    :render: matrix-plot
..    :layout: column-4
..    :width: 300

..    This plot shows the number of intervals with or without motif
..    and their location.

.. .. report:: Motifs.AnnotationsPeakVal
..    :render: matrix-plot
..    :layout: column-4
..    :width: 300

..    This plot shows the number of intervals at a certain location
..    together with the binding strength (:term:`peakval`)

.. Number of motifs per interval
.. =============================

.. The following table shows stats on the number of motifs per interval.

.. .. report:: Motifs.MastNumberOfMotifs
..    :render: table
..    :transform: stats

..    Number of motifs per interval

.. The following table shows histograms with the number of motifs per interval
.. for each motif and dataset.

.. .. report:: Motifs.MastNumberOfMotifs
..    :render: table
..    :transform: histogram
..    :tf-bins: arange(0,20,1)

..    Number of motifs per interval


.. Distance from peak
.. ++++++++++++++++++

.. The following plots show distance of motifs within intervals from the
.. interval peak, the position with the largest number of reads.

.. A strong motif displays a sigmoidal curve, while a weak/unspecific
.. motif creates a diagonal.

.. .. report:: Motifs.MastMotifLocation
..    :render: line-plot
..    :transform: histogram
..    :as-lines:
..    :tf-aggregate: normalized-total,cumulative
..    :layout: column-3
..    :width: 200

..    Location of motifs within intervals. If several
..    motifs are within an interval, the midpoint
..    of all motifs is used. The x-axis shows the
..    distance of the motif to the peak.

