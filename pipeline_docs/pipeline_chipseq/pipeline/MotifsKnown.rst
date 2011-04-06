Presence of canonical motifs
============================

MAST is used to look for the presence of canonical motifs in the dataset.
Canonical motifs are known beforehand.

Summary
-------

.. report:: Motifs.MastSummary
   :render: table
   :force: 

   Number of intervals with MAST matches.

Location of motifs within intervals
-----------------------------------

Absolute distance
+++++++++++++++++
.. report:: Motifs.MastMotifLocation
   :render: line-plot
   :transform: histogram
   :as-lines:
   :tf-aggregate: normalized-total,cumulative
   :layout: column-3
   :width: 300

   Location of motifs within intervals. If several
   motifs are within an interval, the midpoint
   of all motifs is used. The x-axis shows the
   distance of the motif to the peak.

Relative distance
+++++++++++++++++

The following plots show the relative location of motifs within
intervals. These should follow a sigmoidal
curve with the turning point at 0.

.. report:: Motifs.MastMotifLocationMiddle
   :render: line-plot
   :transform: histogram
   :as-lines:
   :tf-aggregate: normalized-total,cumulative
   :layout: column-3
   :width: 300

   Location of motifs within intervals. If several
   motifs are within an interval, the midpoint
   of all motifs is used. The x-axis shows the
   distance of the motif to the middle of the 
   interval.

Control intervals
+++++++++++++++++

The following plots show the relative location of motifs within
*control* intervals. These plots should show a straight line.

.. report:: Motifs.MastControlLocationMiddle
   :render: line-plot
   :transform: histogram
   :as-lines:
   :tf-aggregate: normalized-total,cumulative
   :layout: column-3
   :width: 300

   Location of motifs within *control* intervals.
   If several motifs are within an interval, the midpoint
   of all motifs is used. The x-axis shows the
   distance of the motif to the peak.

Motifs and peak strength
------------------------

.. report:: Motifs.MastPeakValWithMotif
   :render: line-plot
   :groupby: slice
   :as-lines:
   :layout: column-3
   :width: 300
  
   Proportion of intervals with a certain peakvalue or higher
   that contain a motif.

Motifs and interval locations
-----------------------------

.. report:: Motifs.AnnotationsMotifs
   :render: matrix-plot
   :layout: column-4
   :width: 300

   This plot shows the number of intervals with or without motif
   and their location.

.. report:: Motifs.AnnotationsPeakVal
   :render: matrix-plot
   :layout: column-4
   :width: 300

   This plot shows the number of intervals at a certain location
   together with the binding strength (:term:`peakval`)

Number of motifs per interval
-----------------------------

The following table shows stats on the number of motifs per interval.

.. report:: Motifs.MastNumberOfMotifs
   :render: table
   :transform: stats

   Number of motifs per interval

The following table shows histograms with the number of motifs per interval
for each motif and dataset.

.. report:: Motifs.MastNumberOfMotifs
   :render: table
   :transform: histogram
   :tf-bins: arange(0,20,1)

   Number of motifs per interval

MAST Evalue curves
------------------

In order to assess the validity of a match and motif,
same-sized segments on either side of each interval was submitted
to MAST. 

.. report:: Motifs.MastEvalues                                                                                                                                                                                                               
   :render: line-plot                                                                                                                                                                                                                        
   :transform: histogram                                                                                                                                                                                                                     
   :tf-aggregate: normalized-total,cumulative                                                                                                                                                                                                
   :logscale: x                                                                                                                                                                                                                              
   :groupby: track                                                                                                                                                                                                                           
   :as-lines:                                                                                                                                                                                                                                
   :layout: column-5
   :width: 200                                                                                                                                                                                                                         
                                                                                                                                                                                                                                             
   Cumulative distribution of Evalues in intervals 
   and in control regions.

MAST FDR curves
---------------

The following plots fdr values against evalues. The
fdr has been calibrated using the control intervals.

.. report:: Motifs.MastFDR
   :render: line-plot
   :logscale: x
   :yrange: 0,1.1
   :groupby: track
   :as-lines:
   :layout: column-4
   :width: 300

   Fdr values against evalues.

MAST curves
-----------

.. report:: Motifs.MastCurve                                                                                                                                                                                                               
   :render: line-plot                                                                                                                                                                                                              
   :as-lines:
   :layout: column-3
   :width: 300
   :yrange: 0,
                                                                                                                                                                                                                                          
   MAST curves

ROC curves
----------

The ROC curves test several :term:`interval` features for 
their ability to enrich for intervals that contain a motif.

The table below lists the area-under-the-curve (AUC). Higher values
are better values. Values of around 0.5 indicate no predictive value, and
less than 0.5 indicate bad predictors.

.. report:: Motifs.MastAUC                                                                                                                                                                                                                   
   :render: matrix
   :format: %5.2f 

   Table with AUC values.

And here are the ROC curves:

.. report:: Motifs.MastROC
   :render: line-plot
   :as-lines:
   :layout: column-3
   :width: 300

   ROC curves for MAST motifs and interval selection.
