=========
Encode QC
=========

This page summarizes versus ENCODE quality control metrics.

The analysis follows the suggestion by Kundaje et al. (submitted)
using cross-correlation analysis for ChIP-Seq QC. More information
is `here <http://code.google.com/p/phantompeakqualtools/>`_.

Summary
=======

.. report:: PeakCalling.SPPQuality
   :render: table
   
   SPP Quality metrics

.. report:: PeakCalling.SPPQuality
   :render: interleaved-bar-plot
   :slices: nsc,rsc
   :function: 1.0, 1.1

   SPP Quality Metrics. The thresholds
   for good data for nsc (1.1) and rsc (1.0)
   are shown.

Diagnostic plots
================

.. report:: Tracker.TrackerImages
   :render: gallery-plot
   :glob: export/quality/*.pdf
   	     
   Diagnostic plots for all tracks

Full table
==========

bamfile
	tagAlign/BAM filename
mapped_reads
	effective sequencing depth i.e. total number of mapped reads
	in input file

estFragLen
	comma separated strand cross-correlation peak(s) in decreasing order of correlation.
          The top 3 local maxima locations that are within 90% of the
	  maximum cross-correlation value are output.
          In almost all cases, the top (first) value in the list
	  represents the predominant fragment length.

corr_estFragLen
	comma separated strand cross-correlation value(s) in
	decreasing order (estFragLen follows the same order)

phantomPeak
	read length/phantom peak strand shift

corr_phantomPeak
	Correlation value at phantom peak

argmin_corr
	strand shift at which cross-correlation is lowest
min_corr
	minimum value of cross-correlation

nsc 
    Normalized strand cross-correlation coefficient 
    (NSC) = corrEstFraglen / min_corr

rsc
    Relative strand cross-correlation coefficient 
    (RSC) = (corrEstFragLen - min_corr) / (corr_phantomPeak -
    min_corr)

quality
    Quality tag based on thresholded RSC (codes:
    -2:veryLow, -1:Low, 0:Medium, 1:High, 2:veryHigh)

.. report:: PeakCalling.SPPQuality
   :render: table
   
   SPP Quality metrics
