===
IDR
===

IDR is a method for automated thresholding in high-throughput
experiments. The IDR method uses the reproducibility between 
replicates in order to estimate the level at which noise overwhelms
the signal.

The following text from `here
<https://sites.google.com/site/anshulkundaje/projects/idr>` explains
it in more detail:

Reproducibility is essential to reliable scientific discovery in
high-throughput experiments. The IDR (Irreproducible Discovery Rate)
framework is a unified approach to measure the reproducibility of
findings identified from replicate experiments and provide highly stable
thresholds based on reproducibility. Unlike the usual scalar measures
of reproducibility, the IDR approach creates a curve, which
quantitatively assesses when the findings are no longer consistent
across replicates. In layman's terms, the IDR method compares a pair
of ranked lists of identifications (such as ChIP-seq peaks). These
ranked lists should not be pre-thresholded i.e. they should provide
identifications across the entire spectrum of high
confidence/enrichment (signal) and low confidence/enrichment
(noise). The IDR method then fits the bivariate rank distributions
over the replicates in order to separate signal from noise based on a
defined confidence of rank consistency and reproducibility of
identifications i.e the IDR threshold.

The method was developed by Qunhua Li and Peter Bickel's group and is
extensively used by the ENCODE and modENCODE projects. A publication
describing the statistical details of the IDR framework is 
http://www.stat.washington.edu/qli/IDR-AOAS.pdf

In our pipeline, we call peaks with SPP and run the IDR analysis on
these peaks.

.. report:: Tracker.TrackerImages
   :render: gallery-plot
   :glob: export/idr/*.pdf
   	     
   IDR diagnostic plots for all tracks

   
