===============
Reproducibility
===============

This section examines how transcript compare across biological replicates.
Shown are the pearson correlation coefficients of expressed transfrags for each 
pairwise comparison of tracks. Transfrags missing in one experiment but
present in another are ignored for computing the correlation, but are counted
in the absence/presence calls.

See :ref:`CuffCompareCodes` for an explanation of the codes.

Summary
=======

The following tables lists the proportion of transfrags that
are present in each pairwise comparison.

.. report:: Genemodels.TransfragReproducibility
   :render: interleaved-bar-plot
   :layout: column-5
   :width: 200
   :groupby: track
   
   Summary of reproducibility

   *pcalled*  percent of transfrags present in both samples
   *correlation* correlation coefficient

.. report:: Genemodels.TransfragReproducibility
   :render: table
   
   Summary of reproducibility

   *pcalled*  percent of transfrags present in both samples
   *correlation* correlation coefficient

Full table
==========

+------------------------------+---------------------------------------------------------------------------------+
|Column                        |Content                                                                          |
+------------------------------+---------------------------------------------------------------------------------+
|pairs                         |Number of transfrags in :term:`aggregate` track                                  |
+------------------------------+---------------------------------------------------------------------------------+
|both_null                     |transfrags absent in both tracks                                                 |
+------------------------------+---------------------------------------------------------------------------------+
|null1                         |transfrags absent in first track                                                 |
+------------------------------+---------------------------------------------------------------------------------+
|null2                         |transfrags absent in second track                                                |
+------------------------------+---------------------------------------------------------------------------------+
|not_null                      |transfrags present in both tracks                                                |
+------------------------------+---------------------------------------------------------------------------------+
|one_null                      |transfrags present in one track, but not the other                               |
+------------------------------+---------------------------------------------------------------------------------+
|coeff                         |pearson correlation coefficient of FPKM values of expressed transfrags           |
+------------------------------+---------------------------------------------------------------------------------+


The following table lists the correlation of expression.

See :ref:`CuffCompare`.

.. report:: Genemodels.TransfragCorrelation
   :render: table
   
   Correlation of expression levels.  


