====================
Transcript modelling
====================

This section collections parameters helping to model transcription.
It shows various plots with respect to the reference gene set.

Intronic expression
===================

The following plot examines the correlation between
exonic expression and intronic expression. A good correlation
between intronic read depth with exonic read depth hints
at unprocessed mRNA transcripts.

.. .. report:: Reference.IntronicExonicReadDepth
..    :render: r-smooth-scatter-plot
..    :xrange: 1,
..    :yrange: 1,
..    :groupby: ungrouped
..    :layout: column-3
..    :width: 300
..    :logscale: xy

..    Smoothed scatter plots of exonic maximum read depth versus intronic 
..    maximum read depth.
   
.. UTR extension
.. =============

.. The following plots show a distribution of the UTR extension per data set and gene.

.. .. report:: Reference.UTRExtension
..    :render: r-smooth-scatter-plot
..    :layout: column-3
..    :width: 300	

..    Predicted length of UTR versus length of known UTR.

UTR densities
-------------

The following plots show the amount of extension of the UTR. It shows
the expression level of the last exon in a protein coding gene of the
:term:`reference` gene set together with the maximum coverage
in windows upstream and downstream for the first and last exon, respectively.

Expression values are maximum read counts.

+-------+----------------------------------------------+
|Section|Content                                       |
+-------+----------------------------------------------+
|raw    |raw counts scaled by loc10                    |
+-------+----------------------------------------------+
|scaled |counts scaled by expression level of exons    |
+-------+----------------------------------------------+
|fit    |fitting results                               |
+-------+----------------------------------------------+

.. .. report:: Reference.UTRReadDensityTable
..    :render: user

..    UTR extension






