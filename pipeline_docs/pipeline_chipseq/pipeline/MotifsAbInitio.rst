=======================
Ab initio motif finding
=======================

MEME was used to find motifs within each data set.

To increase the signal/noise ratio, 200 bp segments centerered around the 10% of intervals
with the highest peak value are used. Repeats are masked to reduce the change of finding low
complexity repeats. MEME finds the top 10 motifs, permitting any number of
motifs within the sequence collection. 

The following table shows the number of sequences submitted for the motif search:

.. report:: Motifs.MemeRuns
   :render: table

   Meme results.

The following table lists all the motifs found.

.. report:: Motifs.MemeResults
   :render: table

   Meme results, showing for each motif the width,
   the number of times it is found, its evalue
   and its information content.

Comparison to known motifs
==========================

The following table lists the results of comparing the :term:`master motif`
against all motifs in the ab-initio motif search.

.. report:: Motifs.TomTomResults
   :render: table

   TomTom results. The :term:`master motif` is scanned against 
   the motifs returned by MEME.


Sequence composition
====================

The following boxplot shows the distribution G+C content in the
sequences used for motif discovery.

.. report:: Motifs.MemeInputSequenceComposition
   :render: box-plot
   :slices: pGC

   Distribution of G+C content in sequences used for motif discovery.
