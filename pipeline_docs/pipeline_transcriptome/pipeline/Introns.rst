*************************************
Analysis of introns and splice sites
*************************************

Number of introns
=================

.. report:: Trackers.IntronsPerTranscript
   :transform: histogram
   :render: line-plot
   :tf-aggregate: normalized-total
   :as-lines:
   :tf-range: 0,,1

   Distribution of the number of introns per transcript model

.. report:: Trackers.IntronsPerTranscript
   :render: table
   :transform: stats
 
   Statistics of introns per transcript model

Splice sites in transcripts have been examined whether
they contain splice motifs and whether introns boundaries
correspond to the reference gene set.

.. _SpliceSites:

Splice sites
============

The following plots count the number of introns that have known
splice motifs.

.. report:: Trackers.SpliceSites
   :render: interleaved-bar-plot

   Counts of introns with splice motifs.

.. report:: Trackers.SpliceSites
   :render: table
   
   Table with counts of introns with splice motifs.

.. report:: Trackers.SpliceSites
   :render: matrix
   :transform-matrix: normalized-row-total

   Table with normalized counts of introns with splice motifs.

Intron boundaries
=================

The following plots count the number of introns that are
identical to the reference gene set.

+--------------------+--------------------+
|column              |content             |
+--------------------+--------------------+
|total               |introns in          |
|                    |transcript models   |
+--------------------+--------------------+
|found               |introns in          |
|                    |transcript models   |
|                    |present also in     |
|                    |reference gene set  |
+--------------------+--------------------+
|missed              |introns in          |
|                    |transcript models   |
|                    |not present in      |
|                    |reference gene set  |
+--------------------+--------------------+
|both                |both intron         |
|                    |boundaries match    |
+--------------------+--------------------+
|one                 |only one intron     |
|                    |boundary matches    |
+--------------------+--------------------+
|none                |no intron boundary  |
|                    |matches             |
+--------------------+--------------------+
|exon_skipping       |intron contains an  |
|                    |exon in the         |
|                    |reference gene set  |
+--------------------+--------------------+


.. report:: Trackers.IntronBoundaries
   :render: interleaved-bar-plot

   Histogram of intron boundaries.

.. report:: Trackers.IntronBoundaries
   :render: interleaved-bar-plot
   :transform-matrix: normalized-row-total

   Histogram of intron boundaries - normalized.

.. report:: Trackers.IntronBoundaries
   :render: table

   Table comparing the accuracy of intron/exon boundaries.

.. report:: Trackers.IntronBoundaries
   :render: matrix
   :transform-matrix: normalized-row-total

   Table comparing the accuracy of intron/exon boundaries.

Splice motifs
=============

This section lists the type of splice motifs found in 
transcript models.

.. report:: Trackers.SpliceMotifs
   :render: interleaved-bar-plot

   Counts of major splice motifs in introns.

.. report:: Trackers.SpliceMotifs
   :render: interleaved-bar-plot
   :transform-matrix: normalized-row-total

   Proportions of major splice motifs in introns.

.. report:: Trackers.SpliceMotifs
   :render: table

   Table with introns with splice motifs - absolute counts.

.. report:: Trackers.SpliceMotifs
   :render: matrix
   :transform-matrix: normalized-row-max

   Table with introns with splice motifs - normalized counts.
