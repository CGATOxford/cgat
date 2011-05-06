==================
Transcripts models
==================

Classification of transcript models
====================================

Transcripts are compared to a reference gene set. 

Transcript classes
------------------

.. report:: Genemodels.TranscriptClassCountsSummaryByClass
   :render: interleaved-bar-plot
   :layout: column-3
   :width: 300

   Class of transcripts

Protein coding classes
----------------------

The following plot lists the various classes of transcripts overlapping protein coding genes.

.. report:: Genemodels.TranscriptClassCounts
   :render: interleaved-bar-plot
   :slices: protein_coding
   :layout: column-3
   :width: 300

   Overlap with protein coding transcripts

Overlap per transcript source
-----------------------------

.. report:: Genemodels.TranscriptClassCountsSummaryBySource
   :render: interleaved-bar-plot
   :layout: column-3
   :width: 300

   Overlap with transcripts of a certain source

Cuffcompare classes
===================

The tables below show the absolute and relative frequencies of transfrags in each of these classes.
See :ref:`CuffCompareCodes` for an explanation of the codes.

.. report:: Genemodels.GeneModelsCodes
   :render: matrix
   :transform-matrix: normalized-row-total

   Table with location codes of transcripts

.. report:: Genemodels.GeneModelsCodes
   :render: table

   Table with location codes of transcripts

Number of expressed transfrags
------------------------------

The following tables show how many transfrags are present within replicates
within each experiment.

.. report:: Genemodels.GeneModelsSharedTransfrags
   :render: table
   :force:

   Table with number of experiments per transfrags

Number of expressed loci
------------------------

.. report:: Genemodels.GeneModelsSharedLoci
   :render: table
   :force:

   Table with number of experiments per locus

