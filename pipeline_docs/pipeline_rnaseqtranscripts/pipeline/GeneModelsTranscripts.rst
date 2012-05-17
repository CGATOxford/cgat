==================
Transcripts models
==================

.. _Classification of transcript models:

Classification of transcript models
====================================

Transcripts are compared to a reference gene set. 

Transcript classes
------------------

The following table lists the class of transcript models 
ignoring any directional information. The class is derived
from the full :term:`reference gene set`.

.. report:: Genemodels.TranscriptClassCountsSummaryByClass
   :render: interleaved-bar-plot
   :layout: column-3
   :width: 300
   :transform: select
   :tf-fields: ntranscripts
   :transform-matrix: normalized-row-total

   Class of transcripts

.. report:: Genemodels.TranscriptClassCountsSummaryByClass
   :render: table
   :transform: select
   :tf-fields: ntranscripts

   Class of transcripts

Protein coding classes
----------------------

This section lists the various classes of transcript models overlapping protein 
coding genes in the :term:`reference gene set` only.

.. report:: Genemodels.TranscriptClassCounts
   :render: interleaved-bar-plot
   :slices: protein_coding
   :layout: column-3
   :width: 300
   :transform: select
   :tf-fields: ntranscripts
   :transform-matrix: normalized-row-total

   Overlap with protein coding transcripts

.. report:: Genemodels.TranscriptClassCounts
   :render: table
   :slices: protein_coding
   :transform: select
   :tf-fields: ntranscripts

   Overlap with protein coding transcripts

Overlap per transcript source
-----------------------------

This section lists the various sources that transcript models overlap in
the :term:`reference gene set`.

.. report:: Genemodels.TranscriptClassCountsSummaryBySource
   :render: interleaved-bar-plot
   :layout: column-3
   :width: 300
   :transform: select
   :tf-fields: ntranscripts
   :transform-matrix: normalized-row-total

   Overlap with transcripts of a certain source

.. report:: Genemodels.TranscriptClassCountsSummaryBySource
   :render: table
   :transform: select
   :tf-fields: ntranscripts

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
   :groupby: all
   :force:

   Table with number of experiments per transfrags

Number of expressed loci
------------------------

.. report:: Genemodels.GeneModelsSharedLoci
   :render: table
   :groupby: all
   :force:

   Table with number of experiments per locus

