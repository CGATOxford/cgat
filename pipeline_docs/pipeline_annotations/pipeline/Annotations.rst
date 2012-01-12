====================
Genomic annotations
====================

Using the protein coding transcripts of an ENSMEBL gene set, the
genome is divided and annotated into non-overlapping segments.

Briefly, transcripts of all ENSEMBL protein coding genes are merged per gene. If two genes 
overlap, the longer one is taken.

	
.. report:: AnnotationReport.GenomicFeatureCoverage
   :render: pie-plot
   
   Overall coverage of genome by genomic features

.. report:: AnnotationReport.ChromosomeFeatureCoverage
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total

   Coverage of chromosomes by genomic features
