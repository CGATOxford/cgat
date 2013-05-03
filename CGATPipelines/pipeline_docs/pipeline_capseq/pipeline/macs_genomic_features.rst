=========================
MACS Genomic Features
=========================

Intervals were annotated using a reference gene set. Briefly, transcripts 
of all ENSEMBL protein coding genes were merged per gene. The gene transcriptional start site (TSS)
was taken as the first base of the merged gene. 
Intervals were classified according to their overlap with these features. Classification was hierarchical, 
such that if an interval overlapped two types of feature (for example from overlapping gene models) then 
that interval was annotated with the term that is highest in the tree. The hierarchy is outlined in the table below.

+---------------+---------------------------------------------------------------------------------+
|Term           | Definition                                                                      |
+---------------+---------------------------------------------------------------------------------+
|TSS            |Located within 1kb of a transcriptional start site of a protein-coding gene      |
+---------------+---------------------------------------------------------------------------------+
|Upstream       |Not any of the above and within 5kb upstream of a protein-coding gene            |
+---------------+---------------------------------------------------------------------------------+
|Gene           |Overlapping a gene region (intron/exon/utr) but not a TSS region                 |
+---------------+---------------------------------------------------------------------------------+
|Downstream     |Not any of the above and within 5kb downstream of a protein-coding gene          |
+---------------+---------------------------------------------------------------------------------+
|Intergenic     |None of the above. At least 5kb from the nearest protein coding gene             |
+---------------+---------------------------------------------------------------------------------+


Per Ensembl Protein-coding Transcript
======================================

.. report:: macs_genomic_features.mergedIntervalEnsemblTranscriptOverlap
   :render: matrix 

   Summary of all genomic annotations based on Ensembl transcript TSS

.. report:: macs_genomic_features.mergedIntervalEnsemblTranscriptOverlap
   :render: interleaved-bar-plot

   Chart of all genomic annotations based on Ensembl transcript TSS

.. report:: macs_genomic_features.mergedIntervalEnsemblTranscriptOverlap
   :render: pie-plot
   :layout: column-2

   Chart of all genomic annotations based on Ensembl transcript TSS


Per Ensembl Protein-coding Gene
======================================

.. report:: macs_genomic_features.mergedIntervalEnsemblGeneOverlap
   :render: matrix 

   Summary of all genomic annotations based on Ensembl gene TSS

.. report:: macs_genomic_features.mergedIntervalEnsemblGeneOverlap
   :render: interleaved-bar-plot

   Chart of all genomic annotations based on Ensembl gene TSS

.. report:: macs_genomic_features.mergedIntervalEnsemblGeneOverlap
   :render: pie-plot
   :layout: column-2

   Chart of all genomic annotations based on Ensembl gene TSS
   

Genomic Repeats
===============

The following plots show the number of binding intervals that overlap repeats.

.. report:: macs_genomic_features.RepeatOverlap
   :render: table
   :transform-matrix: normalized-row-total

   Proportion of intervals overlapping repeats

.. report:: macs_genomic_features.RepeatOverlap
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total

   Proportion of intervals overlapping repeats




