======================
Predicted CpG Islands
======================

Genomic Annotation of Predicted CGIs
=====================================

Intervals were classified according to their overlap with these features. Classification was hierarchical, 
such that if an interval overlapped two types of feature (for example from overlapping gene models) then 
that interval was annotated with the term that is highest in the tree. The hierarchy is outlined in the table below.

+---------------+---------------------------------------------------------------------------------+
|Term           | Definition                                                                      |
+---------------+---------------------------------------------------------------------------------+
|TSS            |Per Transcript: Located within 1kb of a transcriptional start site of a          |
|               |protein-coding transcript                                                        |
|               |Per Gene: the gene transcriptional start site (TSS) interval was taken as the    |
|               |region from the first TSS to the last per gene + 1kb either side                 |
+---------------+---------------------------------------------------------------------------------+
|Gene           |Overlapping a gene region (intron/exon/utr) but not a TSS region                 |
+---------------+---------------------------------------------------------------------------------+
|Upstream       |Not any of the above and within 5kb upstream of a protein-coding gene            |
+---------------+---------------------------------------------------------------------------------+
|Downstream     |Not any of the above and within 5kb downstream of a protein-coding gene          |
+---------------+---------------------------------------------------------------------------------+
|Intergenic     |None of the above. At least 5kb from the nearest protein coding gene             |
+---------------+---------------------------------------------------------------------------------+

Per Transcript TSS Annotation
-------------------------------

.. report:: predicted_cgis.cgiTranscriptTSSOverlap
   :render: table

   Overlap between predicted CGIs and annotated protein-coding TSSs

.. report:: predicted_cgis.cgiTranscriptOverlap
   :render: matrix 

   Summary of all genomic annotations

.. report:: predicted_cgis.cgiTranscriptOverlap
   :render: interleaved-bar-plot

   Chart of all genomic annotations

.. report:: predicted_cgis.cgiTranscriptOverlap
   :render: pie-plot

   Chart of all genomic annotations

Per Gene TSS Annotation
-------------------------------

.. report:: predicted_cgis.cgiGeneTSSOverlap
   :render: table

   Overlap between predicted CGIs and TSS

.. report:: predicted_cgis.cgiGeneOverlap
   :render: matrix 

   Summary of all genomic annotations

.. report:: predicted_cgis.cgiGeneOverlap
   :render: interleaved-bar-plot

   Chart of all genomic annotations

.. report:: predicted_cgis.cgiGeneOverlap
   :render: pie-plot

   Chart of all genomic annotations

CpG Density of Predicted CGIs
==============================

.. report:: predicted_cgis.CGI_CpGObsExp
   :render: line-plot
   :transform: histogram
   :groupby: all
   :as-lines:

   Distribution observed/expected CpGs (expected = nC*nG/length)

GC Content of Predicted CGIs
=============================

.. report:: predicted_cgis.CGI_GCContent
   :render: line-plot
   :transform: histogram
   :groupby: all
   :as-lines:

   Distribution of GC content


