===================
Predicted CGIs
===================

Computationally predicted CpG islands were annotated using a reference gene set. Briefly, transcripts 
of all ENSEMBL protein coding genes were merged per gene. If two genes 
overlap, the longer one is taken. With the resultant list, 
genomic segments are annotated as exon, CDS, UTR and flank, intronic or intergenic. Next,
intervals are classified according to their overlap with these segments. As intervals
are frequently much longer than genomic segments, the overlap for most
classes is only required to be partial. The classes definitions are:

+---------------+---------------------------------------------------------------------------------+
|Term           | Definition                                                                      |
+---------------+---------------------------------------------------------------------------------+
|cds            |mostly part of a CDS. These are intervals within a protein-coding exon           |
+---------------+---------------------------------------------------------------------------------+
|utr            |mostly part of UTR. These are intervals within the UTR of a gene                 |
+---------------+---------------------------------------------------------------------------------+
|intergenic     |mostly intergenic. These are intervals within an intergenic region               |
|               |i.e. more than 5kb from the closest exon                                         |
+---------------+---------------------------------------------------------------------------------+
|upstream       |not any of the above and partly upstream of a gene. These are intervals that     |
|               |might overlap part of the UTR or the 5kb segment before to the 5'-terminal       |
|               |exon of a gene.                                                                  |
+---------------+---------------------------------------------------------------------------------+
|downstream     |not any of the above and partly downstream of a gene. These are intervals        |
|               |that might overlap part of the UTR or the 5kb segment after to the 3'-terminal   |
|               |exon of a gene.                                                                  |
+---------------+---------------------------------------------------------------------------------+
|intronic       |not any of the above and partly intronic. Note that these could also include     |
|               |promotors of short alternative transcripts that skip one or more of the first    |
|               |exons.                                                                           |
+---------------+---------------------------------------------------------------------------------+
|ambiguous      |none of the above                                                                |
+---------------+---------------------------------------------------------------------------------+

Genomic Annotation of Predicted CGIs
-------------------------------------

.. report:: predicted_cgis.cgiAnnotations
   :render: matrix 

   Summary of all genomic annotations

.. report:: predicted_cgis.cgiAnnotations
   :render: interleaved-bar-plot

   Chart of all genomic annotations

.. report:: predicted_cgis.cgiAnnotations
   :render: pie-plot

   Chart of all genomic annotations

Predicted CGIs Overlapping TSS
------------------------------

.. report:: predicted_cgis.cgitssoverlap
   :render: table

   Overlap between predicted CGIs and TSS

CpG Density of Predicted CGIs
-------------------------------

.. report:: predicted_cgis.CGI_CpGObsExp2
   :render: line-plot
   :transform: histogram
   :groupby: all
   :as-lines:

   Distribution observed/expected CpGs (expected = nC*nG/length)

GC Content of Predicted CGIs
------------------------------

.. report:: predicted_cgis.CGI_GCContent
   :render: line-plot
   :transform: histogram
   :groupby: all
   :as-lines:

   Distribution of GC content


