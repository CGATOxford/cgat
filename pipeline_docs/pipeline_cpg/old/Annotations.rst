===================
Genomic Annotations
===================

Intervals are annotated using a reference gene set. Briefly, transcripts 
of all ENSEMBL protein coding genes are merged per gene. If two genes 
overlap, the longer one is taken. With the resultant list of exons, 
genomic segments are annotated as exon, CDS, UTR and flank, intronic, intergenic. Next,
intervals are classified according to their overlap with these segments. As intervals
are usually much longer than the expected binding site, the overlap for most
classes is only required to be partial. The classes are:

+---------------+---------------------------------------------------------------------------------+
|Term           | Definition                                                                      |
+---------------+---------------------------------------------------------------------------------+
|cds            |mostly part of a CDS. These are intervals within a protein-coding exon           |
+---------------+---------------------------------------------------------------------------------+
|utr            |mostly part of UTR. These are intervals within the UTR of a gene                 |
+---------------+---------------------------------------------------------------------------------+
|intergenic     |mostly intergenic. These are intervals within an intergenic region               |
|               |i.e. more than 1kb from the closest exon                                         |
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


.. report:: Annotations.AllAnnotations
   :render: matrix 

   Summary of all genomic annotations

.. report:: Annotations.AllAnnotations
   :render: interleaved-bar-plot

   Chart of all genomic annotations

.. report:: Annotations.AnnotationsBases
   :render: matrix
   :transform-matrix: normalized-row-max

   Bases overlapping genomic annotations

