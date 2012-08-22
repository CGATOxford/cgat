===========================================
Replicated CAPseq Interval Genomic Features
===========================================

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

.. report:: macs_replicated_interval_genomic_features.replicatedIntervalTranscriptOverlap
   :render: matrix 

   Summary of all genomic annotations

.. report:: macs_replicated_interval_genomic_features.replicatedIntervalTranscriptOverlap
   :render: interleaved-bar-plot

   Chart of all genomic annotations

.. report:: macs_replicated_interval_genomic_features.replicatedIntervalTranscriptOverlap
   :render: pie-plot
   :layout: column-2

   Chart of all genomic annotations

Per Gene TSS Annotation
-------------------------------

.. report:: macs_replicated_interval_genomic_features.replicatedIntervalGeneOverlap
   :render: matrix 

   Summary of all genomic annotations

.. report:: macs_replicated_interval_genomic_features.replicatedIntervalGeneOverlap
   :render: interleaved-bar-plot

   Chart of all genomic annotations

.. report:: macs_replicated_interval_genomic_features.replicatedIntervalGeneOverlap
   :render: pie-plot
   :layout: column-2

   Chart of all genomic annotations

