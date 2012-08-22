=====================================
CAPseq Liver Testes Unique Intervals
=====================================

The following table presents the number of intervals that were unique to each dataset:

.. report:: macs_replicated_liver_testes_unique_intervals.replicatedUniqueIntervals
   :render: table

   Intervals unique to an individual dataset

Genomic Location
-----------------

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

.. report:: macs_replicated_liver_testes_unique_intervals.replicatedUniqueIntervalTranscriptOverlap
   :render: matrix 

   Summary of all genomic annotations

.. report:: macs_replicated_liver_testes_unique_intervals.replicatedUniqueIntervalTranscriptOverlap
   :render: interleaved-bar-plot

   Chart of all genomic annotations

.. report:: macs_replicated_liver_testes_unique_intervals.replicatedUniqueIntervalTranscriptOverlap
   :render: pie-plot
   :layout: column-2

   Chart of all genomic annotations

Per Gene TSS Annotation
-------------------------------

.. report:: macs_replicated_liver_testes_unique_intervals.replicatedUniqueIntervalGeneOverlap
   :render: matrix 

   Summary of all genomic annotations

.. report:: macs_replicated_liver_testes_unique_intervals.replicatedUniqueIntervalGeneOverlap
   :render: interleaved-bar-plot

   Chart of all genomic annotations

.. report:: macs_replicated_liver_testes_unique_intervals.replicatedUniqueIntervalGeneOverlap
   :render: pie-plot
   :layout: column-2

   Chart of all genomic annotations

Length
------

The following plot shows the distribution of interval length for each set.

.. report:: macs_replicated_liver_testes_unique_intervals.replicatedUniqueIntervalLengths
   :render: line-plot
   :transform: histogram
   :groupby: all
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:

   Distribution of interval lengths

Average Coverage
----------------

The following plot shows the distribution of average interval coverage for each set.
The average coverage is the average number of reads covering the bins that constitute the interval.

.. report:: macs_replicated_liver_testes_unique_intervals.replicatedUniqueIntervalAverageValues
   :render: line-plot
   :transform: histogram
   :groupby: all
   :tf-range: 0,50
   :tf-aggregate: normalized-total,reverse-cumulative
   :as-lines:

   Distribution of the average number of reads
   matching to bins within an interval.

Maximum Coverage
----------------

The following plot shows the maximum interval coverage for each set.
The maximum coverage is the maximum number of reads falling into the
bins that constitute an interval. The interval peak is the position with the maximum
number of reads.

.. report:: macs_replicated_liver_testes_unique_intervals.replicatedUniqueIntervalPeakValues
   :render: line-plot
   :transform: histogram
   :groupby: all
   :tf-range: 0,100
   :tf-aggregate: normalized-total,reverse-cumulative
   :as-lines:

   Distribution of the number of reads at the peak within an interval.
   The distribution list the proportion of intervals of a certain peak
   value or more.

Fold Change
-----------

The following plot shows the fold change over control (input) for each set.

.. report:: macs_replicated_liver_testes_unique_intervals.replicatedUniqueIntervalFoldChange
   :render: line-plot
   :transform: histogram
   :groupby: all
   :tf-range: 0,100
   :tf-aggregate: normalized-total,reverse-cumulative
   :as-lines:

   Distribution of fold enrichment for interval compared to control.

Closest TSS
-----------

The following plot shows the distance to the closest transcriptional start site for each set.

.. report:: macs_replicated_liver_testes_unique_intervals.replicatedUniqueIntervalTSS
   :render: line-plot
   :transform: histogram
   :groupby: all
   :xrange: 0,100000
   :yrange: 0,1
   :tf-range: 0,1000000,100
   :tf-aggregate: normalized-total,cumulative
   :as-lines:

   Distribution of distance to the closest transcriptional start site


CpG Density
-----------

The following plot shows the distribution of CpG density for each set.

.. report:: macs_replicated_liver_testes_unique_intervals.replicatedUniqueIntervalCpGDensity
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total
   :groupby: all
   :as-lines:

   Distribution of CpG density


CpG Observed/Expected
----------------------

The following plots show the distribution of observed/expected CpGs for each set.
The expected number of CpG dinucleotides was calculated as the product of the number of C and G nucleotides 
in the interval divided by the interval length as in Emboss cpgplot.
The control dataset was generated by taking an interval of the same size 10kb upstream of the CpG island.

.. report:: macs_replicated_liver_testes_unique_intervals.replicatedUniqueIntervalCpGObsExp
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total
   :groupby: all
   :as-lines:

   Distribution observed/expected CpGs (expected = nC*nG/length)


GC Content
------------

The following plot shows the distribution of GC content for each set.

.. report:: macs_replicated_liver_testes_unique_intervals.replicatedUniqueIntervalGCContent
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total
   :groupby: all
   :as-lines:

   Distribution of GC content


