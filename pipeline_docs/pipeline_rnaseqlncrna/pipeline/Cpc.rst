

===============================
Coding potential calculations
===============================

The coding potential calculator is run on a transcript level basis in order to identify
transcripts have potentail open reading frames. Filtering of coding transcripts is performed on
the level of the gene i.e if any transcript for a particular gene shows evidence for an ORF then
the all transcripts for that gene are filtered out.


Number of predicted coding and non-coding transcripts
======================================================

.. report:: CPC.CPC
   :render: table
   
   Number of predicted coding and non-coding transcripts


.. report:: CPC.CPC
   :render: interleaved-bar-plot
   

   Number of predicted coding and non-coding transcripts




Characteristics of coding and noncoding transcripts
====================================================

The coding potential calculator can be somewhat biased towards longer transcripts. We therefore investigate
some of the properties of both sets of transcripts that are predicted by the CPC to be coding and non-coding.


.. report:: CPC.CPCLength
   :render: table
   

   average length of transcripts and genes of predicted coding and non-coding transcripts


.. report:: CPC.CPCLengthDistribution
   :render: line-plot
   :transform: histogram
   :layout: column-2
   :as-lines:

   length distribution


.. report:: CPC.CPCLengthCumFreq
   :render: line-plot
   :as-lines:
   :layout: column-2   


   length cumulative frequency


.. report:: CPC.CPCExonNumber
   :render: table
   

   number of exons


.. report:: CPC.CPCExonNumber
   :render: interleaved-bar-plot
   

   number of exons


.. report:: CPC.CPCScoreCorellation
   :render: scatter-plot
   :regression: 1

   scatterplot of the relationship between transcript length and coding potential score



Classification of coding and non-coding transcripts
====================================================

.. report:: CPC.CPCClass
   :render: pie-plot
   :layout: column-2 

   classification of coding and non-coding transcripts



