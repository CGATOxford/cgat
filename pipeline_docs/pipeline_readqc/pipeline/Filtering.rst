==================
Processing summary
==================

The following plot shows the reads input and output after read
processing.

.. report:: ReadqcReport.ProcessingSummary
   :render: interleaved-bar-plot
   :slices: ninput,noutput
   
   Summary of processing

.. report:: ReadqcReport.ProcessingSummary
   :render: interleaved-bar-plot
   :slices: percent_output
   
   Percent of reads output after processing

.. report:: ReadqcReport.ProcessingSummary
   :render: table
   :force:

   Summary of processing   

.. report:: ReadqcReport.ProcessingDetails
   :render: table
   :force:

   Summary of processing   

.. .. report:: ReadqcReport.FilteringSummary
..    :render: interleaved-bar-plot
..    :slices: processed_reads,unchanged_reads

..    Summary of filtering	

.. .. report:: ReadqcReport.FilteringSummary
..    :render: table

..    Summary of filtering	


