******************
Contamination
******************

This part of the pipeline examines the possibility of contamination.

.. report:: Trackers.ContaminationCoverage
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total

   Proportions of transcripts with >1 or exactly 1 read matching.

If repeats are assumed to be non-transcribed, the proportion of bases
due to genomic contamination can be estimated from extrapolating from
the proportion of transcribed bases repeats. 

.. report:: Trackers.ContaminationRepeats
   :render: table

   Estimates of contamination based on overlap with repeats

However the assumption that non-coding transcripts contain no repeats is not necessarily
correct. For example, genomic contaminants should contain no introns, but usually there are
transcript models that overlap repeats and are spliced.

