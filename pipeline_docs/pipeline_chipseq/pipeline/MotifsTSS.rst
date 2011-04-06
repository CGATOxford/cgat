==========================================
Combining motifs and location of intervals
==========================================

The following plots show the number of intervals
with or without motifs that overlap or don't 
overlap transcription start sites.

With overlap
------------

.. report:: MotifsTSS.MotifsOverlappingTSS
   :render: stacked-bar-plot
   :layout: column-3

   Number of intervals overlapping transcription start
   sites and whether or not they contain motifs.

.. report:: MotifsTSS.MotifsOverlappingTSS
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total
   :layout: column-3

   Number of intervals overlapping transcription start
   sites and whether or not they contain motifs.

No overlap
----------

.. report:: MotifsTSS.MotifsNonOverlappingTSS
   :render: stacked-bar-plot
   :layout: column-3

   Number of intervals not overlapping transcription start
   sites and whether or not they contain motifs.

.. report:: MotifsTSS.MotifsNonOverlappingTSS
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total
   :layout: column-3

   Number of intervals not overlapping transcription start
   sites and whether or not they contain motifs.
