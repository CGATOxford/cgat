=================================
Gapped motif analysis
=================================

Motif discovery
===============

GLAM2 was run to detect motifs that might contain gaps. The
following table lists the motifs found by GLAM2. A good result
should contain several similar motifs among the top results.

.. report:: Motifs.Glam2Results
   :render: table
   
   Table with GLAM2 results

Motif detection
===============

GLAM2SCAN is used to look for the presence of motifs in the dataset.

Summary
-------

The following table lists the number of intervals that can be called
positive based on scannig with GLAM2SCAN. 

.. report:: Motifs.GlamSummary
   :render: table
 
   Number of intervals with Glam matches. Cutoffs are set such that
   the FDR is 0.2.

GLAM FDR curves
---------------

In order to assess the validity of a match and motif,
same-sized segments on either side of each interval was submitted
to GLAM2SCAN.

.. report:: Motifs.GlamFDR
   :render: line-plot                                                                                                                                                                                                                        
   :groupby: track                                                                                                                                                                                                                           
   :as-lines:                                                                                                                                                                                                                                
   :layout: column-5
   :yrange: 0,1.1
   :width: 200                                                                                                                                                                                                                         
                                                                                                                                                                                                                                             
   Plot of the FDR against the score.																													    

GLAM Score curves
------------------

In order to assess the validity of a match and motif,
same-sized segments on either side of each interval was submitted
to GLAM2SCAN.

.. report:: Motifs.GlamScores
   :render: line-plot                                                                                                                                                                                                                        
   :transform: histogram                                                                                                                                                                                                                     
   :tf-aggregate: normalized-total,cumulative                                                                                                                                                                                                
   :groupby: track                                                                                                                                                                                                                           
   :as-lines:                                                                                                                                                                                                                                
   :layout: column-5
   :width: 200                                                                                                                                                                                                                         
                                                                                                                                                                                                                                             
   Cumulative distribution of scores in intervals 
   and in control regions.
