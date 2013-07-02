===============
Medips analysis
===============

The R package MEDIPS is a tool to calculate methylation profiles from Medip-Seq
data. This page summarizes some quality statistics from that package.

Quality control
===============

Saturation
----------

The following plots examine saturation. From the manual:

"""
The saturation analysis addresses the question, whether the number of input
regions is sufficient to generate a saturated and reproducible methylation pro-
file of the reference genome. The main idea is that an insufficent number of
short reads will not result in a saturated methylation profile. Only if there is a
sufficient number of short reads, the resulting genome wide methylation profile
will be reproducible by another independent set of a similar number of
short reads.
"""

Basically, the reads are split into two sets of equal size. From each
set, an increasing fraction is selected and the correlation of the
genome-wide methylation profile is computed. 

For the estimated saturation, the set of reads is artifically doubled 
and the analysis is performed as before.

.. report:: Medips.MedipsPlots        
   :render: table
   :layout: column-3
   :width: 200
   :slices: saturation
   
   Saturation analysis. The plots show the genome-wide correlation of
   methylation profiles for increasing numbers of reads. The red curve
   stops at the total number of reads in the experiment, the blue
   curve at half that number.

Coverage
--------

The following plots examine the coverage of CpG islands achievd. From
the manual:

"""
The main idea of the coverage analysis is to test the number of CpGs
covered by the given short reads and to have a look at the depth of
coverage. 

For the coverage analysis, the total set of available regions is
divided into distinct random subsets of equal size, where the number of subsets is
determined by the parameter no_iterations (default=10). The coverage analysis
iteratively selects an increasing number of subsets and and tests how many
CpGs are covered by the available regions. Moreover, it is tested how many
CpGs are covered at least 1x, 2x, 3x, 4x, 5x, and 10x.

"""

.. report:: Medips.MedipsPlots        
   :render: table
   :layout: column-3
   :width: 200
   :slices: cpg_coverage
   
   CpG coverage

Calibration
-----------

.. report:: Medips.MedipsPlots        
   :render: table
   :layout: column-3
   :width: 200
   :slices: calibration
   
   Calibration fit
