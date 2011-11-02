===============
Medips analysis
===============

Medips is a tool to calculate methylation profiles from Medip-Seq
data.

Quality control

Saturation
----------

The following plots examine saturation.

.. report:: Medips.MedipsPlots        
   :render: table
   :layout: column-3
   :width: 200
   :slices: saturation
   
   Saturation analysis

Coverage
--------

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
