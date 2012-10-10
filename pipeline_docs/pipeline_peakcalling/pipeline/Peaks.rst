===========
Peak shapes
===========

This page plots the read density under peaks.

.. report:: Intervals.PeakShapeTracker
   :render: matrixNP-plot
   :width: 300
   :nolabel-cols:
   :nolabel-rows:
   :slices: peak_height
   :max-rows: 0
   :max-cols: 0
   :layout: column-2

   Peak shapes sorted by peak height.

.. report:: Intervals.PeakShapeTracker
   :render: matrixNP-plot
   :width: 300
   :nolabel-cols:
   :nolabel-rows:
   :slices: peak_width
   :transform-matrix: normalized-row-total
   :max-rows: 0
   :max-cols: 0
   :layout: column-2

   Peak shapes sorted by peak width and
   normalized by height.
