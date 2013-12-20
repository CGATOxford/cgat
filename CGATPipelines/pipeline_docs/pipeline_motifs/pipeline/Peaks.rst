==========
Peak shape
==========

.. report:: Peaks.PeakShapeTracker
   :render: matrixNP-plot
   :nolabel-cols:
   :nolabel-rows:
   :slices: peak_height
   :max-rows: 0
   :max-cols: 0
   :layout: column-2

   Peak shapes sorted by peak height.

.. report:: Peaks.PeakShapeTracker
   :render: matrixNP-plot
   :nolabel-cols:
   :nolabel-rows:
   :slices: peak_width
   :transform-matrix: normalized-row-max
   :max-rows: 0
   :max-cols: 0
   :layout: column-2

   Peak shapes sorted by peak width and
   normalized by height.
