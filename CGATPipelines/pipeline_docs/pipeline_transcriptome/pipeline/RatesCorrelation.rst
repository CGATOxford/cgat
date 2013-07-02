**************************************
Full correlation matrix of all rates
**************************************

Statistics
----------

.. report:: Trackers.RatesAll
   :render: table
   :transform: correlation
   :groupby: track

Correlation coefficients
------------------------

.. report:: Trackers.RatesAll
   :render: matrix
   :transform: correlation,select 
   :tf-fields: coefficient
   :format: %6.4f  

PValues
-------

.. report:: Trackers.RatesAll
   :render: matrix
   :transform: correlation,select 
   :tf-fields: logpvalue
   :format: %6.4f  

