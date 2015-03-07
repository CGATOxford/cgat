.. _mappingpipeline:

===========================
Short Read Mapping Pipeline
===========================

Contents:

.. toctree::
   :maxdepth: 2

   pipeline/Methods.rst
   pipeline/Status.rst
   pipeline/Mapping.rst
   pipeline/MappingSummary.rst
   pipeline/MappingContext.rst
   pipeline/MappingAlignmentStatistics.rst
   pipeline/MappingComplexity.rst
   pipeline/ReferenceCoverage.rst

.. need to sort out variables. Need to be in conf.py
.. but there should be a generic way to push updates
.. to all reports.

.. ifconfig:: "tophat" in PARAMS['mappers']

   .. toctree::
      pipeline/MappingTophat.rst

.. ifconfig:: "star" in PARAMS['mappers']

   .. toctree::
      pipeline/MappingStar.rst

.. ifconfig:: "tophat" in PARAMS['mappers'] or "star" in PARAMS['mappers'] or "gsnap" in PARAMS['mappers']

   .. toctree::
      pipeline/RNASeq.rst
 
===================
Errors and Warnings
===================

.. errorlist::

.. warninglist::

===================
Pipeline Status
===================

.. report:: Tracker.TrackerImages
   :render: gallery-plot
   :glob: export/pipeline.svg

   Pipeline status

======================
Pipeline Documentation
======================

.. automodule:: pipeline_mapping
   :members:
   :inherited-members:
   :show-inheritance:
  
==================
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

