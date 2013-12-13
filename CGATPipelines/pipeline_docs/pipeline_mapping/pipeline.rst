.. Test documentation master file, created by
   sphinx-quickstart on Mon Mar 23 15:27:57 2009.

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
   python/Trackers.rst

.. ifconfig:: "tophat" in MAPPERS

   .. toctree::
      pipeline/MappingTophat.rst

.. ifconfig:: "star" in MAPPERS

   .. toctree::
      pipeline/MappingStar.rst

.. ifconfig:: "tophat" in MAPPERS or "star" in MAPPERS or "gsnap" in MAPPERS

   .. toctree::
      pipeline/Validation.rst
 
.. errorlist::

.. warninglist::
  
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


