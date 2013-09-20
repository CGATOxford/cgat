===========================
The CGAT code collection
===========================

The CGAT code collection has grown out of the work in comparative
genomics by the Ponting group in the last decade. Now, CGAT_ has added
functionality to do next-generation sequencing analysis.

The CGAT Code collection has two components. The first component
is a collection of scripts, the CGAT tools. The second component
is a collection of pipelines. While both components are part of this
collection, we are currently concentrating on publishing the CGAT
tools.

For questions, please subscribe and contact us at the 
`CGAT user group
<https://groups.google.com/forum/?fromgroups#!forum/cgat-user-group>`_.

CGAT Tools
==========

The CGAT tools can be installed from pypi::

   pip install cgat

To use CGAT Tools, use the ``cgat`` front-end. For example, to
strip sequence and quality information from a bam_ file, type:

   cgat bam2bam --strip=sequence < in.bam > out.bam

Documentation of CGAT tools is available 
`here <http://www.cgat.org/~andreas/documentation/cgat/cgat.html#cgat>`_.

CGAT Pipelines
==============

We have developed numerous pipelines in comparative genomics
and NGS analysis. The pipelines are generally available and should
be fairly portable. Some documentation of the pipelines is 
`here <http://www.cgat.org/~andreas/documentation/cgat/Pipelines.html#pipelines>`_.

Note that we currently are not able to fully support and document the 
pipelines. They are under continuous development and changing rapidly.
However, they might give some ideas or building blocks when developing
your own pipelines.

.. _bam: http://en.wikipedia.org/wiki/SAMtools
.. _CGAT: http://www.cgat.org


