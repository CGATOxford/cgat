=============
Documentation
=============

Overview
========

CGAT scripts and modules use sphinx_ for documentation. The philosophy
is to maintain documentation and code together. Thus, most
documentation will be kept inside the actual scripts and modules,
supported by overview documents explaining usage and higher level
concepts.

Building the documentation
==========================

CGAT's documentation lives in the :file:`doc` directory of the
repository. To build the documentation, enter the :file:`doc`
directory and type::

   make html 

The output will be in the directory :file:`_build/html`. 

.. note::
   Each script, module and pipeline needs to be importable,
   i.e, the following must work::

       python -c "import pipeline_mapping"

   Especially in pipelines some care is necessary to avoid failing
   with an error if no input or configuration files are present.

Writing documentation
=====================

sphinx_ documentation is written in `Restructured Text`_. A useful
primer is `here <http://sphinx-doc.org/rest.html>`_.

Some specifics for the CGAT code base are:

* Refering to a separate script can be done using the ``:doc:``
  directive, for example::

     :doc:`scripts/bed2summary`

  Note that the path relative to the current directory needs to
  be supplied.
   
* Glossary terms (``:term:``) are defined in
  :file:`glossary.rst`.
  
Adding documentation
====================

In order to add a new script, module or pipeline to the documentation documement,
perform the following steps.

Here, we will be adding the script :file:`bed2summary.py` to
the documentation.

1. Create a file :file:`doc/scripts/bed2summary.rst` with the
   following contents::

      .. automodule:: bed2summary

      .. program-output:: python ../scripts/bed2summary.py --help

   This will build the documentation within the bed2summary script
   and add the command line help to the document.

2. Add an entry to :file:`doc/scripts.rst`. For example::

       .. toctree::

          scripts/bed2summary

   Please add your script to the toctree of an existing group.

3. For scripts that are part of the CGAT code collection, also add an
   entry into :file:`doc/CGATReference.rst`.

Adding a module or pipeline is similar to adding a script, except that:

1. the :file:`.rst` file should be in :file:`doc/modules` or
   :file:`doc/pipelines`, respectively.

2. The entry needs to be added to :file:`modules.rst` or 
   :file:`CGATPipelines.rst`, respectively.

3. no ``program-output`` is necessary.

Requisites
==========

Building the documentation requires the following components:

sphinx_
   The documenation building system.

sphinxcontrib-programoutput_
   Adding command line output to documenation.

.. _sphinxcontrib-programoutput: https://pypi.python.org/pypi/sphinxcontrib-programoutput
.. _Restructured Text: http://docutils.sourceforge.net/rst.html
