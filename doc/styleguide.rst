===========
Style Guide
===========

Coding style
============

This style guide lays down coding conventions in the CGAT repository.
For new scripts, follow the guidelines below. 

As the repository has grown over years and several people contributed,
the style between scripts can vary. For older scripts, follow the style within a
script/module. If you want to apply the newer style, make consistent
changes across the script.

In general, we want to adhere to the following conventions:

    * *Variable names* are lower case throughout with underscores to
      separate words, such as ``peaks_in_interval = 0``

    * *Function names* start with a lower case character and a
        verb. Additional words start in upper case, such as
    	``doSomethingWithData()``

    * *Class names* start with an upper case character, additional words
        start again in upper case, such as ``class AFancyClass():``

    * *Class methods* follow the same convention as functions, such as
    	``self.calculateFactor()``

    * *Class attributes* follow the same convention as variables, such
        as ``self.factor``

    * *Global variables* - in the rare cases they are used, are upper case
        throughout such as ``DEBUG=False``

    * *Module names* should start with an uppercase letter, for example,
        ``TreeTools.py`` in order to distinguish them from built-in
	and third-party python modules.

    * *Script names* are lower-case throughout with underscores to
      separate words, for example, ``bam2geneprofile.py`` or
      ``join_table.py``.

    * *Cython extensions* to scripts (via pyximport) should be put
      into the script name starting with an underscore. For example,
      The extensions to ``bam2geneprofile.py`` are in
      ``_bam2geneprofile.pyx``.

For new scripts, use the template :file:`script_template.py`.

The general rule is to write easily readable and maintainable
code. Thus, please

   * document code liberally and accurately
   * make use of whitespaces and line-breaks to break long statements
      into easily readable statements.

In case of uncertainty, follow the python style guides as much as
possible. The relevant documents are:

   * `PEP0008 - Style Guide for Python Code <http://www.python.org/dev/peps/pep-0008/>`_
   * `PEP0257 - Docstring Conventions <http://www.python.org/dev/peps/pep-0257/>`_

In terms of writing scripts, we follow the following conventions:

   * Each script should define the ``-h`` and ``--help`` options to
     give command line help usage.

   * For tabular output, scripts should output :term:`tsv` formatted
     tables. In these tables, records are separated by new-line
     characters and fields by tab characters. Lines with comments are started
     by the ``#`` character and are ignored. The first uncommented line
     should contain the column headers. For example::

        # This is a comment
	gene_id	length
	gene1	1000
	gene2	2000
     	# Another comment

Where to put code
=================

Different parts of the code base go into separate directories.

Scripts
   Scripts are python code that contains a main() function and
   are intended to be executed. Scripts go into the directory 
   :file:`/scripts`

Modules
  Modules contain supporting code and are imported by scripts or
  other modules. Modules go into the directory :file:`/CGAT`.

Pipelines
  Pipeline scripts and modules go into the directory :file:`/CGATPipelines`.

Pipelines
================

All components of a pipeline should go into the :file:`CGATPipelines`
directory. The basic layout of a pipeline is::

   CGATPipelines/pipeline_example.py
                /PipelineExample.py
                /PipelineExample.R
                /pipeline_example/pipeline.ini
                                 /conf.py
                                 /sphinxreport.ini
   

pipeline_example.py
    The main pipeline code. Pipelines start with the word ``pipeline``
    and follow the conventions for *script names*, all lower case with
    underscores separating words.

pipeline_example/pipeline.ini
    Default values for pipeline configuration values.

pipeline_example/conf.py
    Configuration script for sphinxreport.

pipeline_example/sphinxreport.ini
    Configuration script for sphinxreport.

pipeline_docs/pipeline_example
    Sphinxreport for pipeline.

PipelineExample.py
    Python utility methods and classes specific to this pipeline. Once
    methods and classes are shared between pipelines, consider moving
    them to a separate module.

PipelineExample.R
    R utility functions specific to this pipeline.

* Make sure that the pipeline.ini file exists and contains example/default
  values with annotation.

* Make sure that the pipeline can be imported from any directory,
  especially those not containing any data files or configuration
  files. This is important for the documentation of the pipeline
  to be built.

Other guidelines
================

* Only add source code and required data to the repository. Do
  not add .pyc files, backup files created by your editor or other 
  files.

* In order to build documentation, each script and module needs to
  be importable from anywhere. 
  












