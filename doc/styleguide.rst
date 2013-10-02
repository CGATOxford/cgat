.. _styleguide:

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

   * `PEP0257 - Docstring Conventions
     <http://www.python.org/dev/peps/pep-0257/>`_

For documenting CGAT code, we follow the conventions for documenting
python code:

   * `Python Developer's guide <http://docs.python.org/devguide/documenting.html>`_

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

   * Scripts should follow the 
     `unix philosophy <http://en.wikipedia.org/wiki/Unix_philosophy>`_.
     They should concentrate on one task and do it well. Ideally,
     the major input and output can be read from and written to standard
     input and standard output, respectively. 

   * The names of scripts should be meaningful. Most of our scripts
     perform data transformation of one kind of another, these are
     often called ``a2b.py``. The distinctions can be subtle.
     Examples are:
     
     :doc:`scripts/gtf2gtf`
        Input is :term:`gtf`, output is :term:`gtf`. This script
        manipulates gene sets (filtering, merging, ...).

     :doc:`scripts/gtf2gff`
        Input is :term:`gtf`, output is :term:`gff`. This script
	takes gene sets and changes the hierarchical description
	within a :term:`gtf` file to the flat description of features
	in a :term:`gff` file. For example, this script can define
	gene territories, regulatory domains or genomic annotations
	based on a gene set.
 
     :doc:`scripts/bed2gff`
        Input is :term:`bed`, output is :term:`gff`. As both 
      	formats describe intervals in the genome, this script
        basically does a conversion between the two formats.

     Quite a few scripts contain the ``2table`` or ``2stats``. These
     compute, respectively, properties or summary statistics for
     entries in a file. For example:

     :doc:`scripts/gtf2table` 
         Input is :term:`gtf`. For each gene or transcript, compute
	 selected properties. If there are 10,000 genes in the input,
	 the output table will contain 10,000 rows.
	 
     :doc:`scripts/gff2stats`
         Input is :term:`gff`. Compute summary statistics across
	 all features in the file. Here, aggregate sizes or similar
	 by feature type or name per chromosome. No matter if there
	 are 10,000 or 100,000 interval is the input, the output 
	 will be have the same number of rows.

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

* In order to build documentation, each script, module and pipeline needs to
  be importable. Thus, make sure that when your pipeline depends on
  specific files, it does not fail when imported but not executed.
  
Documentation
=============

Writing doc-strings
-------------------

Functions should be documented through their doc-string using
restructured text. For example::

    def computeValue( name, method, accuracy=2):  

        :param name: The name to use.
    	:type name: str.
	:param method: method to use.
	:type state: choice of ('empirical', 'parametric')
	:param accuracy:
	:type accuracy: integer
	:returns:  int -- the value
	:raises: AttributeError, KeyError

Writing documentation for scripts
---------------------------------

Please follow the example in :doc:`cgat_script_template` for
documenting scripts. In addition, please pay attention to the following:

* Declare input data types for genomic data sets in optparse using 
  the `metavar` keyword. For example::

      parser.add_option( "--extra-intervals", dest = "extra_intervals",
                      metavar="bed", help = "..." )

  Setting the type permits the script to be integrated into workflow
  sytemns such as galaxy_.

* Please provide a meaningful example in the command line help.

* Be verbose. Something that is not documented within a script
  will not be used.

* Add meaningful tags to your scripts (``:Tags:``) so that they can be grouped into
  categories. Please choose from the following controlled
  vocabulary. If needed, additional terms can be added to this list.

  * Genomics
  * NGS
  * MultipleAlignment
  * GenomeAlignment
  * Intervals
  * Genesets
  * Sequences
  * Statistics
  * Summary
  * Variants
  * Protein













