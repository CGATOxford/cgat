.. _PipelineReports:

================
Pipeline reports
================

CGAT pipelines use SphinxReport_ to report the outcome of a pipeline
run. Conceptually, the workflow is that a CGAT pipeline creates data
and uploads it into a database. SphinxReport_ then creates a report
from the database.

Background
==========



Advanced topics
===============

Linking to other reports
------------------------

Conditional content
===================

The ifconfig_ extension allows to include content depending on configuration
values. To use this extension you will need to modify
:file:`conf.py`. The example below shows the modifications implemented
in :doc:`pipelines/pipeline_mappping` to permit the conditional
inclusion of sections of the report depending on the mapper chosen::

    # add sphinx.ext.ifconfig to the list of extensions
    extensions.append( 'sphinx.ext.ifconfig' )
    
    # define a new configuration variable
    ################################################################
    # Add custom configuration variables for ifconfig extension
    def setup(app):
    	app.add_config_value('MAPPERS', '', True)

    # Set the value of custom configuration variables
    import CGAT.Pipeline as P
    P.getParameters(
	["%s/pipeline.ini" % os.path.splitext(__file__)[0],
	     "../pipeline.ini",
	     "pipeline.ini" ] )

    MAPPERS = P.asList( P.PARAMS["mappers" ] )
    
The thus defined and set custom configuration value ``MAPPERS`` can
now be used inside an rst document::

   .. toctree::
      :maxdepth: 2

      pipeline/Methods.rst
      pipeline/Status.rst
      pipeline/Mapping.rst
      pipeline/MappingSummary.rst
      pipeline/MappingContext.rst
      pipeline/MappingAlignmentStatistics.rst
      pipeline/MappingComplexity.rst
      .. ifconfig:: "tophat" in MAPPERS:
         pipeline/MappingTophat.rst
      .. ifconfig:: "star" in MAPPERS:
	 pipeline/MappingStar.rst
      .. ifconfig:: "tophat" in MAPPERS or "star" in MAPPERS or "gsnap" in MAPPERS
	 pipeline/Validation.rst

Referring to other reports
--------------------------

The intersphinx_ extension permits referring to other
sphinxreport documents. To use this extension, add the following to
your :file:`conf.py` configuration file::

    extensions = [ ..., 'sphinx.ext.intersphinx', ...]

    intersphinx_mapping = {'<identifier>': ('<>', None) }

where identifier is a suitable identifier and ``absolute path name to html`` is
the absolute location of the html build of the sphinx document you want
to refer to. This directory should contain the :file:`objects.inv` file. The
file is automatically created by sphinx, but sphinx needs to be run at least
once.

To refer to the other documentation, type::

   :ref:`My link to another documentation <identifier:label>`

Where ``label`` needs to be a valid identifier in the referred to document.


.. _intersphinx: http://sphinx-doc.org/ext/intersphinx.html
.. _ifconfig: http://sphinx-doc.org/ext/ifconfig.html



