.. _WritingReports:

========================
Writing pipeline reports
========================

CGAT pipelines use SphinxReport_ to report the outcome of a pipeline
run. Conceptually, the workflow is that a CGAT pipeline creates data
and uploads it into a database. SphinxReport_ then creates a report
from the database.

Background
==========

.. todo::

   Some text here about why sphinxreport

Advanced topics
===============

Conditional content
-------------------

The ifconfig_ extension allows to include content depending on configuration
values. To use this extension you will need to modify
:file:`conf.py`. The example below shows the modifications implemented
in :doc:`pipelines/pipeline_mapping` to permit the conditional
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

   .. ifconfig:: "tophat" in MAPPERS

      .. toctree::
	 pipeline/MappingTophat.rst

   .. ifconfig:: "star" in MAPPERS

      .. toctree::
	 pipeline/MappingStar.rst

   .. ifconfig:: "tophat" in MAPPERS or "star" in MAPPERS or "gsnap" in MAPPERS

      .. toctree::
	 pipeline/Validation.rst

Note that ``.. ifconfig`` needs to be a first level directive and
can not be include into another directive such as ``.. toctree``.

Referring to other reports
--------------------------

The intersphinx_ extension permits referring to other
sphinxreport documents. To use this extension you will need to modify
your :file:`conf.py` configuration file. For example::

    # add sphinx.ext.ifconfig to the list of extensions
    extensions.append( 'sphinx.ext.intersphinx' )

    # add mapping information
    intersphinx_mapping = {
       'readqc': ('/ifs/projects/proj013/readqc/report/html', None) ,
       'mapping1': ('/ifs/projects/proj013/mapping1/report/html', None),
       'mapping2': ('/ifs/projects/proj013/mapping2/report/html', None),
	}
	
This will link to three other reports. The three reports are
abbreviated as ``readqc``, ``mapping1`` and ``mapping2``. The paths
need to be the absolute location of the html build of the sphinx
documents you created previously. These directories should contain a
:file:`objects.inv` file which is usually automatically created by sphinx.

To refer to the other documentation, type::

    :ref:`My link to another documentation <identifier:label>`

``label`` is a valid identifier in the referred to
document. For example::

    :ref:`ReadQC <readqc:readqcpipeline>`

	ReadQC pipeline - fastqc

    :ref:`Unique Mapping  <mapping1:mappingpipeline>`

	Mapping pipeline - short read mapping with bwa. Only
	uniquely mapping reads are kept.

    :ref:`Non-unique mapping <mapping2:mappingpipeline>`

	Mapping pipeline - short read mapping with bwa with same
	parameters as above, but all reads are kept. 

.. _intersphinx: http://sphinx-doc.org/ext/intersphinx.html
.. _ifconfig: http://sphinx-doc.org/ext/ifconfig.html



