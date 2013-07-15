This section provides some background on CGAT pipelines.

Background
============

There really are two types of pipelines. In ``production pipelines`` the inputs are usually
the same every time the pipeline is run and the output is known beforehand. For example, 
read mapping and quality control is a typical pipeline. These pipelines can be well optimized
and can be re-used with little change in configuration.

``analysis pipelines`` control scientific analyses and are much more in a state of flux. 
Here, the input might change over time as the analysis expands and the output will change
with every new insight or new direction a project takes. It will be still a pipeline as long as
the output can be generated from the input without manual intervention. These pipelines leave
less scope for optimization compared to ``production pipelines`` and adapting a pipeline to
a new project will involve significant refactoring.

In CGAT, we are primarily concerned with ``analysis pipelines``, though we have some 
``production pipelines`` for common tasks.

There are several ways to build pipelines. For example, there are generic workflow
systems like `taverna <http://www.taverna.org.uk>`_ which even provide GUIs for connecting
tasks. A developer writes some glue code permitting the output of one application to
be used as input for another application. Also, there are specialized workflow systems 
for genomics, for example `galaxy <http://galaxy.psu.edu>`_, which allows you to save and share
analyses. New tools can be added to the system and new data imported easily for example
from the UCSC genome browser.

Flexibility
   There always new tools and insights. A pipeline should be ultimately 
   flexible and not constraining us in the things we can do.

Scriptability
   The pipeline should be scriptable, i.e, the whole pipeline can be run within
   another pipeline. Similarly, parts of a pipeline can be duplicated to process 
   several data streams in parallel. This is a crucial feature in genome studies
   as a single analysis will not permit making inferences by itself. For example,
   consider you find in ChIP-Seq data from a particular transcription factor that
   it binds frequently in introns. You will need to run the same analysis on 
   data from other transcription factors in order to assess if intronic binding is
   remarkable.

Reproducibility
   The pipeline is fully automated. The same inputs and configuration will produce
   the same outputs.

Reusability
   The pipeline should be able to be re-used on similar data, maybe only requiring 
   changes to a configuration file.

Archivability
   Once finished, the whole project should be able to archived without too many
   major dependencies on external data. This should be a simple process and hence
   all project data should be self-contained. It should not involve going through 
   various directories or databases to figure out which files and tables belong
   to a project or a project depends on.

There probably is not one toolset to satisfy all these criteria.. We use the following 
tools to build a pipeline:

   * ruffus_ to control the main computational steps
   * sqlite_ to store the results of the computational steps
   * sphinxreport_ to visualize the data in the sqlite database

.. _ruffus: http://www.ruffus.org.uk/
.. _sqlite: http://www.sqlite.org/
.. _sphinxreport: http://code.google.com/p/sphinx-report/
