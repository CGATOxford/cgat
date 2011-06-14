=======================
Building CGAT pipelines
=======================

This section describes how CGAT pipelines can be constructed using the
:mod:`Pipeline` module. The Pipeline.py module contains a
variety of useful functions for pipeline construction.

See the script :file:`pipeline_template.py` in the :term:`source directory`
as a starting point for pipelines.

.. _PipelineCommands:

Running commands within tasks
=============================

To run a command line program within a pipeline task, build a statement and call the
:meth:`Pipeline.run` method::

   @files( '*.unsorted', suffix('.unsorted'), '.sorted')
   def sortFile( infile, outfile ):

       statement = '''sort %(infile)s > %(outfile)s'''
       P.run()

On calling the :meth:`Pipeline.run` method, the environment of the caller
is examined for a variable called ``statement``. The variable is subjected
to string substitution from other variables in the local namespace. In the
example above, ``%(infile)s`` and ``%(outfile)s`` are substituted with the
values of the variables ``infile`` and ``outfile``, respectively.

The same mechanism also permits setting configuration parameters, for example::

   @files( '*.unsorted', suffix('.unsorted'), '.sorted')
   def sortFile( infile, outfile ):

       statement = '''sort -t %(tmpdir)s %(infile)s > %(outfile)s'''
       P.run()

will automatically substitute the configuration parameter ``tmpdir``
into the command. See ConfigurationValues_ for more on using configuration
parameters.

The pipeline will stop and return an error if the command exits with an error code.

If you chain multiple commands, only the return value of the last command
is used to check for an error. Thus, if an upstream command fails, it will go unnoticed.
To detect these errors, insert the ``checkpoint`` statement between commands. For example::

   @files( '*.unsorted.gz', suffix('.unsorted.gz'), '.sorted)
   def sortFile( infile, outfile ):

       statement = '''gunzip %(infile)s %(infile)s.tmp; 
                      checkpoint;
		      sort -t %(tmpdir)s %(infile)s.tmp > %(outfile)s;
		      checkpoint;
		      rm -f %(infile)s.tmp
       P.run()

Of course, the statement aboved could be executed more efficiently using pipes::

   @files( '*.unsorted.gz', suffix('.unsorted.gz'), '.sorted.gz')
   def sortFile( infile, outfile ):

       statement = '''gunzip < %(infile)s 
		      | sort -t %(tmpdir)s 
		      | gzip > %(outfile)s'''
       P.run()

The pipeline inserts code automatically to check for error return codes if multiple commands are
combined in a pipe.

Running commands on the cluster
-------------------------------

In order to run commands on cluster, use ``to_cluster=True``.

To run the command from the previous section on the cluster::

   @files( '*.unsorted.gz', suffix('.unsorted.gz'), '.sorted.gz')
   def sortFile( infile, outfile ):

       to_cluster = True
       statement = '''gunzip < %(infile)s 
		      | sort -t %(tmpdir)s 
		      | gzip > %(outfile)s'''
       P.run()

The pipeline will automatically create the job submission files, submit
the job to the cluster and wait for its return. 

Pipelines will use the command line options ``--cluster-queue``,
``--cluster-priority``, etc. for global job control. For example, to
change the priority when starting the pipeline, use::

   python <pipeline_script.py> --cluster-priority=-20

To set job options specific to a task, you can define additional variables::

   @files( '*.unsorted.gz', suffix('.unsorted.gz'), '.sorted.gz')
   def sortFile( infile, outfile ):

       to_cluster = True
       job_queue = 'longjobs.q'
       job_priority = -10
       job_options= "-pe dedicated 4 -R y" 
 
       statement = '''gunzip < %(infile)s 
		      | sort -t %(tmpdir)s 
		      | gzip > %(outfile)s'''
       P.run()

The above statement will be run in the queue ``longjobs.q`` at a priority of ``-10``.
Additionally, it will be executed in the parallel environment ``dedicated`` with at
least 4 cores.

Array jobs can be controlled through the ``job_array`` variable::

   @files( '*.in', suffix('.in'), '.out')
   def myGridTask( infile, outfile ):

       job_array=(0, nsnps, stepsize)
   
       statement = '''grid_task.bash %(infile)s %(outfile)s
          > %(outfile)s.$SGE_TASK_ID 2> %(outfile)s.err.$SGE_TASK_ID
       '''
       P.run()


Note that the :file:`grid_task.bash` file must be grid engine aware. This means
it makes use of the :env:`SGE_TASK_ID``, :env:`SGE_TASK_FIRST``, :env:`SGE_TASK_LAST` and 
:env:`SGE_TASK_STEPSIZE` environment variables to select the chunk of data it wants
to work on.

The job submission files are files called `tmp*` in the :term:`working directory`.
These files will be deleted automatically. However, the files will remain after 
aborted runs to be cleaned up manually.

.. _PipelineTracks:

Tracks
======

.. _PipelineDatabases:

Databases
=========

Loading data into the database
------------------------------

:mod:`Pipeline.py` offers various tools for working with databases. By default,
it is configured to use an sqlite3 database in the :term:`working directory` called :file:`csvdb`.

Tab-separated output files can be loaded into a table using the :meth:`Pipeline.load`
function. For example::

   @transform( 'data_*.tsv.gz', suffix('.tsv.gz'), '.load' )
   def loadTables( infile, outfile ):
      P.load( infile, outfile )

The task above will load all tables ending with ``tsv.gz`` into the database
Table names are given by the filenames, i.e, the data in :file:`data_1.tsv.gz` will
be loaded into the table :file:`data_1`.

The load mechanism uses the script :file:`csv2db.py` and can be configured using
the configuration options ``database`` and ``csv2db_options``. Additional options
can be given via the optional *options* argument::

   @transform( 'data_*.tsv.gz', suffix('.tsv.gz'), '.load' )
   def loadTables( infile, outfile ):
      P.load( infile, outfile, "--index=gene_id" )

Connecting to a database
------------------------

To use data in the database in your tasks, you need to first connect to the 
database. It helps to encapsulate the connection in a separate function. For 
example::

    def connect():
	dbh = sqlite3.connect( PARAMS["database"] )
	statement = '''ATTACH DATABASE '%s' as annotations''' % (PARAMS["annotations_database"])
	cc = dbh.cursor()
	cc.execute( statement )
	cc.close()

	return dbh

The above function will connect to the database. It will also attach a secondary database
``annotations``.

The following example illustrates how to use the connection::

    @transform( ... )
    def buildCodingTranscriptSet( infile, outfile ):

	dbh = connect()

	statement = '''SELECT DISTINCT transcript_id FROM transcript_info WHERE transcript_biotype = 'protein_coding' '''
	cc = dbh.cursor()
	transcript_ids = set( [x[0] for x in cc.execute(statement)] )
	...

.. _PipelineReports:

Reports
=======

build_report
update_report
sphinxreport.ini
conf.py

Including other pipelines


.. _ConfigurationValues:

Configuration values
====================

Setting up configuration values
--------------------------------

Pipelines are configured via a configuration script. The
following snippet can be included at the beginning of a pipeline
to set it all up::

   # load options from the config file
   import Pipeline as P
   P.getParameters( 
          ["%s.ini" % __file__[:-len(".py")],
	  "../pipeline.ini",
	  "pipeline.ini" ] )
   PARAMS = P.PARAMS

Configuration parameters will be read first from the 
file named :file:`pipeline_<pipeline_name>.ini` in the :term:`source directory`.
These sets all configuration values to default paramteres.

Next, the file :file:`../pipeline.ini` will be read (if it exists) and configuration
values that are specific to a certain project will overwrite default values.

Finally, run specific configuration will be read from the file :file:`pipeline.ini`
in the :term:`working directory`.

The method :meth:`Pipeline.getParameters` reads parameters and updates a global
dictionary of parameter values. It automatically guesses the type of parameters
in the order of ``int()``, ``float()`` or ``str()``.

If a configuration variable is empty (``var=``), it will be set to ``None``.

Configuration values from another pipeline can be added in a separate namespace::

   PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
   		                                 "pipeline_annotations.py" )

The statement above will load the parameters from a `pipeline_annotations`_ project
in :file:`annotations_dir`.

Using configuration values
--------------------------

Configuration values are accessible via the :var:`PARAMS` variable. The :var:`PARAMS`
variable is a dictionary mapping configuration parameters to values. Keys are in the
format ``section_parameter``. For example, the key ``bowtie_threads`` will provide the
configuration value of::

   [bowtie]
   threads=4

In a script, the value can be accessed via ``PARAMS["bowtie_threads"]``. 

Undefined configuration values will throw a :class:`ValueError`. To test if
a configuration variable exists, use::

   if 'bowtie_threads' in PARAMS: pass
      
To test, if it is unset, use::

   if 'bowie_threads' in PARAMS and not PARAMS['botwie_threads']: pass

Task specific parameters
------------------------

Task specific parameters can be set by creating a task specific section in
the :file:`pipeline.ini`. The task is identified by the output filename.
For example, given the following task::

   @files( '*.fastq', suffix('.fastq'), '.bam')
   def mapWithBowtie( infile, outfile ):
      ...

and the files :files:`data1.fastq` and :file:`data2.fastq` in the :term:`working directory`,
two output files :file:`data.bam` and :file:`data2.bam` will be created on executing
``mapWithBowtie``. Both will use the same parameters. To set parameters specific to the 
execution of :file:`data1.fastq`, add the following to :file:`pipeline.ini`::

   [data1.fastq]
   bowtie_threads=16

This will set the configuration value ``bowtie_threads`` to 16 when using the command
line substitution method in :meth:`Pipeline.run`. To get an task-specific parameter values in
a python task, use::

   @files( '*.fastq', suffix('.fastq'), '.bam')
   def mytask( infile, outfile ):
       MY_PARAMS = P.substituteParameters( locals() )
       
Thus, task specific are implemented generically using the :meth:`Pipeline.run` mechanism,
but pipeline authors need to explicitely code for track specific parameters.

.. _PipelineDocumentation:

Documentation
=============

Modular pipelines
=================

Exporting data
==============
