==================================
Importing CGAT scripts into galaxy
==================================

General Preparation
====================

Add :file:`/ifs/devel/cgat` to :env:`PYTHONPATH`.

Make sure that extensions have been built::

   python setup.py develop --multi-version

The following directories are important:

galaxy-dist
    Location of the galaxy distribution

cgat-xml
    CGAT directory within the galaxy distribution. Create by typing::
    
        mkdir <galaxy-dist>/tools/cgat

cgat-scripts
     The CGAT scripts directory.


Adding a script manually
========================

The following instructions describe the steps necessary to add a cgat
script to galaxy. 

For example, we want to publish the :file:`bam2stats.py`
script. First, create a file in :file:`<galaxy-dist>/tools/cgat` called
:file:`bam2stats.xml` with the following contents::

    <tool id="bam2stats.py" name="Compute Stats from BAM file">
      <description>Compute stats for a bam file</description>
	<command
          interpreter="python">/ifs/devel/cgat/scripts/bam2stats.py -v 0 &lt; $input &gt; $output
	</command>
	<inputs>
	   <param format="bam" name="input" type="data" label="Source file"/>
	</inputs>
	<outputs>
	   <data format="tabular" name="output" />
	</outputs>
      <help>
      Compute statistics for a bam file.
      </help>
    </tool>

Add an entry to :file:`tool_conf.xml` for the script::

  <section name="CGAT Tools" id="cgat_tools">
    <tool file="cgat/bam2stats.xml" />
  </section>

After restarting galaxy, the :file:`bam2stats` command should now be
visible in the ``CGAT`` section.

Automatic conversion of scripts
===============================

The CGAT tool collection contains a script called :doc:`cgat2rdf` that can create
and xml file for inclusion into galaxy. To create a wrapper for
:doc:`bam2stats`, run::

    python <cgat-scripts>cgat2rdf.py --format=galaxy <cgat-scripts>bam2stats.py > <cgat-xml>bam2stats.xml

As before, add an entry to :file:`tool_conf.xml` for the script.

For automatted conversion, a few rules need to be followed (see below).

Writing galaxy compatible scripts
---------------------------------

CGAT scripts have generally a call interface that is compatible with
galaxy and can thus be easily integrated. However, to make automatic
conversion as easy as possible, conforming to a few coding conventions
help.

1. Assign a metavar type to command line options of genomic file
   formats. For example::

      parser.add_option("-b", "--bam-file", dest="bam_files", type="string", metavar="bam",
                        help="filename with read mapping"
                             " information. Multiple files can be "
			     " submitted in a comma-separated list"  )

2. Use Experiment.OptionParser instead of optparse.OptionParser. The
   former has some extensions that make creating galaxy xml files
   easier. In particular, Experiment.OptionParser permits supplying
   a list of ','-separated values to options that accept multiple
   values.

3. Follow the CGAT script naming convention. If possible, scripts
   should be named ``<format_in>2<format_out>.py``. Formats can
   be mapped to other types in :doc:`cgat2rdf`. For example,
   ``stats`` and ``table`` are both mapped to the format ``tabular``.

   
	
