.. _CGATInstallation:

=========================
Installation instructions
=========================

The section below describes how to install the CGAT scripts. Please
note that we can not test our code on all systems and configurations
out there. If something does not work, please try a
:ref:`CleanInstall`.

Quick installation
==================

Pre-install dependencies
------------------------

Installing CGAT should be straight-forward::

   pip install cgat

CGAT depends on numerous other python packages which themselves might require
manual intervention. Please see :ref:`ManualInstallation` for a
step-by-step installation approach.

Initialization
--------------

In order to run pipelines and code directly from the CGAT script
repository, you need to perform the following initializations::

   python setup.py develop --multi-version

This will compile all the extension modules without installing 
anything. To use, add the CGAT directory to ``$PYTHONPATH``
environment variable::

   export PYTHONPATH=$PYTHONPATH:/location/to/cgat

You might also want to run the script::

   python scripts/cgat_build_extensions.py 

to test if all the scripts with associated cython_ code compile
cleanly.

.. _ManualInstallation:

Manual installation
===================

CGAT depends on numerous other python packages not all of which
might install cleanly. Here, we give some more detailed instructions.
Generally we recommend when troubleshooting CGAT installation to do so
within a virtual environment. To create a clean environment, type::

    virtualenv --no-site-packages cgat-python
    source cgat-python/bin/activate

Now, download the list of required packages::

    wget https://raw.github.com/CGATOxford/cgat/master/requires.txt

To install the required basic packages::

    pip install -r requires.txt

Also, bx-python needs to be installed. The current version on pypi is
currently out of date, so to install, do::

    pip install https://bitbucket.org/james_taylor/bx-python/get/tip.tar.bz2

If all of that works, installing the CGAT tools should now be
straight-forward::

    pip install cgat

If you continue having problems with the installation please try the
:ref:`CleanInstall` guide.

Troubleshooting
---------------

Some packages will require additional system-level packages to 
be installed. The following depencies might cause problems:

PyGreSQL
    requires postgres-devel

PyGTK
    not installable via setuptools, install separately.

biopython
    pip occasionally fails for biopython. If so, try installing 
    manually.

.. _GalaxyInstallation:

Installing in Galaxy
====================

CGAT tools can be used through the `galaxy`_ framework. In order
to set up the CGAT tool box in you own galaxy_ instance, use the 
:file:`cgat2rdf.py` script.

The sequence of commands is:

1. Install Galaxy

2. Install CGAT 

3. Run the `cgat2rdf.py` script (see :doc:`scripts/cgat2rdf`) to create an xml file for inclusion into
   galaxy_. For example, to create a wrapper for `bam2stats.py` (see :doc:`scripts/bam2stats`), run,
   where ``cgat-xml`` is the location of tool xml files within galaxy_::

       python <cgat-scripts>cgat2rdf.py --format=galaxy <cgat-scripts>bam2stats.py > <cgat-xml>bam2stats.xml

4. Add an entry to :file:`tool_conf.xml` for the script within the
   galaxy_ distribution::

      <section name="CGAT Tools" id="cgat_tools">
          <tool file="<cgat-xml>/bam2stats.xml" />
      </section>


A list of galaxy compatible scripts is in file :file:`galaxy.list`. This file is part of the
CGAT repository and can be used to create all wrappers in one go::

   cat galaxy.list
   | cgat2rdf.py
        --source-dir=<cgat-scripts>  --input-regex="(.*).py"
	--output-pattern=<galaxy-xml>/%s.xml --format=galaxy

Within galaxy_, CGAT scripts will use samtools_ formatted genomic
sequences, which are located in the ``sam_fa_indexes`` galaxy_ resource.

