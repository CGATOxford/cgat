.. _CGATInstallation:

=========================
Installation instructions
=========================

The section below describes how to install the CGAT scripts. We distinguish between two different installation
types: production and development. The former refers to a well tested subset of scripts published here_, and is
the recommended installation. The latter refers to the whole collection of scripts developed at CGAT, which
may contain code under active development.

Please note that we can not test our code on all systems and configurations out there so please bear with us.

Quick installation
==================

Install using Conda
-------------------

The preferred method to install the CGAT code collection is using the installation script, which uses conda_.

Here are the steps::

        # download installation script:
        curl -O https://raw.githubusercontent.com/CGATOxford/cgat/master/install-CGAT-tools.sh

        # see help:
        bash install-CGAT-tools.sh

        # install set of production scripts (well tested):
        bash install-CGAT-tools.sh --production [--location </full/path/to/folder/without/trailing/slash>]

        # or go for the latest development version:
        bash install-CGAT-tools.sh --devel [--location </full/path/to/folder/without/trailing/slash>]

        # enable the conda environment as requested by the installation script:
        source </full/path/to/folder/without/trailing/slash>/conda-install/bin/activate cgat-s

        # finally, please run the cgatflow command-line tool to check the installation:
        cgat --help

The installation script will put everything under the specified location. It needs 5 GB of disk space
and it takes about 10 minutes to complete. The aim of the script is to provide a portable installation
that does not interfere with the existing software. As a result, you will have a conda environment
working with the CGAT scripts which can be enabled on demand according to your needs.


Install using pip
-----------------

You can also use pip_ to install the CGAT scripts. To go down this route, please type::

   pip install cgat

However, CGAT depends on numerous other python packages which themselves might require
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

The CGAT installation requires setuptools version 1.1 or higher
to be installed. If your system has no setuptools installed, or
an old version, please install setuptools_ first by::

   wget https://bitbucket.org/pypa/setuptools/raw/bootstrap/ez_setup.py -O - | python

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

.. note::
   The order in which packages are installed matters. The order	
   in :file:`requires.txt` should work, but pip might ignore that. To
   install requirements in order, try the following::
      
       while read line ; do echo $line > x ; pip install -r x; rm x; done < requires.txt

If all of that works, installing the CGAT tools should now be
straight-forward::

    pip install cgat

Troubleshooting
---------------

Some packages will require additional system-level packages to 
be installed. The following depencies might cause problems:

PyGreSQL
    requires postgres-devel

PyGTK
    not installable via setuptools_, install separately.

biopython_
    pip occasionally fails for biopython_. If so, try installing 
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

3. Run the `cgat2rdf.py` script (see :doc:`scripts/cgat2rdf`) to
   create an xml file for inclusion into galaxy_. For example, to
   create a wrapper for `bam2stats.py` (see :doc:`scripts/bam2stats`),
   run, where ``cgat-xml`` is the location of tool xml files within
   galaxy_::

       python <cgat-scripts>cgat2rdf.py --format=galaxy <cgat-scripts>bam2stats.py > <cgat-xml>bam2stats.xml

4. Add an entry to :file:`tool_conf.xml` for the script within the
   galaxy_ distribution::

      <section name="CGAT Tools" id="cgat_tools">
          <tool file="<cgat-xml>/bam2stats.xml" />
      </section>


A list of galaxy compatible scripts is in file
:file:`galaxy.list`. This file is part of the CGAT repository and can
be used to create all wrappers in one go::

   cat galaxy.list
   | cgat2rdf.py
        --source-dir=<cgat-scripts>  --input-regex="(.*).py"
	--output-filename-pattern=<galaxy-xml>/%s.xml --format=galaxy

Within galaxy_, CGAT scripts will use samtools_ formatted genomic
sequences, which are located in the ``sam_fa_indexes`` galaxy_
resource.

.. _setuptools: https://pypi.python.org/pypi/setuptools
.. _biopython: http://biopython.org/
.. _conda: https://conda.io
.. _pip: https://pypi.python.org/pypi/CGAT
.. _here: https://doi.org/10.1093/bioinformatics/btt756
