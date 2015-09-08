.. _testing:

=======
Testing
=======

This module describes the implementation of tests for the CGAT
code collection. The CGAT testing includes 

* Regression testing of CGAT scripts
* Testing CGAT code for style conformance to pep8
* Testing CGAT code to be importable
* Unit testing of CGAT modules

Regression testing of scripts
=============================

Scripts are regression tested by comparing the expected output with
the latest output. The tests are implemented in the script
:file:`test_scripts.py`.

This script collects tests from subdirectories in the :file:`tests`
directory. Each test is named by the name of the script it tests.

Adding a new test manually
--------------------------

To add a new test for a CGAT script, create a new :term:`test
directory` in the directory :file:`tests`. The name of the :term:`test
directory` has to correspond to the name of the script the tests will
tested.

In this directory, create a file called :file:`tests.yaml`. This file is
in :term:`yaml` format, a simple text-based format to describe nested data
structures.

The :file:`tests.yaml` file contains the descriptions of the
individual tests to run. Each test is a separate data structure in
this file. The fields are:

options
	Command line options for running the test. If you need to
	provide additional files as input, use the ``%DIR%`` place
	holder for the :term:`test directory`.

stdin
	Filename of file to use as stdin to the script. If no stdin is
	required, set to ``null`` or omit. 

outputs
	A list of output files obtained by running the script that
	should be compared to the list of files in ``references``.
	``stdout`` signifies the standard output.

references
	A list of expected output files. The order of ``outputs`` and
	``references`` should be the same. The reference files are
	expected to be found in the directory :term:`test directory`
	and thus need not prefixed with a directory place holder.

description
	A description of test.

To illustrate, we will be creating tests for the scripts
``fasta2counts.py``. First we create the :term:`test directory`
:file:`tests/fasta2counts.py`. Next we create a file
:file:`tests/fasta2counts.py/tests.yaml` with the following content::

   basic_test:
       outputs: [stdout]
       stdin: null 
       references: [test1.tsv]
       options: --genome-file=<DIR>/small_genome

``basic_test`` is the name of the test. There is no standard input
and the output of the script goes to stdout. Stdout will be compared to
the file :file:`test1.tsv`. The script requires the ``--genome-file``
option, which we supply in the ``options`` field. The ``<DIR>`` prefix
will be expanded to the directory that contains the file
:file:`tests.yaml`.

Finally, we create the required input and reference files in the
:term:`test directory`. Our directory structure looks thus::

   |___tests
     |___fasta2counts.py
     | |___small_genome.fasta
     | |___small_genome.idx
     | |___test1.tsv
     | |___tests.yaml

Multiple tests per script can be defined by adding additional data structures in
the :file:`tests.yaml` file.

Please write abundant tests, but keep test data to a minimum. Thus,
instead of running on a large bam file, create stripped down versions
containing only relevant data that is sufficient for the test at hand.

Re-use test data as much as possible. Some
generic test data used by multiple tests is in the :file:`tests/data`
directory. 

Creating a test
---------------

The script :file:`tests/setup_test.py` can be used to set up 
a testing stub. For example::

   python tests/setup_test.py scripts/bam2bam.py

will add a new test for the script :file:`bam2bam.py`.

The script will create a new testing directory for each script passed
on the command line and create a simple :file:`tests.yaml` file. The
basic test will simply call a script to check if starts without error
and returns a version string.

Running tests
-------------

The CGAT code collection runs both under nose or py.test. In order to
run the tests on CGAT scripts under the nose framework, type::

   nosetests tests/test_scripts.py

In order to get more information, type::

   nosetests -v tests/test_scripts.py

To run individual tests, edit the file
:file:`tests/test_scripts.yaml`. In order to restrict testing to
a single script, for example ``beds2counts.py``, add the following::

   restrict:
         regex: beds2counts.py

To run the tests using py.test, type::

   py.test tests/test_scripts.py

Testing for style
=================

All of CGAT python code are tested for pep8_ conformance using the
pep8 tools. Not all pep8 rules are enforced, though we aim for
increasing compatibility with pep8. Please see also the
:ref:`styleguide`.

The testing is controlled by the script :file:`test_style.py`.
In order to run the tests, type::

   nosetests tests/test_scripts.py

We have also added a test that will scan all the command line options
used in CGAT script against a white- and black-list of
acceptable/unacceptable option names. The purpose of this test is to
ensure consistency between scripts. To run this test, type::

   nosetests tests/test_commandline.py

This test is based on a list of acceptable/unacceptable options in
:file:`tests/option_list.tsv` that is within the repository. The list
has been created by the script :file:`cgat_get_option_list.py` and
been manually annotated. Errors are flagged if a deprecated option
is used in a script or an unregistered option is encountered. To
update option list, type::

   python scripts/cgat_get_option_list.py --in-place --options-tsv-file=tests/option_list.tsv 

Testing for import
==================

In order for documentation to be built or scripts to be 
usable by the ``cgat`` frontend, scripts need to be importable
from anywhere. Importability might fail if a script or module
executes statements on improt or rely reading from input or
configuration files that are not present or have non-sensical
values.

The testing is controlled by the script :file:`test_style.py`.
In order to run the tests, type::

   nosetests tests/test_import.py

Testing modules
===============

There are some unit tests for specific functions in modules, but
a testing regime has not been formalized.

.. _pep8: http://legacy.python.org/dev/peps/pep-0008/

Code coverage
=============

Code coverage for modules can be computed if the
`python coverage <https://pypi.python.org/pypi/coverage>`_ module has
been installed. To compute coverage, use::

   nosetests --with-coverage --cover-package=CGAT --cover-package=scripts tests/test_scripts.py

or using py.test::

   py.test -s tests/test_scripts.py --cov=`pwd` >& out

