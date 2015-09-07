.. _developing:

=========================
Contributing to CGAT code
=========================

We encourage everyone who uses parts of the CGAT code collection to
contribute. Contributions can take many forms: bugreports, bugfixes,
new scripts and pipelines, documentation, tests, etc. All
contributions are welcome.

Checklist for new scripts/modules
=================================

Before adding a new scripts to the repository, please check if the
following are true:

1. The script performs a non-trivial task. If a one-line command line
   entry using standard unix commands can give the same effect, avoid
   adding a script to the repository.

2. The script has a clear purpose. Scripts should follow the 
   `unix philosophy <http://en.wikipedia.org/wiki/Unix_philosophy>`_.
   They should concentrate on one task and do it well. Ideally,
   the major input and output can be read from and written to standard
   input and standard output, respectively. 

3. The script follows the naming convention of CGAT scripts. 

4. The scripts follows the :ref:`styleguide`.

5. The script implements the ``-h/--help`` options. Ideally, the
   script has been derived from :file:`scripts/cgat_script_template.py`.

6. The script can be imported. Ideally, it imports without performing
   any actions or writing output.

7. The script is well documented and the documentation has been added
   to the CGAT documentation. There should be an entry in
   :file:`doc/scripts.rst` and a file
   :file:`doc/scripts/newscript.py`.

8. The script has at least one test case added to :file:`tests` - and
   the test works (see :ref:`Testing`).

Building extensions
===================

Using pyximport_, it is (relatively) straight-forward to add optimized
C-code to python scripts and, for example, access pysam internals and
the underlying samtools library. See for example :doc:`scripts/bam2stats`.

To add an extension, the following needs to be in place:

1. The main script (:file:`scripts/bam2stats.py`). The important lines
   in this script are::

      try:
          import pyximport
          pyximport.install()
          import _bam2stats
      except ImportError:
          import CGAT._bam2stats as _bam2stats

   The snippet first attempts to build and import the extension by
   setting up pyximport_ and then importing the cython module as
   :file:`_bam2stats`.  In case this fails, as is the case for an
   installed code, it looks for a pre-built extension (by
   ``setup.py``) in the CGAT pacakge.
 
2. The cython implementation :file:`_bam2stats.pyx`. This script
   imports the pysam API via::

      from csamtools cimport *

   This statement imports, amongst others, :class:`AlignedRead` into
   the namespace. Speed can be gained from declaring variables. For
   example, to efficiently iterate over a file, an
   :class:`AlignedRead` object is declared::

      # loop over samfile
      cdef AlignedRead read
      for read in samfile:
          ...

3. A :file:`pyxbld` providing pyximport_ with build information.
   Required are the locations of the samtools and pysam header
   libraries of a source installation of pysam plus the
   :file:`csamtools.so` shared library. For example::

     def make_ext(modname, pyxfilename):
	 from distutils.extension import Extension
	 import pysam, os
	 dirname = os.path.dirname( pysam.__file__ )[:-len("pysam")]
	 return Extension(name = modname,
			  sources=[pyxfilename],
			  extra_link_args=[ os.path.join( dirname,
			  	 "csamtools.so")],
			  include_dirs =  pysam.get_include(),
			  define_macros = pysam.get_defines() )


If the script :file:`bam2stats.py` is called the first time,
pyximport_ will compile the cython_ extension :file:`_bam2stats.pyx`
and make it available to the script. Compilation requires a working
compiler and cython_ installation.  Each time :file:`_bam2stats.pyx`
is modified, a new compilation will take place.

pyximport_ comes with cython_.

Writing recipes
===============

Recipes are short use cases demonstrating the use of one or more
CGAT utilities to address a specific problem.

Recipes should be written as ipython_ notebooks. The recipe notebooks
are stored in the :file:`recipes` directory in the repository. Each
recipe is within its individual directory.  This minimizes
interference between each document, but also means that currently each
notebook needs a separate notebook server to be developped.

To build all recipes, type::

    cd recipes
    make html
    make clean

This will build html files that are deposited in the docs directory.

The last cleaning up step is important in order to remove large files created
during the notebook execution.

.. note::
   The commands above require the runipy python module. To install,
   type::

       pip install runipy

Data for recipes can be made available in www.cgat.org/downloads/public/cgat/recipes.
Ideally, recipes should make use of publicly available data sets such
as ENCODE.

Attempt to add a plot to the end of a recipe, using
R commands to create the plot within the notebook.

Writing pipelines
=================

Best practice for CGAT pipelines:

1. All non-trivial code should be extracted to modules or scripts.

2. Modules should not access PARAMS dictionary directly, but
   parameters should be passed to the function.

3. Important processing steps where different external tools could
   potentially be employed the design of the module classes should be
   carefully considered to ensure consistent input and output file
   formats for different tools. PipelineMapping provides a good
   example for this.

4. All production pipelines should include tests for consistency which
   can be run automatically.

5. Where appropriate pipelines should include a small test dataset
   with published results for comparison. This dataset can be run on
   each pipeline run and included in the pipeline report where it can
   be used as a pipeline control.

6. Periodic code review meetings where interested parties can agree of
   major changes to production pipelines and associated modules – to
   be arranged as required.

7. The best way to manage pipeline improvements is by individuals
   using pipelines taking responsibility for incremental
   improvement. As best practice fellows should announce plans to
   modify particular pipelines and modules on the CGAT members list to
   avoid duplication of effort. Fellows should log the changes that
   they make in a change log and document both modules and pipelines
   in detail.

8. Add a section with Requirements to all pipeline scripts and tools.
   Only add them in files where the actual dependency arises, see
   :doc:`modules/Requirements`.
   



