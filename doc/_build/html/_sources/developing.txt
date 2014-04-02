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

1. The main script (:file:`scripts/bam2stats.py`). The important lines in this script
   are::

      try:
          import pyximport
          pyximport.install()
          import _bam2stats
      except ImportError:
          import CGAT._bam2stats as _bam2stats

   The snippet first attempts to build and import the extension by
   setting up pyximport_ and then importing the cython module
   as :file:`_bam2stats`. 
   In case this fails, as is the case for an installed code, it 
   looks for a pre-built extension (by ``setup.py``) in the CGAT 
   pacakge.
 
2. The cython implementation :file:`_bam2stats.pyx`. This script imports the pysam API via::

      from csamtools cimport *

   This statement imports, amongst others, :class:`AlignedRead` into the namespace. Speed can be
   gained from declaring variables. For example, to efficiently iterate
   over a file, an :class:`AlignedRead` object is declared::

      # loop over samfile
      cdef AlignedRead read
      for read in samfile:
          ...

3. A :file:`pyxbld` providing pyximport_ with build information. 
   Required are the locations of the samtools and pysam header libraries 
   of a source installation of pysam plus the :file:`csamtools.so` 
   shared library. For example::

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


If the script :file:`bam2stats.py` is called the first time, pyximport_ will 
compile the cython_ extension :file:`_bam2stats.pyx` and make it available 
to the script. Compilation requires a working compiler and cython_ installation.
Each time :file:`_bam2stats.pyx` is modified, a new compilation will take place.

pyximport_ comes with cython_.



