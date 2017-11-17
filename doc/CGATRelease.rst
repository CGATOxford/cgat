=============
Release Notes
=============

The code collection is under continuous improvement and the 
latest code can always be found in the code repository.
Nevertheless, we occasionally prepare releases. Notes on
each release are below.

Release 0.3.1
=============

* bugfix cgat CLI
* added scripts to help find R and Python dependencies: scripts/cgat_deps_R.sh, scripts/cgat_deps_python.sh
* modified install script to use conda environment files instead of conda meta-packages; https://github.com/CGATOxford/cgat/pull/343
* update test_style.py to reflect new repository structure; https://github.com/CGATOxford/cgat/pull/345
* update test_import.py to reflect new repository structure; https://github.com/CGATOxford/cgat/commit/decc5ec0e93b4d2c05b6fa3793b5a73caf1fa8b9
* added script to perform Gene Set Enrichment Analysis; https://github.com/CGATOxford/cgat/pull/344
* updated scripts to work with Python 3.6; https://github.com/CGATOxford/cgat/pull/346
* bugfix bed2bed script; https://github.com/CGATOxford/cgat/pull/349; https://github.com/CGATOxford/cgat/issues/347
* bugfix GTF.py module file; https://github.com/CGATOxford/cgat/issues/350
* various bugfixes for Python 3 problems when running pipeline_annotations in python 3, plus a modification of the new "assembly report" section to allow genomes without a standard assembly report (e.g. yeast) to be processed; https://github.com/CGATOxford/cgat/pull/348
* refactor CGAT/PipelineGWAS.py to CGAT/GWAS.py; https://github.com/CGATOxford/cgat/pull/354
* removed unused imports with autoflakes; https://github.com/CGATOxford/cgat/pull/355
* removed basestr, use str instead in Python 3; https://github.com/CGATOxford/cgat/pull/359
* added test for csv2db; https://github.com/CGATOxford/cgat/pull/360
* fixed rename_transcripts test for gtf2gtf.py; https://github.com/CGATOxford/cgat/pull/361
* updated installation; https://github.com/CGATOxford/cgat/pull/364 ; https://github.com/CGATOxford/cgat/pull/373
* updated documentation; https://github.com/CGATOxford/cgat/pull/365 ; https://github.com/CGATOxford/cgat/pull/367
* replaced bx.bbi.bigwig_file with pyBigWig ; https://github.com/CGATOxford/cgat/pull/369
* updated Expression.py to iterate over pandas dataframe; https://github.com/CGATOxford/cgat/pull/370
* updated Expression.py to work with DESeq2 > 1.16 ; https://github.com/CGATOxford/cgat/pull/371
* moved bed_vs_bed.py to obsolete ; https://github.com/CGATOxford/cgat/pull/372

Release 0.3.0
=============

First release of the code working in Python 3 only.

Release 0.2.7
=============

This release is compatible with pysam-0.11.1 and is the last
version of the code compatible with Python 2.7. Next releases
will be in Python 3 only.

Release 0.2.6
=============

This release is compatible with pysam-0.10.0 and it is also half
way through the upgrade to Python 3 compatible code.

Release 0.2.5
=============

Minor bugfix release to ensure compatibility with pysam >= 0.8.4

Release 0.2.4
=============

Minor bugfix release to ensure compatibility with pysam 0.8.4

Release 0.2.3
=============

Minor release to ensure compatibility with pysam-0.8.1

Release 0.2.2
=============

We have reviewed the command line options in all of the CGAT
scripts and have changed them to make them more consistent
between tools and more informative on the command line. This
means that existing scripts that call CGAT tools need to be
updated. A table with all options used in scripts and how
they have or have not changed 
`here <https://github.com/CGATOxford/cgat/blob/master/tests/option_list.tsv>`_.

There is a script called `cgat_refactor.py
<https://github.com/CGATOxford/cgat/blob/master/refactor/cgat_refactor.py>`_
to facilicate refactoring the repository. Instructions on how to use
the script are given at the command line help.

There are also several bugfixes and new features.

Release 0.2.1
=============

Bugfix release

* Fixed a variety of bugs
* Increased test coverage
* Updated documentation
* Added paired read counting to gtf2table

Release 0.2
===========

* release for CGAT manuscript - fixed various installation issues

Release 0.1.9
=============

* alignlib incompatibility fixed
* various bugfixes

Release 0.1.8
=============

* OS X compatibility release
* various bugfixes

Contributors
============

The following people have contributed to the CGAT Code collection:

* Andreas Heger
* Antonio Berlanga-Taylor
* Martin Dienstbier
* Nicholas Ilott
* Jethro Johnson
* Katherine Fawcett
* Stephen Sansom
* David Sims
* Ian Sudbery
* Hu Xiaoming
* Lesheng Kong
* Michael Morgan
* Thomas Smith
* Katherine Brown
* Charlotte George
* Adam Cribbs
* Hania Pavlou
* Reshma Nibhani
* Sebastian Luna-Valero

3rd party code
==============

The CGAT code collection has been made possible by the many developers
in the bioinformatics and python community that have made their code
available for sharing. The code collection includes some snippets of
code taken from elsewhere for convenience, most notably:

1. IGV.py from Brent Petersen 
   https://github.com/brentp/bio-playground/blob/master/igv/igv.py

2. Nested containment list from the Pygr project
   http://code.google.com/p/pygr/

3. SVGdraw.py was written by Fedor Baart & Hans de Wit

4. list_overlap.py from Brent Petersen
   https://github.com/brentp/bio-playground/blob/master/utils/list_overlap_p.py

5. Iterators.py from an unknown source.

Licence
=======

The CGAT code is released under the new BSD licence::

    Copyright (c) 2013, Andreas Heger, MRC CGAT

    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

	Redistributions of source code must retain the above copyright
	notice, this list of conditions and the following disclaimer.
	Redistributions in binary form must reproduce the above copyright
	notice, this list of conditions and the following disclaimer in
	the documentation and/or other materials provided with the
	distribution.  Neither the name of the Medical Research Council nor the
	names of its contributors may be used to endorse or promote
	products derived from this software without specific prior written
	permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

