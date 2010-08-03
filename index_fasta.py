################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
'''
index_fasta.py - Index fasta formatted files 
============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script indexes one or more :term:`fasta` formatted files into a 
database that can be used by :mod:`IndexedFasta` for quick acces to
sequences.

By default, the database itself is a :term:`fasta` formatted file
itself. Compression methods are available to conserve disk space,
though they do come at a performance penalty.

See also http://pypi.python.org/pypi/pyfasta for a similar
implementation.

Usage
-----

Example::

   python index_fasta.py oa_ornAna1_softmasked /net/cpp-mirror/ucsc/ornAna1/bigZips/ornAna1.fa.gz > oa_ornAna1_softmasked.log

Type::

   python index_fasta.py --help

for command line help.

Documentation
-------------

Code
----

'''
import IndexedFasta

if __name__ == "__main__":
    IndexedFasta.main()
