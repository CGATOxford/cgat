################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
list_overlap.py - compute overlap between lists
===============================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

taken from https://github.com/brentp/bio-playground/blob/master/utils/list_overlap_p.py
find the probability that as high as `shared_genes` is random
given the number of genes: `A_genes`, `B_genes` drawn from `total_genes`
e.g.:

$ %prog shared_genes total_genes A_genes B_genes

or:

$ %prog 10 30000 345 322
0.0043679470685

gives the probability that, 2 random gene subsets (chosen from 30000 genes)
of length 345 and 322 would share at least 10 genes by chance.


This can also be called with 3 files of gene_names:

$ %prog all_genes.txt A_genes.txt B_genes.txt

A_genes.txt and B_genes.txt are then intersected to get the numbers for
the hypergeometric test. Each file must be a single column containing
the gene-name. The comparison *is* case sensitive.


See: http://www.nslij-genetics.org/wli/pub/ieee-embs06.pdf

Usage
-----

Example::

   python cgat_script_template.py --help

Type::

   python cgat_script_template.py --help

for command line help.

Documentation
-------------

Code
----

'''

import optparse
import sys
import os.path as op
import scipy.stats as ss
import CGAT.Experiment as E

def hypergeom(m, n, n1, n2, p=False):
    """
>>> hypergeom(1, 1000, 1000, 1000) # has to be shared.
1.0

>>> all(hypergeom(i, 1000, 1000, 1000) == 1.0 for i in range(100))
True

>>> hypergeom(1, 30000, 20, 20)
0.013253396616299651

>>> hypergeom(2, 30000, 20, 20)
7.9649366037104485e-05

>>> hypergeom(11, 30000, 20, 20)
4.516176321800458e-11

>>> hypergeom(10, 30000, 20, 20) # very low prob.
4.516176321800458e-11

>>> hypergeom(20, 30000, 20, 20) # very low chance that all are shared.
4.516176321800458e-11

"""
    if m <= 0: return 1.0
    mmin = m - 1
    mmax = min(n1, n2)
    return ss.hypergeom.cdf(mmax, n, n1, n2) - ss.hypergeom.cdf(mmin, n, n1, n2)

def with_genes(fftot, ffa, ffb, asfile=True):
    """
    given 3 genelists, calculate the p-value of the shared
    genes between fa and fb that are drawn from ftot.
    """

    if asfile:
        ftot = frozenset(f.strip() for f in open(fftot) if f.strip())
        fa = frozenset(f.strip() for f in open(ffa) if f.strip())
        fb = frozenset(f.strip() for f in open(ffb) if f.strip())
    else:
        fa, fb, ftot = frozenset(ffa), frozenset(ffb), frozenset(fftot)
    
    n1, n2 = len(fa), len(fb)
    m = len(fa.intersection(fb))
    n = len(ftot)

    if asfile:
        print "A : %-32s:%-5i" % (ffa, n1)
        print "B : %-32s:%-5i" % (ffb, n2)
        print "total : %-32s:%-5i" % (fftot, n)
        print "shared: %-32s:%-5i" % (' ', m)
    else:
        print "A : %-32s:%-5i" % ("set A", n1)
        print "B : %-32s:%-5i" % ("set B", n2)
        print "total : %-32s:%-5i" % ("total", n)
        print "shared: %-32s:%-5i" % (' ', m)
   
    return hypergeom(m, n, n1, n2)


def main():
    p = E.OptionParser(__doc__,
                       version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", )
    opts, args = p.parse_args()
    if (len(args) not in (3, 4)):
        sys.exit(not p.print_help())
    if len(args) == 4 and not all(a.isdigit() for a in args):
        print >>sys.stderr, "four arguments must be integers"
        sys.exit(not p.print_help())
    elif len(args) == 3 and not all(op.exists(f) for f in args):
        sys.exit(not p.print_help())

    if len(args) == 4:
        args = map(long, args)
        m, n, n1, n2 = args
        result = hypergeom(m, n, n1, n2)
        #print type(result)
        print result
    else:
        tot_genes, a_genes, b_genes = map(str.strip, args)
        print with_genes(tot_genes, a_genes, b_genes)



if __name__ == "__main__":
    #import doctest
    #if doctest.testmod(optionflags=doctest.ELLIPSIS).failed == 0:
    #    main()
    main()

