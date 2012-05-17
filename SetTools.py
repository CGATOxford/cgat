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
SetTools.py - Tools for working on sets
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Many of the functions in this module precede the :py:class:`set` datatype
in python.

Code
----

'''

import itertools
import numpy

def combinations( list_of_sets ):
    '''create all combinations of a list of sets

    returns a list of tuples ( set_composition, union, intersection )
    '''

    sr = range(len(list_of_sets))
    results = []
    for l in range(1, len(list_of_sets)):
        for combination in itertools.combinations( sr, l ):
            union = list_of_sets[combination[0]].union( *[list_of_sets[x] for x in combination[1:]])
            inter = list_of_sets[combination[0]].intersection( *[list_of_sets[x] for x in combination[1:]])

            results.append( (combination, union, inter ) )
    return results

def writeSets( outfile, list_of_sets, labels = None ):
    '''output a list of sets as a tab-separated file.

    *labels* is a list of set labels.
    '''
    
    all_ids = list_of_sets[0].union( *list_of_sets[1:] ) 
    if not labels: labels = range( len(list_of_sets))

    outfile.write("id\t%s\n" % "\t".join(map(str,labels)))

    for i in sorted(list(all_ids)):
        outfile.write( "%s\t%s\n" % (i,
                                  "\t".join( map(str, [ [0,1][i in x] for x in list_of_sets ] ) ) ))
    
def unionIntersectionMatrix( list_of_sets ):
    '''build union and intersection of a list of sets.

    return a matrix with the upper diagonal the union and the lower diagonal the intersection.
    '''

    l = len(list_of_sets)
    matrix = numpy.zeros( (l, l), dtype = numpy.int )
    for x in range( l ):
        xx = list_of_sets[x]
        for y in range( x ):
            yy = list_of_sets[y]
            union = xx.union( yy )
            inter = xx.intersection( yy )
            matrix[x][y] = len(union)
            matrix[y][x] = len(inter)

    return matrix

#------------------------------------------------------------------------
def CompareSets( set1, set2):
    """returns the union and the disjoint members of two sets. The sets have
    to be sorted.
    """

    unique1 = []
    unique2 = []
    common = []    

    x1 = 0
    x2 = 0
    while 1:

        if x1 >= len(set1) or x2 >= len(set2): break
        
        if set2[x2] == set1[x1]:
            common.append( set1[x1] )
            x1 += 1
            x2 += 1
            continue

        
        if set1[x1] < set2[x2]:
            unique1.append(set1[x1])
            x1 += 1
            continue
        
        if set1[x1] > set2[x2]:
            unique2.append( set2[x2] )
            x2 += 1
            
    if x2 < len(set2):
        unique2 += set2[x2:]
    elif x1 < len(set1):
        unique1 += set1[x1:]        
    
    return unique1, unique2, common
        
#------------------------------------------------------------------------
def MakeNonRedundant( set ):

    set.sort()

    last = -1
    result = []
    for x in set:
        if x != last:
            result.append(x)
        last = x
    
    return result

## The following is derived from the python recipe http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/410685
## by Zoran Isailovski 

def getAllCombinations( *sets ):

    if not sets: return []
    F = MakeListComprehensionFunction('F', len(sets) )
    return F(*sets)

def MakeListComprehensionFunction (name, nsets):
    """Returns a function applicable to exactly <nsets> sets.
    The returned function has the signature
    F(set0, set1, ..., set<nsets>)
    and returns a list of all element combinations as tuples.
    A set may be any iterable object.
    """
    if nsets <= 0:
        source = 'def %s(): return []\n' % name
    else:
        constructs = [ ('set%d'%i, 'e%d'%i, 'for e%d in set%d'%(i,i))
                       for i in range(nsets) ]
        a, e, f = map(None, *constructs)
        ##e.reverse() # <- reverse ordering of tuple elements if needed
        source = 'def %s%s:\n   return [%s %s]\n' % \
                 (name, _tuplestr(a), _tuplestr(e), ' '.join(f))
    scope = {}
    exec source in scope
    return scope[name]
    
def _tuplestr(t):
    if not t: return '()'
    return '(' + ','.join(t) + ',)'

def xuniqueCombinations(items, n):
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xuniqueCombinations(items[i+1:],n-1):
                yield [items[i]]+cc


