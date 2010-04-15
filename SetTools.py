####
####
##
## Project PythonTools
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: SetTools.py 2784 2009-09-10 11:41:14Z andreas $
##
##
####
####


import array
import numpy

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


