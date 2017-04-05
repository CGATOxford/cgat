'''SetTools.py - Tools for working on sets
==========================================

Some of the functions in this module precede the :py:class:`set`
datatype in python.

Reference
---------

'''

import itertools
import numpy


def combinations(list_of_sets):
    '''create all combinations of a list of sets

    >>> combinations([set((1,2)), set((2,3))])
    [((0,), set([1, 2]), set([1, 2])), ((1,), set([2, 3]), set([2, 3]))]
    >>> combinations([set((1,2)), set((2,3)), set((3,4))])
    [((0,), set([1, 2]), set([1, 2])), ((1,), set([2, 3]), set([2, 3])), ((2,), set([3, 4]), set([3, 4])), ((0, 1), set([1, 2, 3]), set([2])), ((0, 2), set([1, 2, 3, 4]), set([])), ((1, 2), set([2, 3, 4]), set([3]))]

    Returns
    -------
    result : list
        The resut is a list of tuples containing (set_composition, union,
        intersection)
    '''

    sr = list(range(len(list_of_sets)))
    results = []
    for l in range(1, len(list_of_sets)):
        for combination in itertools.combinations(sr, l):
            union = list_of_sets[combination[0]].union(
                *[list_of_sets[x] for x in combination[1:]])
            inter = list_of_sets[combination[0]].intersection(
                *[list_of_sets[x] for x in combination[1:]])

            results.append((combination, union, inter))
    return results


def writeSets(outfile, list_of_sets, labels=None):
    '''output a list of sets as a tab-separated file.

    This method build a list of all items contained across all sets
    and outputs a matrix of 0's and 1's denoting set membership. The
    items are in the table rows and the sets are in the table columns.

    Arguments
    ---------
    outfile : File
        File to write to
    list_of_sets : list
        The list of sets to output
    labels : list
        List of labels(column names)

    '''

    all_ids = list_of_sets[0].union(*list_of_sets[1:])
    if not labels:
        labels = list(range(len(list_of_sets)))

    outfile.write("id\t%s\n" % "\t".join(map(str, labels)))

    for i in sorted(list(all_ids)):
        outfile.write("%s\t%s\n" % (
            i,
            "\t".join(map(str, [[0, 1][i in x] for x in list_of_sets]))))


def unionIntersectionMatrix(list_of_sets):
    '''build union and intersection matrix of a list of sets.

    >>> unionIntersectionMatrix([set((1,2)), set((2,3))])
    array([[0, 1],
           [3, 0]])
    >>> unionIntersectionMatrix([set((1,2)), set((2,3)), set((3,4))])
    array([[0, 1, 0],
           [3, 0, 1],
           [4, 3, 0]])

    Arguments
    ---------
    list_of_sets : list
        The list of sets to work with.

    Returns
    -------
    matrix : numpy.matrix
        The matrix is a list of lists. The upper diagonal of the
        matrix contains the size of the union of two sets and the
        lower diagonal the intersection of two sets.

    '''

    l = len(list_of_sets)
    matrix = numpy.zeros((l, l), dtype=numpy.int)
    for x in range(l):
        xx = list_of_sets[x]
        for y in range(x):
            yy = list_of_sets[y]
            union = xx.union(yy)
            inter = xx.intersection(yy)
            matrix[x][y] = len(union)
            matrix[y][x] = len(inter)

    return matrix


def getAllCombinations(*sets):
    """generate all combination of elements from a collection of sets.

    This method is derived from a python recipe by Zoran Isailovski:
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/410685

    >>> getAllCombinations(set((1,2)), set((2,3)), set((3,4)))
    [(1, 2, 3), (1, 2, 4), (1, 3, 3), (1, 3, 4), (2, 2, 3), (2, 2, 4), (2, 3, 3), (2, 3, 4)]
    """
    if not sets:
        return []
    F = _makeListComprehensionFunction('F', len(sets))
    return F(*sets)


def _makeListComprehensionFunction(name, nsets):
    """Returns a function applicable to exactly <nsets> sets.  The
    returned function has the signature F(set0, set1, ..., set<nsets>)
    and returns a list of all element combinations as tuples.  A set
    may be any iterable object.
    """
    if nsets <= 0:
        source = 'def %s(): return []\n' % name
    else:
        constructs = [('set%d' % i, 'e%d' % i, 'for e%d in set%d' % (i, i))
                      for i in range(nsets)]
        a, e, f = map(None, *constructs)
        # e.reverse() # <- reverse ordering of tuple elements if needed
        source = 'def %s%s:\n   return [%s %s]\n' % \
                 (name, _tuplestr(a), _tuplestr(e), ' '.join(f))
    scope = {}
    exec(source, scope)
    return scope[name]


def _tuplestr(t):
    if not t:
        return '()'
    return '(' + ','.join(t) + ',)'


def xuniqueCombinations(items, n):
    """Return a list of unique combinations of items in list.

    >>> list(xuniqueCombinations([1, 2, 3], 1))
    [[1], [2], [3]]
    >>> list(xuniqueCombinations([1, 2, 3], 2))
    [[1, 2], [1, 3], [2, 3]]
    >>> list(xuniqueCombinations([1, 2, 3], 3))
    [[1, 2, 3]]
    """

    if n == 0:
        yield []
    else:
        for i in range(len(items)):
            for cc in xuniqueCombinations(items[i + 1:], n - 1):
                yield [items[i]] + cc


def compareLists(list1, list2):
    """returns the union and the disjoint members of two lists.

    .. note:: Deprecated
        Use python sets instead.

    Returns
    -------
    unique1 : set
        Elements unique in set1
    unique2 : set
        Elements unique in set2
    common: set
        Elements in both lists.

    """

    unique1 = []
    unique2 = []
    common = []

    set1 = sorted(list1)
    set2 = sorted(list2)

    x1 = 0
    x2 = 0
    while 1:

        if x1 >= len(set1) or x2 >= len(set2):
            break

        if set2[x2] == set1[x1]:
            common.append(set1[x1])
            x1 += 1
            x2 += 1
            continue

        if set1[x1] < set2[x2]:
            unique1.append(set1[x1])
            x1 += 1
            continue

        if set1[x1] > set2[x2]:
            unique2.append(set2[x2])
            x2 += 1

    if x2 < len(set2):
        unique2 += set2[x2:]
    elif x1 < len(set1):
        unique1 += set1[x1:]

    return unique1, unique2, common
