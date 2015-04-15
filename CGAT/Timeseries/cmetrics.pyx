import CGAT.Experiment as E
import numpy as pynp
cimport numpy as np

def consensus_metrics(array):
    '''Cythonised attempt at consensus clustering
    metrics.'''

    cdef int g
    cdef int i
    cdef int j
    cdef int k
    cdef int l
    cdef int n

    idx = [y for y in xrange(array.shape[0])]
    cidx = [p for p in xrange(array.shape[1])]
    n = len(cidx)
    g = len(idx)
    
    # create memory views of empty arrays to fill
    a_array = pynp.zeros((n, n), dtype=pynp.uint64)
    b_array = pynp.zeros((n, n), dtype=pynp.uint64)
    c_array = pynp.zeros((n, n), dtype=pynp.uint64)
    d_array = pynp.zeros((n, n), dtype=pynp.uint64)

    cdef int [:,:] mem_view = array
    cdef unsigned long [:,:] agree = a_array
    cdef unsigned long [:,:] disagree = b_array
    cdef unsigned long [:,:] agree1 = c_array
    cdef unsigned long [:,:] agree2 = d_array
    cdef unsigned long perms = 0

    E.info("Iteratively counting clustering overlap")
    for i from 0 <= i < n:
        for j from 0 <= j < n:
            perms += 1
            for k from 0 <= k < g:
                for l from 0 <= l < g:
                    if mem_view[k][i] == mem_view[l][i] and mem_view[k][j] == mem_view[l][j]:
                        agree[i][j] += 1
                    elif mem_view[k][i] != mem_view[l][i] and mem_view[k][j] == mem_view[l][j]:
                        agree2[i][j] += 1
                    elif mem_view[k][i] == mem_view[l][i] and mem_view[k][j] != mem_view[l][j]:
                        agree1[i][j] += 1
                    else:
                        disagree[i][j] += 1

    E.info("Counting finished: %i clustering combinations counted" % perms)

    return adjusted_rand_index(agree, disagree, agree1, agree2, n)


def adjusted_rand_index(agree_mat, disagree_mat, agree1_mat, agree2_mat, n):
    '''
    calculate the rand and adjusted rand index over combinations of clusterings
    return matrix.
    '''

    # calculate rand index                                                  
    cdef int h
    cdef int t
    cdef float a = 0.0
    cdef float b = 0.0
    cdef float c = 0.0
    cdef float d = 0.0
    cdef float adjrand = 0.0
    cdef float rand = 0.0
    cdef float n_choose_2 = 0.0
    cdef float sum_n = 0.0
    cdef float max_index = 0.0

    adjrand_array = pynp.zeros((n, n), dtype=pynp.float64)
    rand_array = pynp.zeros((n, n), dtype=pynp.float64)  
    cdef np.float64_t [:,:] adjrand_mat = adjrand_array
    cdef np.float64_t [:,:] rand_mat = rand_array

    for h from 0 <= h < n:
        for t from 0 <= t < n:
            a = agree_mat[h][t]
            b = disagree_mat[h][t]
            c = agree1_mat[h][t]
            d = agree2_mat[h][t]

            # the adjusted rand index is of the form:
            # [index - expected index]/[max index - expected index]
            # a = both agree
            # b = both disagree
            # c = cluster1 agree, cluster2 disagree
            # d = cluster1 disagree, cluster2 agree

            sum_n = (a + b + c + d)
            n_choose_2 = (sum_n*(sum_n-1.0))/2.0
            exp = ((a+d)*(a+c))/n_choose_2
            max_index = (((a+d) + (a+c))/2.0)
            adjrand = (a - exp)/(max_index - exp)
            adjrand_mat[h][t] = adjrand

            # the rand index is of the form:
            # (a+b)/(a+b+c+d)
            rand = (a + b)/(a + b + c + d)
            rand_mat[h][t] = rand

    # coerce memory view to numpy array

    return (pynp.asarray(adjrand_array), pynp.asarray(rand_array))


