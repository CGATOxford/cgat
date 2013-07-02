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
PamMatrices.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import string
import re

##-------------------------------------------------------------------------------------------------------------
## data section

# translate pam matrix
new_labels = "ACDEFGHIKLMNPQRSTVWY"
old_labels = "ARNDCQEGHILKMFPSTWYV"

# note: pam-matrices are in the literature cited as Mij, Mij gives the probability that j mutates into i.
# Therefore rows and columns are exchanged here.

pam1 = """
         Ala  Arg  Asn  Asp  Cys  Gln  Glu  Gly  His  Ile  Leu  Lys  Met  Phe  Pro  Ser  Thr  Trp  Tyr  Val
           A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V
Ala A  9867     2    9   10    3    8   17   21    2    6    4    2    6    2   22   35   32    0    2   18
Arg R     1  9913    1    0    1   10    0    0   10    3    1   19    4    1    4    6    1    8    0    1
Asn N     4     1 9822   36    0    4    6    6   21    3    1   13    0    1    2   20    9    1    4    1
Asp D     6     0   42 9859    0    6   53    6    4    1    0    3    0    0    1    5    3    0    0    1
Cys C     1     1    0    0 9973    0    0    0    1    1    0    0    0    0    1    5    1    0    3    2
Gln Q     3     9    4    5    0 9876   27    1   23    1    3    6    4    0    6    2    2    0    0    1
Glu E    10     0    7   56    0   35 9865    4    2    3    1    4    1    0    3    4    2    0    1    2
Gly G    21     1   12   11    1    3    7 9935    1    0    1    2    1    1    3   21    3    0    0    5
His H     1     8   18    3    1   20    1    0 9912    0    1    1    0    2    3    1    1    1    4    1
Ile I     2     2    3    1    2    1    2    0    0 9872    9    2   12    7    0    1    7    0    1   33
Leu L     3     1    3    0    0    6    1    1    4   22 9947    2   45   13    3    1    3    4    2   15
Lys K     2    37   25    6    0   12    7    2    2    4    1 9926   20    0    3    8   11    0    1    1
Met M     1     1    0    0    0    2    0    0    0    5    8    4 9874    1    0    1    2    0    0    4
Phe F     1     1    1    0    0    0    0    1    2    8    6    0    4 9946    0    2    1    3   28    0
Pro P    13     5    2    1    1    8    3    2    5    1    2    2    1    1 9926   12    4    0    0    2
Ser S    28    11   34    7   11    4    6   16    2    2    1    7    4    3   17 9840   38    5    2    2
Thr T    22     2   13    4    1    3    2    2    1   11    2    8    6    1    5   32 9871    0    2    9
Trp W     0     2    0    0    0    0    0    0    0    0    0    0    0    1    0    1    0 9976    1    0
Tyr Y     1     0    3    0    3    0    1    0    4    1    1    0    0   21    0    1    1    2 9945    1
Val V    13     2    1    1    3    2    2    3    3   57   11    1   17    1    3    2   10    0    2 9901
"""

pam250 = """
         Ala  Arg  Asn  Asp  Cys  Gln  Glu  Gly  His  Ile  Leu  Lys  Met  Phe  Pro  Ser  Thr  Trp  Tyr  Val
          A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V
Ala A    13    6    9    9    5    8    9   12    6    8    6    7    7    4   11   11   11    2    4    9
Arg R     3   17    4    3    2    5    3    2    6    3    2    9    4    1    4    4    3    7    2    2
Asn N     4    4    6    7    2    5    6    4    6    3    2    5    3    2    4    5    4    2    3    3
Asp D     5    4    8   11    1    7   10    5    6    3    2    5    3    1    4    5    5    1    2    3
Cys C     2    1    1    1   52    1    1    2    2    2    1    1    1    1    2    3    2    1    4    2
Gln Q     3    5    5    6    1   10    7    3    7    2    3    5    3    1    4    3    3    1    2    3
Glu E     5    4    7   11    1    9   12    5    6    3    2    5    3    1    4    5    5    1    2    3
Gly G    12    5   10   10    4    7    9   27    5    5    4    6    5    3    8   11    9    2    3    7
His H     2    5    5    4    2    7    4    2   15    2    2    3    2    2    3    3    2    2    3    2
Ile I     3    2    2    2    2    2    2    2    2   10    6    2    6    5    2    3    4    1    3    9
Leu L     6    4    4    3    2    6    4    3    5   15   34    4   20   13    5    4    6    6    7   13
Lys K     6   18   10    8    2   10    8    5    8    5    4   24    9    2    6    8    8    4    3    5
Met M     1    1    1    1    0    1    1    1    1    2    3    2    6    2    1    1    1    1    1    2
Phe F     2    1    2    1    1    1    1    1    3    5    6    1    4   32    1    2    2    4   20    3
Pro P     7    5    5    4    3    5    4    5    5    3    3    4    3    2   20    6    5    1    2    4
Ser S     9    6    8    7    7    6    7    9    6    5    4    7    5    3    9   10    9    4    4    6
Thr T     8    5    6    6    4    5    5    6    4    6    4    6    5    3    6    8   11    2    3    6
Trp W     0    2    0    0    0    0    0    0    1    0    1    0    0    1    0    1    0   55    1    0
Tyr Y     1    1    2    1    3    1    1    1    3    2    2    1    2   15    1    2    2    3   31    2
Val V     7    4    4    4    4    4    4    4    5    4   15   10    4   10    5    5    5   72    4   17
"""


##-------------------------------------------------------------------------------------------------------------
def ExponentiateMatrix( matrix ):

    l = len(matrix)
    new_matrix = []
    for row in range(0, l):
        new_matrix.append([0.0] * l)

        for col in range(0,l):
            v = 0
            for x in range(0,l):
                v = v + matrix[row][x] * matrix[col][x]
            new_matrix[row][col] = v
    return new_matrix

def NormalizeMatrix( matrix ):

    l = len(matrix)    
    # normalize matrix:
    for row in range(0,l):
        total = 0.0
        for col in range(0,l):
            total += matrix[row][col]
        for col in range(0,l):
            matrix[row][col] = float(matrix[row][col]) / total


def MultiplyMatrices( matrix_a, matrix_b ):

    l = len(matrix_b)
    m = len(matrix_a[0])
    new_matrix = []
    for row in range(0, l):
        new_matrix.append([0.0] * l)

        for col in range(0,m):
            v = 0
            for x in range(0,m):
                v = v + matrix_a[row][x] * matrix_b[x][col]
            new_matrix[row][col] = v
    return new_matrix

def PrintMatrix( matrix ):
    
    # print matrix
    for row in range(0,len(matrix)):
        print new_labels[row], string.join(map(lambda x: "%6.4f " % x,matrix[row]), "")


def CreateMatrix( matrix_string ):

    new_indices = {}
    for x in range(0,len(new_labels)):
        new_indices[new_labels[x]] = x

    # create matrix
    matrix = []
    for i in range(0,len(new_labels)):
        matrix.append([0.0] * len(new_labels))

    rx = re.compile("\s+") 

    row = 0
    for line in string.split(pam1, "\n")[3:]:

        data = map(string.atoi, rx.split(line)[2:])

        if len(data) <> len(old_labels):
            continue

        for col in range(0,len(old_labels)):
            #exchange row and col
            matrix[new_indices[old_labels[col]]][new_indices[old_labels[row]]] = data[col]

        row = row + 1    

    return matrix

def GetIdentityMatrix():

    matrix = []
    for i in range(0,len(new_labels)):
        matrix.append([0.0] * len(new_labels))

    for i in range(0,len(new_labels)):
        matrix[i][i] = 1.0
        
    return matrix

def GetMatrix( matrix_number ):

    current_matrix = CreateMatrix( pam1 )
    NormalizeMatrix( current_matrix )

    result_matrix = GetIdentityMatrix()

    x = matrix_number
    
    while x > 0:
        if x & 1:
            result_matrix = MultiplyMatrices(result_matrix, current_matrix)
        x = x >> 1
        current_matrix = ExponentiateMatrix(current_matrix)

    return result_matrix
##-------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    
    matrix = GetMatrix( 100 )
    PrintMatrix( matrix )







