"""general functions for working with dot-files.

Version: $Id: Dots.py 2784 2009-09-10 11:41:14Z andreas $
Author:  Andreas Heger (2003)
"""

import string
import alignlib

def ReadFromFile( file, min_dot_score = 0 ):
    """read dots from a file.
    """
    
    dots_list = file.readlines()
    dots = alignlib.makeAlignataMatrixRow()

    max_row = 0
    max_col = 0

    has_prefix = None
    checked_prefix = None

    for dot in dots_list:
        if dot == "\n": continue
        if dot[0] == "#": continue
        if dot[0] == "d": continue        

        data = string.split(dot, "\t")

        if len(data) < 2: raise "parsing error in line", dot
        
        if not checked_prefix:
            if string.find(data[0], "-"): has_prefix = 1
            checked_prefix = 1

        row, col = data[:2]
        
        if has_prefix:
            row=string.split(row, "-")[1]
            col=string.split(col, "-")[1]            

        row = string.atoi(row)
        col = string.atoi(col)
            

        if len(data) >= 3:
            score = string.atof( data[2] )
            if score < min_dot_score:
                continue
        else:
            score = 0
            

        max_row = max( max_row, row)
        max_col = max( max_col, col)            

        dots.addPairExplicit( row, col, score)

    return max_row, max_col, dots

