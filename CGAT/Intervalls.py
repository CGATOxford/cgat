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
Intervalls.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
def CombineIntervallsLarge( intervalls ):
    """combine intervalls. Overlapping intervalls
    are concatenated into larger intervalls.
    """
    if not intervalls:
        return []

    new_intervalls = []
    
    intervalls.sort()
    first_from, last_to = intervalls[0]
    
    for this_from, this_to in intervalls[1:]:
        if this_from > last_to:
            new_intervalls.append( (first_from, last_to ) )
            first_from, last_to = this_from, this_to
            continue

        if last_to < this_to:
            last_to = this_to
    
    new_intervalls.append( ( first_from, last_to ))

    return new_intervalls

#----------------------------------------------------------------
def ComplementIntervalls( intervalls, first = None, last = None):
    """complement a list of intervalls with intervalls not
    in list.
    """

    if not intervalls:
        if first and last:
            return [(first,last)]
        else:
            return []
        
    new_intervalls = []
    
    intervalls.sort()
    last_from, last_to = intervalls[0]

    if first != None and first < last_from:
        new_intervalls.append( (first, last_from-1) )

    for this_from, this_to in intervalls:
        if this_from > last_to+1:
            new_intervalls.append( (last_to+1, this_from-1 ) )            

        last_from = this_from
        last_to = max(last_to, this_to)

    if last and last > last_to+1:
        new_intervalls.append( (last_to+1, last))

    return new_intervalls

#----------------------------------------------------------------
def AddComplementIntervalls( intervalls, first = None, last = None):
    """complement a list of intervalls with intervalls not
    in list and return both.
    """

    return intervalls + ComplementIntervalls( intervalls, first, last)

#----------------------------------------
def CombineIntervallsDistance( intervalls, min_distance ):
    """combine a list of non-overlapping intervalls,
    and merge those that are less than a certain
    distance apart.
    """
    
    if not intervalls: return []

    new_intervalls = []
    
    intervalls.sort()
    
    first_from, last_to = intervalls[0]
    
    for this_from, this_to in intervalls[1:]:
        
        if this_from - last_to - 1 >= min_distance:
            new_intervalls.append( (first_from, last_to ) )
            first_from = this_from

        last_to = this_to
    
    new_intervalls.append( ( first_from, last_to ))

    return new_intervalls

#----------------------------------------
def DeleteSmallIntervalls( intervalls, min_length ):
    """combine a list of non-overlapping intervalls,
    and delete those that are too small.
    """
    
    if not intervalls: return []

    new_intervalls = []

    for this_from, this_to in intervalls:
        if (this_to - this_from + 1) >= min_length:
            new_intervalls.append( (this_from, this_to ) )
            
    return new_intervalls

#----------------------------------------------------------------
def CombineIntervallsOverlap( intervalls ):
    """combine intervalls.
    Overlapping intervalls are reduced to their intersection.

    first_from, last_to contain region of current maximum overlapping segment.
    max_right is maximum extension of any sequence overlapping with current
    overlapping segment.
    
    """
    if not intervalls:
        return []

    new_intervalls = []
    
    intervalls.sort()
    
    biggest_from, smallest_to = intervalls[0]
    biggest_to = smallest_to

    tos = []

    for this_from, this_to in intervalls[1:]:

        # no overlap: write everything and reset
        if this_from > biggest_to:
            
            new_intervalls.append( (biggest_from, smallest_to ) )
##             print "-->", biggest_from, smallest_to
            biggest_from = this_from
            smallest_to = this_to
            biggest_to = this_to
##             print this_from, this_to, "biggest_from=", biggest_from, "smallest_to=", smallest_to, "biggest_to=", biggest_to
            tos = []
            continue

        # no overlap with common region: write and reset overlap to new start
        elif this_from > smallest_to:
            new_intervalls.append( (biggest_from, smallest_to ) )
##             print "-->", biggest_from, smallest_to
            biggest_from = this_from
            tos = filter( lambda x, y = biggest_from: x > y, tos )
            if tos:
                smallest_to = min(tos)
            else:
                smallest_to = biggest_to
            
        biggest_to = max( biggest_to, this_to )
        biggest_from = max( biggest_from, this_from)
        smallest_to = min( smallest_to, this_to)

        tos.append( this_to )
        
##         print this_from, this_to, "biggest_from=", biggest_from, "smallest_to=", smallest_to, "biggest_to=", biggest_to, tos

    new_intervalls.append( ( biggest_from, smallest_to ))

    return new_intervalls

#----------------------------------------------------------------
def RemoveIntervallsContained( intervalls ):
    """
    remove intervalls that are fully contained in another.

    [(10, 100), (20, 50), (70, 120), (130, 200), (10, 50), (140, 210), (150, 200)]

    results:
    
    [(10, 100), (70, 120), (130, 200), (140, 210)]   
    """
    if not intervalls:
        return []

    new_intervalls = []
    
    intervalls.sort()
    last_from, last_to = intervalls[0]
    
    for this_from, this_to in intervalls[1:]:
        # this is larger:
        if this_from <= last_from and this_to >= last_to:
            last_from,last_to = this_from,this_to
            continue

        # last is larger
        if last_from <= this_from and last_to >= this_to:
            continue

        # no complete overlap
        new_intervalls.append( (last_from, last_to ) )

        last_from,last_to = this_from,this_to        
        
    new_intervalls.append( ( last_from, last_to ))
    
    return new_intervalls

#----------------------------------------------------------------
def RemoveIntervallsSpanning( intervalls ):
    """remove intervalls that are full covering
    another, i.e. always keep the smallest.

    [(10, 100), (20, 50), (70, 120), (40,80), (130, 200), (10, 50), (140, 210), (150, 200)]

    result:
    
    [(20, 50), (40, 80), (70, 120), (150, 200)]
    """
    
    if not intervalls: return []
    
    intervalls.sort()
    
    last_intervalls = intervalls
        
    while 1:

        new_intervalls = []
    
        last_from, last_to = last_intervalls[0]
    
        for this_from, this_to in last_intervalls[1:]:
            # print last_from, last_to, this_from, this_to
            # this is larger:
            if this_from <= last_from and this_to >= last_to:
                continue
            
            # last is larger:
            if last_from <= this_from and last_to >= this_to:
                last_from,last_to = this_from,this_to
                continue

            # no complete overlap
            new_intervalls.append( (last_from, last_to ) )
            last_from,last_to = this_from,this_to        
        
        new_intervalls.append( ( last_from, last_to ))

        if len(last_intervalls) == len(new_intervalls): break

        last_intervalls = new_intervalls
    
    return new_intervalls

#----------------------------------------------------------------
def ShortenIntervallsOverlap( intervalls, to_remove ):
    """shorten intervalls, so that there is no
    overlap with another set of intervalls.

    assumption: intervalls are not overlapping
    
    """
    if not intervalls:
        return []

    if not to_remove:
        return interalls
    
    new_intervalls = []
    
    intervalls.sort()
    to_remove.sort()

    current_to_remove = 0
    
    for this_from, this_to in intervalls:

        for remove_from, remove_to in to_remove:
            #print this_from, this_to, remove_from, remove_to
            if remove_to < this_from: continue
            if remove_from > this_to: continue

            if remove_from <= this_from and remove_to >= this_to:
                this_from = remove_to
                break

            if this_from < remove_from:
                new_intervalls.append( (this_from, remove_from) )
                # print "adding", this_from, remove_from
            
            this_from = max(this_from, remove_to)

            if this_to < this_from: break
            
        if this_to > this_from:
            # print "adding", this_from, this_to
            new_intervalls.append( (this_from, this_to) )            


    return new_intervalls

#----------------------------------------------------------------
def CalculateOverlap( intervalls1, intervalls2 ):
    """calculate overlap between intervalls.
    """

    if not intervalls1 or not intervalls2:
        return 0

    intervalls1.sort()
    intervalls2.sort()

    overlap = 0
    x = 0
    y = 0

    while x < len(intervalls1) and y < len(intervalls2):

        xfrom, xto = intervalls1[x]
        yfrom, yto = intervalls2[y]

        if xto < yfrom:
            x += 1
        elif yto < xfrom:
            y += 1
        else:
            overlap += min( xto, yto ) - max(xfrom, yfrom ) + 1
            
            if xto < yto:
                x += 1
            elif yto < yto:
                y += 1
            else:
                x += 1
                y += 1

    return overlap
        
        
    
