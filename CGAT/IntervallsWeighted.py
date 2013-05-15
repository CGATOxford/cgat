####
####
##
## Project PythonTools
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: IntervallsWeighted.py 2784 2009-09-10 11:41:14Z andreas $
##
##
####
####
"""
IntervallsWeigted.py - working with weigted intervals
=====================================================

Work with weighted invervalls. A weighted intervall is
a tuple of the form (from,to,weight).

Funktions in this module take an optional function parameter object fct,
that will give the weight if two intervalls are combined. The default
is to add the weights of intervalls that are combined.

This module is work in progress. Finished are:

CombineIntervallsLarge
RemoveIntervallsSpanning
"""

#----------------------------------------------------------------
def CombineIntervallsLarge( intervalls, fct = lambda x,y: x+y ):
    """combine intervalls. Overlapping intervalls
    are concatenated into larger intervalls.
    """
    if not intervalls:
        return []

    new_intervalls = []
    
    intervalls.sort()
    first_from, last_to, last_weight = intervalls[0]
    
    for this_from, this_to, this_weight in intervalls[1:]:
        if this_from > last_to:
            new_intervalls.append( (first_from, last_to, last_weight ) )
            first_from, last_to, last_weight = this_from, this_to, this_weight
            continue

        if last_to < this_to:
            last_to = this_to

        last_weight = fct(last_weight, this_weight)
    
    new_intervalls.append( ( first_from, last_to, last_weight ))

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
def RemoveIntervallsSpanning( intervalls, fct = lambda x,y: x+y ):
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

        last_from, last_to, last_weight = last_intervalls[0]
        
        for this_from, this_to, this_weight in last_intervalls[1:]:
            # print last_from, last_to, this_from, this_to
            # this is larger:
            if this_from <= last_from and this_to >= last_to:
                continue
            
            # last is larger:
            if last_from <= this_from and last_to >= this_to:
                last_from,last_to,last_weight = this_from,this_to,this_weight
                continue

            # no complete overlap
            new_intervalls.append( (last_from, last_to, last_weight ) )
            last_from,last_to,last_weight = this_from,this_to,this_weight        
        
        new_intervalls.append( ( last_from, last_to, last_weight ))

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


