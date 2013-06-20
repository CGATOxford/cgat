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
ProgressBar.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import sys, os, string, time

class ProgressBar:
    """progress bar class

    adapted from Randy Pargman (2002)
    see http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/168639
    """
    
    def __init__(self, minValue = 0, maxValue = 10, totalWidth=12):
        
        self.progBar = "[]"   # This holds the progress bar string
        self._old_pbar = ""
        self.min = minValue
        self.max = maxValue
        self.span = maxValue - minValue
        self.width = totalWidth
        self.amount = 0       # When amount == max, we are 100% done
        self.updateAmount(0)  # Build progress bar string

    def updateAmount(self, newAmount = 0):

        if newAmount < self.min: newAmount = self.min
        if newAmount > self.max: newAmount = self.max
        self.amount = newAmount

        # Figure out the new percent done, round to an integer
        diffFromMin = float(self.amount - self.min)
        percentDone = (diffFromMin / float(self.span)) * 100.0
        percentDone = round(percentDone)
        percentDone = int(percentDone)

        # Figure out how many hash bars the percentage should be
        allFull = self.width - 2
        numHashes = (percentDone / 100.0) * allFull
        numHashes = int(round(numHashes))

        # build a progress bar with hashes and spaces
        self.progBar = "[" + '#'*numHashes + ' '*(allFull-numHashes) + "]"

        # figure out where to put the percentage, roughly centered
        percentPlace = (len(self.progBar) / 2) - len(str(percentDone))
        percentString = str(percentDone) + "%"

        # slice the percentage into the bar
        self.progBar = self.progBar[0:percentPlace] + percentString + \
                       self.progBar[percentPlace+len(percentString):]

    def draw(self):
        # draw progress bar - but only if it has changed

        if self.progBar != self._old_pbar:
            self._old_pbar = self.progBar
            sys.stdout.write(self.progBar + '\n')
            sys.stdout.flush()

    def __str__(self):
        return str(self.progBar)


if __name__ == "__main__":

    mi, ma= 0, 100
    bar = ProgressBar( 0, 100, 20 )

    for x in range(mi, ma):
        time.sleep(1)
        bar.updateAmount( x )
        bar.draw()
        
