#########################################################################
#
#   MRC FGU Computational Genomics Group
#
#
#   Copyright (C) 2016 Andreas Heger
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
##########################################################################

'''
R.py - Intergrate python and R
===========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This module file contains code to allow integration of python and
R:

   - Redefined R class which enables saving R history.

Usage
-----
'''

from rpy2.robjects import R as RBase


class R_with_History(RBase):
    ''' Redefine the RBase class to allow history to be saved'''
    _instance = None

    def __new__(cls):
        c = RBase.__new__(cls)
        cls._instance = c
        c._history = []
        return cls._instance

    def __call__(self, string):
        self._history.append(string)
        RBase.__call__(self, string)

    def loadImage(self, imageFile):
        self["load.image"](file=imageFile)

    def saveImage(self, imageFile):
        self["save.image"](file=imageFile)

    def saveHistory(self, historyFile, append=False):
        ''' save history '''

        filetype = "w"
        if append:
            filetype = "a"

        with open(historyFile, filetype) as outf:
            outf.write("\n".join(self._history) + "\n")
