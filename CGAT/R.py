'''
R.py - Intergrate python and R
===========================================================

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
