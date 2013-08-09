#####################################################
# classes and functions for dealing with the output
# from nucmer
#####################################################

import os, sys, re


class Coords(object):
    '''
    class for reading, writing etc the output from 
    show-coords - with nucmer
    '''
    def __init__(self):
        '''
        initialise class
        '''
        self.fields = []

    def __call__(self, coordsfile):
        '''
        geric parsing of the file -  skips evry line until data
        starts. User needs to know columns
        '''
        inputfiles = coordsfile.readline()
        program = coordsfile.readline()
        spacer = coordsfile.readline()
        self.fields = coordsfile.readline().split("\t")
        for line in coordsfile.readlines():
            yield line[:-1].split("\t")


        

