"""Library for using suffix arrays.
"""
import os, string, sys

ParamExecutableLookup   = "sary"

def Search( pattern, suffix_arrays ):

    matches = []

    for suffix_array in suffix_arrays:
        statement = "%s %s %s" % (ParamExecutableLookup, pattern, suffix_array)
        matches += map( lambda x: x[:-1], os.popen( statement ).readlines())
        
    return matches
    

