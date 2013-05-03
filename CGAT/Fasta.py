"""Methods for dealing with fasta files.
"""

import string, os, sys, re

class Fasta:

    def __init__ (self,file):
        self.mFile = file
        self.mLastLine = " "
        
    def FetchOne( self ):
        """returns a tuple (description, sequence).

        Returns (None, None) if no more data is there.
        """
        while self.mLastLine != None:
            if self.mLastLine[0] == ">": break
            self.mLastLine = self.mFile.readline()
        else:
            return (None, None)
        
        description = self.mLastLine[1:-1]
        sequence = ""
        while 1:
            line = self.mFile.readline()
            if not line: break
            if line[0] == ">":
                self.mLastLine = line
                break

            sequence += line[:-1]

        return (description, re.sub("\s", "", sequence))
            
            
    
