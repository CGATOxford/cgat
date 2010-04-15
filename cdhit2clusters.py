"""convert cd-hit clusters file to tab separated file.
"""

import sys
import Genomics

if __name__ == "__main__":

    map_rep2mem, map_mem2rep = Genomics.ReadClusters( sys.stdin )

    #  print all members
    for k in map_mem2rep.keys():
        print "%s\t%s" % (k, map_mem2rep[k])

    # print all representatives
    for k in map_rep2mem.keys():
        print "%s\t%s" % (k, k)

