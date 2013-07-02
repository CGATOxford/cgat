#!/home/belgardt/bin/python2.6

"""
return_insert_sizes
--------------------

This script returns the genomic insert size distances between mate-paired reads mapping to the same chromosome.  Hence, these values could in theory be negative if the fragment was so short that sequenced ends overlap.  It assumes it is piped an Illumina .export file::

   zcat whatever.export.gz | $RMAAPATH/map/return_insert_sizes > file.for.further.manipulation

The export file is assumed to be organized with mate pairs next to each other as follows:

#. Read 1234/1
#. Read 1234/2
#. Read 1635/1
#. Read 1635/2
#. etc.

Inner-distances are returned on individual lines as follows:

#. 101
#. 120
#. 115
#. etc.
"""

import sys, collections

# read through a piped export file two lines at a time; if both map to the same chromosome, return the absolute value of the difference between the mapped positions
# print this number - ultimately, list will be piped to uniq -c | sort -n | tail -n 1

def main():

   infile = sys.stdin
   outfile = sys.stdout

   line=infile.readline()
   line2=infile.readline()
   seq_len=len(line.split("\t")[8])

   counts = collections.defaultdict( int )

   while 1:

      data = line.split("\t")
      data2 = line2.split("\t")

      try:
         if data[10]==data[10] and 'chr' in data[10]:
            size = abs(int(data[12])-int(data2[12]))-seq_len
            counts[size] += 1
      except (ValueError, IndexError):
         pass

      line=infile.readline()
      line2=infile.readline()

      if not line or not line2: break

   for key in sorted(counts.keys()):
      outfile.write( "%i\t%i\n" % (key, counts[key]) )
      
if __name__ == "__main__":
   sys.exit( main() )

