#!/home/belgardt/bin/python2.6

"""
filter_exports
--------------

Creates a fastq file of mate-paired reads from an export file::

   $RMAAPATH/filter_exports.py /path/to/flowcell.lane.export.00 /path/to/flowcell.lane.00.fq.1 /path/to/flowcell.lane.00.fq.2 1 True

Note, by default this also creates *.fq.2, which includes the second mate pair.  The module also includes a set of functions that can be used for filtering from an export file for other purposes.  Last two argumenta are the number of bases that should be removed from the right (leave as '0' if none) and if the export file uses new quality scores or not (pipeline 1.3+).

Output files are compressed with gzip.

"""

import sys
import gzip

def solexa_char_to_sanger(__char):
   return chr(ord(__char)-31)

def solexa_str_to_sanger(__qualities):
   __new_quals=""
   for __char in __qualities:
      __new_quals="".join([__new_quals,solexa_char_to_sanger(__char)])
   return __new_quals

def return_fastq_entry(__export_line_array, remove_bases, new_quals):
   __name=":".join([__export_line_array[0], __export_line_array[2], __export_line_array[3], __export_line_array[4], __export_line_array[5]])
   __name="/".join([__name, __export_line_array[7]])
   if new_quals:
      return "".join(['@', __name, '\n', "".join(__export_line_array[8][0:len(__export_line_array[8])-remove_bases]), '\n+\n', "".join(__export_line_array[9][0:len(__export_line_array[9])-remove_bases]), '\n'])
   else:
      return  "".join(['@', __name, '\n', "".join(__export_line_array[8][0:len(__export_line_array[8])-remove_bases]), '\n+\n', "".join(solexa_str_to_sanger(__export_line_array[9])[0:len(__export_line_array[9])-remove_bases]), '\n'])

def main( args = sys.argv ):
   exp_file_nm=args[1]
   fq_file_1=gzip.open(args[2], 'w')
   fq_file_2=gzip.open(args[3], 'w')
   exp_file=gzip.open(exp_file_nm,'rb')
   while True:
      try:
         line1=exp_file.readline().rstrip('\n')
         line2=exp_file.readline().rstrip('\n')
         fq_file_1.write(return_fastq_entry(line1.split('\t'), int(args[4]), bool(args[5]))); fq_file_2.write(return_fastq_entry(line2.split('\t'), int(args[4]), bool(args[5])))
      except:
         break
   fq_file_1.close(); fq_file_2.close()
   return

if __name__ == "__main__":
   sys.exit( main (sys.argv) )

