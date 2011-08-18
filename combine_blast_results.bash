#!/bin/bash

####
####
##
## Project Pairsdb
##
## scripts and makefiles for maintaining and updateable BLAST graph
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <andreas.heger@helsinki.fi>
##
## $Id: combine_blast_results.bash 2041 2008-07-25 09:41:23Z andreas $
##
## combine BLAST results from a cluster run. The blast files are assumed
## to be compressed.

file_output=$1
file_error=$2
dir_tmp=$3

## number of fields in output
num_fields=$4

echo "## `date`: started"

echo "## concatenting files in directory ${dir_tmp} and writing to ${file_output}. Errors go to ${file_error}."
echo "## `date` : combining output."
rm -f ${file_error}

# note: the algorithm has to process the filenames in sorted order,
# so that they don't need to be sorted afterwards

for file in `echo ${dir_tmp}*.blast.gz | sort -k1n`; do 
    echo "## `date` : processing $file."

    ## add a line break at the end, so that incomplete lines do not get
    ## joined.
    ## removed, does not work with gzipped files
    ## echo "" >> ${file}

    ## concatenate files. If you see an incomplete line, append it to
    ## the error file with the prefix missing as we might want to 
    ## redo that sequence again.
    gunzip < ${file} |\
    awk -v FS="\t" -v file=${file} -v file_error=${file_error} -v num_fields=${num_fields} \
        '{ if(NF==num_fields) { print } else { if (NF > 0) { \
              printf("# skipped: %s\n", $1) >> file_error; printf("error in file: %s\n", file) >> file_error; print >> file_error; \
          } }}' | gzip >> ${file_output}

done
  
echo "## `date` : finished"

exit 0





