#!/bin/bash
#
# Execute md5sum on file given as argument.
# 
# If file ends with .gz, uncompress and compute
# from stream.
#
# In both cases, ignore any lines starting with "#"
#
# The md5 sum is prefixed by the file name
#
suffix="${1##*.}"

printf "%s\t" $1

if [[ $suffix = 'gz' ]]; then
    zcat $1 | grep -v '^#' | md5sum
else
    cat $1 | grep -v '^#' | md5sum
fi
