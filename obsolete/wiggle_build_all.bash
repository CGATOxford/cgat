#!/bin/bash
#
# USAGE:
#
# wiggle_build_all.bash src_dir
#
# index all wiggle tracks in a src directory and save in current directory
# This script runs on the cluster.
#
src_dir=$1

for x in ${src_dir}/*.gz; do 
    aa=`basename $x`; 
    a=${aa%.gz};
    echo -e "gunzip < ${x} | bzip2> ${a}.bz2" > tmp_${a}.qsub; 
    echo -e "bzip-table < ${a}.bz2 > ${a}.bz2t" >> tmp_${a}.qsub;
    echo -e "python /home/andreas/t/wiggle_build_index.py ${a}.bz2" >> tmp_${a}.qsub;
    qsub -cwd -q medium_jobs.q -e index.err -o index.log tmp_${a}.qsub;
    rm -f tmp_${a}.qsub;
done 
