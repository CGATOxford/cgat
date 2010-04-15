#!/bin/bash

if [ ! -n "$1" ]; then
    echo "USAGE: wiggle2fasta.bash genome quality"
    exit $E_BADARGS
fi

genome=$1
quality=$2

for x in `cut -f 1 ${genome}.idx`; do
    printf ">%s\n" ${x}; 
    hgWiggle -chr=${x} ${quality} |\
    awk 'BEGIN{ p = -1; } \
	  /^variable/ {next;} \
	  !/^#/ { if ($1 - p != 1) { printf("error: nosequential bases %i %i", p, $1); exit;}  \
	   p = $1; printf("%i ",$2); if( ++l > 20 ) { l = 0; printf ("\n") }} \
         END { printf("\n"); } '; 
done
