#!/usr/bin/env perl
#
# $Id$
#
# mul2fasta.pl < in > out
# 
# convert mul formated malis
# to fasta formatted malis

while (<STDIN>) {
    /^(\S+)\s+(\S+)/;
    print ">$1\n$2\n";
}
