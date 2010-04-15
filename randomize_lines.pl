#!/usr/bin/env perl
# $Id$
#
# randomize_lines.pl < IN > OUT
#
# randomize lines from STDIN
use strict;

my (@lines) = <STDIN>;

srand();
fisher_yates_shuffle( \@lines );    # permutes @array in place


for (@lines) { print $_; };

# generate a random permutation of @array in place
# important: random number is between 1 and i. Choosing
# integer between 1 and N would lead to bias!
sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}
