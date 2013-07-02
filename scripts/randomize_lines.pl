#!/usr/bin/env perl
# $Id$
#
# randomize_lines.pl < IN > OUT
#
# randomize lines from STDIN
#
# If the -h option (header) is given, the first line is kept in place.
# The default is to randomize all lines.
#
use strict;
use Getopt::Std;

my %opts = ();
getopts('h', \%opts);

my $param_keep_header = $opts{'h'} || 0;

my (@lines) = <STDIN>;

# output header
if ($param_keep_header)
{ 
    print $lines[0];
    delete $lines[0];
}

# randomize
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
