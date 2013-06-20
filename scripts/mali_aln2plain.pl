#!/usr/bin/env perl
#
# $Id$
#
# aln2mul.pl < in > out
# 
# convert aln formated malis
# to mul formated malis

my %lines;
<STDIN>;
while (<STDIN>) {
    chop();
    if (/^\s*(\S+)\s+(\S+)/) {
	$lines{$1} .= $2;
    }
}

my $max_length = 0;

foreach (sort keys %lines) {
    my $x = $lines{$_};

    $x =~ s/\.//g;
    $x =~ s/\-//g;    

    my $l = length($x);
    print "1\t" . $lines{$_} . "\t$l\t$_\n";
}

