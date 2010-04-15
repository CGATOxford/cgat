#!/usr/bin/env perl
# $id$
#
# perl regextract.pl regex1 n1 regex2 n2 [flag] < in > out
#
# regular expression extraction
#
# prints everything from stdin to stdout between two regular
# expression that have matched. Delimiters are printed, unless flag is
# set. N1 and n2 specify, how many of the regex matches are skipped

my ($token1, $skip1, $token2, $skip2, $flag) = @ARGV;

while (<STDIN>) {
    if (/$token1/) {
	$skip1--;
	last if ($skip1 < 0);
    }
}

print $_ unless ($flag);

while (<STDIN>) {
    if (/$token2/) {
	$skip2--;
	if ($skip2 < 0) {
	    print $_ unless ($flag);
	    last;
	}
    }
    print $_;
}

## finish pipe
while (<STDIN>) {};

