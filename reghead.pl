#!/usr/bin/env perl
# $id$
#
# perl reghead.pl token [flag] < in > out
#
# regular expression head
#
# prints everything from stdin to stdout until regular expression 
# in token has matched, including the line, where the regular 
# expression matched, unless flag is set.

my ($token, $flag) = @ARGV;

while (<STDIN>) {
    last if (/$token/);
    print $_;
}

print $_ unless ($flag);


