#!/usr/bin/env perl
# $id$
#
# perl regtail.pl token [flag] < in > out
#
# regular expression tail
#
# prints everything from stdin to stdout after regular expression 
# in token has matched, including the line, where the regular 
# expression matched, unless flag is set.

my ($token, $flag) = @ARGV;

while (<STDIN>) {
    last if (/$token/);
}

print $_ unless ($flag);

while (<STDIN>) {
    print $_;
}
