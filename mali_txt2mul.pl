#!/usr/bin/env perl
#
# $Id$
#
# txt2mul.pl < in > out
# 
# convert txt formatted malis (from bioedit)
# to mul formated malis

while(<STDIN>) {
    last if (/\.\.\.\.\|\.\./);
}
<STDIN>;
while(<STDIN>) {
    chop();
    s/~/./g;
    my ($h, $ali) = /(\S+\s+)\s(.+)/;
    $ali =~ s/\s//g;
    print $h . $ali . "\n";
}
