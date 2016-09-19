use strict;

open(IN, "<$ARGV[0]") or die "could not find $ARGV[0]\n"; 
my @a = (<IN>); 
close(IN);

my %a = ();
foreach (@a) {
    chop();
    next if (/^\#/);
    /^(\S+)\s*(.*)/;
    $a{$1} = $2;
}

open(IN, "<$ARGV[1]") or die "could not find $ARGV[1]\n"; 

while (<IN>) {
    /^(\S+)\s*(.*)/;
    if (defined $a{$1}) {
	print $1 . "\t" . $a{$1} . "\t" . $2 . "\n";
    }
}

close(IN);
