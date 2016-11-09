use strict;

open(IN, "<$ARGV[0]") or die "could not find $ARGV[0]\n"; 

my @a = (<IN>); 
close(IN);

my %a = ();
foreach (@a) {
    chop();
    /^(\S+)\s*(.*)/;
    $a{$1} = $2;
    print $_ . "\n";
}

open(IN, "<$ARGV[1]") or die "could not find $ARGV[1]\n"; 

while (<IN>) {
    /^(\S+)\s*(.*)/;
    print $_ unless (defined $a{$1});
}

close(IN);
