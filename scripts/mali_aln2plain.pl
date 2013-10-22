
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

