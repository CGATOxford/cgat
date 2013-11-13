
my %lines;
<STDIN>;
while (<STDIN>) {
    chop();
    if (/^(\S+)\s+(\S+)/) {
	$lines{$1} .= $2;
    }
}

my $max_length = 0;

foreach (sort keys %lines) {
    $max_length = length($_) if (length($_) > $max_length);
}

foreach (sort keys %lines) {
    $lines{$_} =~ s/\-/\./g;
    printf "%-${max_length}s  %s\n",  $_ , $lines{$_};
}
