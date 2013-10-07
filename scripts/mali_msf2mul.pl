use strict;

while(<STDIN>) {
    last if (/\/\//);
}
<STDIN>;

my %lines;
my $max_length = 0;
while(<STDIN>) {
    chop();
    next if (!$_);
    my ($i, $m) = /^(\S+)\s+(.+)$/;
    $lines{$i}.=$m;
    $max_length = length($i) if (length($i) > $max_length);
}

for (keys %lines) {
    my $a = $lines{$_};
    $a =~ s/\s//g;
    printf "%-${max_length}s  %s\n", $_ , $a;
}
