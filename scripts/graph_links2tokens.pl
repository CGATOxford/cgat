use strict;

my %tokens;

while (<STDIN>) {
    chop();
    my ($token1, $token2) = /^(\S+)\t(\S+)/;

    $tokens{$token1} = 1;
    $tokens{$token2} = 1;
}

print join("\n", keys %tokens);
print "\n";
