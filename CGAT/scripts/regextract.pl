
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

