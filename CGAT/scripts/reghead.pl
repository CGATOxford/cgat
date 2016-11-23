
my ($token, $flag) = @ARGV;

while (<STDIN>) {
    last if (/$token/);
    print $_;
}

print $_ unless ($flag);


