
my ($token, $flag) = @ARGV;

while (<STDIN>) {
    last if (/$token/);
}

print $_ unless ($flag);

while (<STDIN>) {
    print $_;
}
