
while (<STDIN>) {
    /^(\S+)\s+(\S+)/;
    print ">$1\n$2\n";
}
