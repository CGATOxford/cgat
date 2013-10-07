
while(<STDIN>) {
    last if (/\.\.\.\.\|\.\./);
}
<STDIN>;
while(<STDIN>) {
    chop();
    s/~/./g;
    my ($h, $ali) = /(\S+\s+)\s(.+)/;
    $ali =~ s/\s//g;
    print $h . $ali . "\n";
}
