
use strict;

my @lines = <STDIN>;

use Getopt::Std;

my %opts = ();getopts('s:g:', \%opts);

my $param_sequences = $opts{'s'} || "";
my $param_gap_char = $opts{'g'} || "-";

open(IN, "<$param_sequences") or die "Could not open file $param_sequences\n";
my %sequences;
while (<IN>) {
    chop();
    my ($nid, $sequence) = split(/\t/);
    $sequences{$nid} = $sequence;
}

my $prefix_inserted = 0;
my $suffix_inserted = 0;

my @suffix;
my @prefix;

foreach (@lines) {
    chop();
    last unless ($_);
    my (@data) = split(/\t/);
    my $ali = $data[1];
    my $nid = $data[3];
    
    if (!$sequences{$nid}) {
	warn "sequence for $nid not found!\n";
	next;
    }
    
    my $seq = $ali;
    $seq =~ s/\W//g;
    $seq =~ tr/[a-z]/[A-Z]/;

    my ($extra) = ($sequences{$nid} =~ /^(\S+)$seq/);
    
    push @prefix, $extra;
    $suffix_inserted += length($extra);
    
}

foreach (@lines) {
    last unless ($_);
    my (@data) = split(/\t/);
    my $ali = $data[1];
    my $nid = $data[3];
    
    my $seq = $ali;
    $seq =~ s/\W//g;
    $seq =~ tr/[a-z]/[A-Z]/;
    my ($extra) = ($sequences{$nid} =~ /^(\S+)$seq/);
    
    push @prefix, $extra;
    $prefix_inserted += length($extra);
    
}

my $index = 0;
my $suffix = 0;
my $prefix = 0;
foreach (@lines) {
    last unless ($_);
    my (@data) = split(/\t/);
    my $nid = $data[3];
    next if (!$sequences{$nid});

    $data[1] = $param_gap_char x $prefix . 
	    $prefix[$index] . 
	    $param_gap_char x ($prefix_inserted - $prefix - length($prefix[$index])) . 
	    $data[1];
    $data[0] = 1;
    print join("\t", @data) . "\n";

    $prefix += length($prefix[$index++]);
}


    
    
