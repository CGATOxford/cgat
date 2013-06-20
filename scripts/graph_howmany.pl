use strict;

use Getopt::Std;

my %opts = (); getopts('v:', \%opts);

my $param_loglevel = $opts{'v'} || 0;

my $lastnid = 0;
my $count = 0;
my $nedges = 0;

my %nids;
my %queries;
my %subjcts;
while (<STDIN>) 
{
    next if (/^\#/);
    next if (/^>/);
    my ($nid1, $nid2) = /^(\S+)\S*\s+(\S+)\S*/;
    $nids{$nid1} = 1;
    $nids{$nid2} = 1;
    $queries{$nid1} = 1;
    $subjcts{$nid2} = 1;
    $nedges += 1;
}

my $nnids = scalar keys %nids;
my $nsubjcts = scalar keys %subjcts;
my $nqueries = scalar keys %queries;

printf("queries\tsbjcts\tvertices\tedges\n");
printf("%i\t%i\t%i\t%i\n", $nqueries, $nsubjcts, $nnids, $nedges);


if (shift(@ARGV)) {
    for (sort { $a <=> $b } keys %nids) {
	print "$_\n";
    }
}




