use strict;
use Getopt::Std;

my %opts = (); getopts('drm', \%opts);

my $param_no_merge = $opts{'m'} || 0;
my $param_direction = $opts{'d'} || 0;
my $param_reverse      = $opts{'r'} || 0;

my ($param_filter) = shift(@ARGV);

open(IN, "<$param_filter") or die "could not open $param_filter\n";

my %links;

while (<IN>) {
    chop();
    my ($token1, $token2, $rest) = /^([^\t]+)\t([^\t]+)\t*(.*)$/;
    if (!$param_no_merge) {
	$rest = "\t" . $rest;
    } else {
	$rest = "";
    }
    $links{"$token1-$token2"} = $rest;
    $links{"$token2-$token1"} = $rest unless ($param_direction);
}    

close(IN);

while (<STDIN>) {
    chop();
    if ( /^\#/ ) { print $_. "\n"; next;}
    if ( /^>/ ) { print $_. "\n"; next;}

    my ($token1, $token2) = /^([^\t]+)\t([^\t+]+)/;

    my $e1 = exists $links{"$token1-$token2"};
    my $e2 = exists $links{"$token2-$token1"};
    
    if ($param_reverse) 
    {
	if ($param_direction)
	{
	    print $_ . "\n" if (!$e1);
	}
	else 
	{
	    print $_ . "\n" if (!$e1) && (!$e2);
	}
    }
    else
    {
	if ($e1)
	{
	    print $_ . $links{"$token1-$token2"} . "\n";
	} 
	else
	{
	    if (!$param_direction) 
	    {
		print $_ . $links{"$token2-$token1"} . "\n" if ($e2);
	    }
	}
    }
}
 
# print join("-\n", keys %links);
# print "\n";


