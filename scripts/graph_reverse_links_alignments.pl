$| = 1;

use strict;
use Getopt::Std;

my %opts = (); getopts('a', \%opts);

my $param_align = $opts{'a'} || 0;

while (<STDIN>) {
    chop();
    my ($token1, $token2, $evalue, 
	$from1, $to1, $ali1,
	$from2, $to2, $ali2, @rest) = split(/\t/);
    
    if ($param_align)
    {
	if ($token1 gt $token2)
	{
	    print "$token2\t$token1\t$evalue\t$from2\t$to2\t$ali2\t$from1\t$to1\t$ali1";
	    if (@rest) { print "\t" . join("\t", @rest); }
	} 
	else 
	{
	    print $_;
	}
    }
    else
    {
	print "$token2\t$token1\t$evalue\t$from2\t$to2\t$ali2\t$from1\t$to1\t$ali1";
	if (@rest) { print "\t" . join("\t", @rest); }
    }

    print "\n";
}





