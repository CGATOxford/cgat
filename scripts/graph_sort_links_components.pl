$| = 1;

use strict;
use Getopt::Std;

my %opts = (); getopts('o:f:s:', \%opts);

my $param_offset = $opts{'o'} || 0;
my $param_format = $opts{'f'} || "fasta";
my $param_sort = $opts{'s'} || "id";
my ($components_filename) = (@ARGV);

## read tokens
open(IN, "<$components_filename") or die "error: could not open $components_filename\n";

print "# OPTIONS:\n";
print "# format:	${param_format}\n";
print "# offset:	${param_offset}\n";
print "# sort:		${param_sort}\n";

my %tokens;
my %sizes;
while (<IN>) {
    next if (/^\#/);
    last if (/^\/\//);
    my ($token, $component) = (/^(\S+)\s+(\S+)/);
    if (!defined $tokens{$token}) {
	$tokens{$token} = {};
    }
    ${$tokens{$token}}{$component + $param_offset} = 1;    
    $sizes{$component+$param_offset} += 1;
}
close(IN);

my $iteration = 0;
my $last_file = "";
my $tmpfile = "tmp_$$";

open(OUT, ">$tmpfile") or die "Could not open temporary file.";

while (<STDIN>) 
{
    $iteration++;

    if ($iteration % 10000 == 0) {
	print "# iteration $iteration\n";
    }

    my ($token1, $token2) = /^(\S+)\s+(\S+)/;
    
    if (defined $tokens{$token1} && defined $tokens{$token2} ) 
    {
	## calculate intersection of both arrays
	my $c1;
	for $c1 (keys %{$tokens{$token1}})
	{
	    if (${$tokens{$token2}}{$c1}) 
	    {
		print OUT "$c1\t$sizes{$c1}\t$_";
	    }
	}
    }
	
}

close OUT;

if ($param_sort eq "id") 
{ 
    open(IN, "sort -k1,1n $tmpfile |") or die "Could not open temporary file.";
}
elsif ($param_sort eq "size")
{
    open(IN, "sort -k2,2n -k1,1n $tmpfile |") or die "Could not open temporary file.";    
}

my $last_id = -1;
while (<IN>) 
{
    my @data = split(/\t/);
    if ($data[0] != $last_id) 
    {
	$last_id = $data[0];
	if ($param_format == "fasta") 
	{
	    print ">cluster# $last_id size= $sizes{$last_id}\n";
	}
    }

    print join("\t", @data[2..$#data]);
}
 
system("rm $tmpfile");



