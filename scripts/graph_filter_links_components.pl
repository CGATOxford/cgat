## $Id$
##
## usage:
## perl filter_links_components.pl [OPTIONS] components < graph.in
##
## filter links according to components. Only links inside the same
## component are kept.
##
$| = 1;

use strict;
use Getopt::Std;

my %opts = (); getopts('o:f:s:', \%opts);

my ($components_filename) = (@ARGV);

## read tokens
open(IN, "<$components_filename") or die "error: could not open $components_filename\n";

my %tokens;

while (<IN>) {
    next if (/^\#/);
    last if (/^\/\//);
    my ($token, $component) = (/^(\S+)\s+(\S+)/);
    $tokens{$token} = $component;    
}
close(IN);

my $iteration = 0;

while (<STDIN>) 
{
    $iteration++;

    if ($iteration % 10000 == 0) {
	print "# iteration $iteration\n";
    }

    my ($token1, $token2) = /^(\S+)\s+(\S+)/;
    
    if (defined $tokens{$token1} && defined $tokens{$token2} &&
	$tokens{$token1} eq $tokens{$token2} ) 
    {
	print;
    }
	
}



