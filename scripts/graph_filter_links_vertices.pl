## filter links with a list of tokens
## usage:
## perl filter_links_tokens.pl tokens < in > out
##
## Options:
## -s: sufficient if one vertex matches
## -r: reverse match: delete if a vertex matches
## $Id$

use strict;
use Getopt::Std;

my %opts = ();getopts('sr', \%opts);

my $param_single  = $opts{'s'} || 0;
my $param_reverse = $opts{'r'} || 0;

my ($tokens_filename) = shift(@ARGV);

## read tokens
open(IN, "<$tokens_filename") or 
    die "error: could not open $tokens_filename\n";

my %tokens;
while (<IN>) {
    my ($token) = (/^(\S+)/);
    $tokens{$token} = 1;
}
close(IN);

while (<STDIN>) {

    my ($token1, $token2) = /^(\S+)\s+(\S+)/;
    my $match = 0;

    if ($param_single) 
    {
	$match= $tokens{$token1} || $tokens{$token2}
    } 
    else 
    {
    	$match = $tokens{$token1} && $tokens{$token2};
    }
    
    if ($param_reverse)
    {
	print $_ unless ($match);
    }
    else
    {
	print $_ if ($match);
    }
}




