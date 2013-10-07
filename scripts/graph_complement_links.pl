use strict;
use Getopt::Long;

my $usage = <<'USAGE';

Usage: perl graph_complement_links.pl [OPTIONS]< links.in > links.out

Takes a graph and adds missing links so that the resulting graph is complete.
Clusters in the graph are separated by '>'.

Options:
    --missing		only print missing links
    --directed		graph is directed
    --substitute #	substitute missing links from file #
USAGE

# optional parameters
my $param_help = undef;
my $param_only_missing = 0;
my $param_directed_graph = 0;
my $param_file_size_limit = 100000;
my $param_substitute_missing = 0;
GetOptions(
	   'help' => \$param_help,
	   'missing:s' => \$param_only_missing,
	   'directed:s' => \$param_directed_graph,
	   'substitute:s' => \$param_substitute_missing,
	   );

die $usage if ($param_help);

use strict;

my %replacement_links;
if ($param_substitute_missing)
{
    open(IN, "<$param_substitute_missing") or die "Could not open file $param_substitute_missing.\n";
    while (<IN>)
    {
	my ($token1, $token2) = /^(\S+)\s+(\S+)/;
	my $key1 = "$token1-$token2";
	$replacement_links{$key1} = $_;
    }
    close(IN);
}

my %tokens;
my %links;
while (<STDIN>) 
{
    if (/^>/) 
    {
	PrintMissingLinks( \%tokens, \%links, \%replacement_links );
	undef %tokens;
	undef %links;
	print $_;
	next;
    }

    my ($token1, $token2) = /^(\S+)\s+(\S+)/;
    my $key1 = "$token1-$token2";
    my $key2 = "$token2-$token1";
    
    $links{$key1} = 1;
    $links{$key2} = 1 unless ($param_directed_graph);
    $tokens{$token1} = 1;
    $tokens{$token2} = 1;
    print $_ unless $param_only_missing;
}

PrintMissingLinks( \%tokens, \%links, \%replacement_links );


sub PrintMissingLinks 
{
    my ($tokens, $links) = @_;

    my @keys = keys %$tokens;
    my ($x, $y);
    for $x (0..$#keys - 1) 
    {
	my $token1=$keys[$x];
	for $y ($x+1..$#keys)
	{
	    my $token2=$keys[$y];
	    my $key1 = "$token1-$token2";		
	    my $key2 = "$token2-$token1";
	    if (!$links{$key1})
	    {
		if ($replacement_links{$key1})
		{
		    print $replacement_links{$key1};
		} 
		elsif (!$param_directed_graph && $replacement_links{$key2})
		{
		    print $replacement_links{$key2};	    
		}
		else
		{
		    print "$token1\t$token2\tmissing\n";
		}
	    }
	    if ($param_directed_graph && !$links{$key2})
	    {
		if ($replacement_links{$key2})
		{
		    print $replacement_links{$key2};    
		}
		else
		{
		    print "$token2\t$token1\tmissing\n";
		}
	    }
	}
    }
}
