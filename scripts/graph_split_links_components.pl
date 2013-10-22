$| = 1;

use strict;
use Getopt::Std;

my %opts = (); getopts('s:p:o:na:', \%opts);

my $param_suffix = $opts{'s'} || "";
my $param_prefix = $opts{'p'} || "";
my $param_offset = $opts{'o'} || 0;
my $param_nosplit = $opts{'n'} || 0;
my $param_pattern = $opts{'a'} || "";
my ($components_filename) = (@ARGV);

## read tokens
open(IN, "<$components_filename") or die "error: could not open $components_filename\n";

print "# OPTIONS:\n";
print "# prefix:	${param_prefix}\n";
print "# suffix:	${param_suffix}\n";
print "# offset:	${param_offset}\n";
print "# nosplit:	${param_nosplit}\n";
print "# pattern:       ${param_pattern}\n";

my %tokens;
while (<IN>) {
    next if (/^\[>#]/);
    last if (/^\/\//);
    my ($token, $component) = (/^(\S+)\s+(\S+)/);
    if (!defined $tokens{$token}) {
	$tokens{$token} = {};
    }
    $component += $param_offset if ($param_offset);
    ${$tokens{$token}}{$component} = 1;    
}
close(IN);

if (!$param_pattern) 
{
    $param_pattern = "${param_prefix}%s${param_suffix}";
}

my $iteration = 0;
my $last_file = "";
while (<STDIN>) 
{
    $iteration++;

    if ($iteration % 10000 == 0) {
	print "# iteration $iteration\n";
    }

    my ($token1, $token2) = /^(\S+)\s+(\S+)/;
    
    if (defined $tokens{$token1} && defined $tokens{$token2}) {
	## calculate intersection of both arrays
	my $c1;
	for $c1 (keys %{$tokens{$token1}}) {
	    if (${$tokens{$token2}}{$c1}) {
		my $file = $param_pattern;
		$file =~ s/%s/$c1/g;
		if ($last_file ne $file)
		{
		    if (!$param_nosplit)
		    {
			if ($last_file ne "") 
			{
			    close OUT;
			}
			open (OUT, ">>$file") or die "could not open $file\n";
		    
			$last_file = $file;
		    }
		}
		if ($param_nosplit) 
		{
		    print $_;
		} 
		else
		{
		    print OUT $_;
		}
	    }
        }
    }
}

close OUT unless ($param_nosplit);
