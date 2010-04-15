use strict;

<STDIN>;

my %ali;
my @keys;
while (<STDIN>) 
{
    chop();
    next unless ($_);
    
    my ($id, $s) = /^(\S+)\s+(.*)/;
    if (!defined $ali{$id}) 
    {
	push @keys, $id;
	$ali{$id} = "";
    }
    $ali{$id} .= $s;
    
}

foreach (@keys) 
{
    $ali{$_} ~= tr/\s//;
    printf(">%s\n%s\n", $key, $ali{$_});
}


