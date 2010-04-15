use strict;

<STDIN>;

my ($num_lines, $width) = /(\d+)\s+(\d+)/;


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
    $ali{$_} =~ s/\s//g;
    printf(">%s\n%s\n", $_, $ali{$_});
}


