use strict;


my @patterns;

my ($filename_patterns) = = @ARGV;

open(IN, "<$filename_patterns") or die "$filename_patterns not found.";

while (<IN>) 
{
    chop();
    ($old, $new) = split(\t);
    push $patterns ($old, $new);
}
close(IN);

while (<STDIN>) 
{
    for pattern in $patterns


}



