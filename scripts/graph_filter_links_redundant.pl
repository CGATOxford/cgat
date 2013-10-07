
use strict;
use Getopt::Std;

my %opts = ();getopts('d', \%opts);

my $param_direction = $opts{'d'} || 0;

my %pairs;

while (<STDIN>) {
  
  my ($token1, $token2) = /^(\S+)\s+(\S+)/;
  
  my $key1 = "$token1-$token2";
  my $key2 = "$token2-$token1";

  if ( !defined $pairs{$key1})
  {
      $pairs{$key1} = 1;
      if ($param_direction)
      {
	  if (!defined $pairs{key2})
	  {
	      $pairs{$key2} = 1;
	      print $_;
	  }
      }
      else
      {
	  print $_;
      }
  }
}




