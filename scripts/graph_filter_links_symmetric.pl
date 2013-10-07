
use strict;

my %pairs;

while (<STDIN>) {
    
  my ($token1, $token2) = /^(\S+)\s+(\S+)/;
  
  my $key1 = "$token1-$token2";
  my $key2 = "$token2-$token1";

  if ( !defined $pairs{$key1} && !defined $pairs{key2} ) {
      print $_;
      $pairs{$key1} = 1;
      $pairs{$key2} = 1;
  }
}




