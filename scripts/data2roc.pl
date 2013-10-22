use strict;
use Getopt::Std;

my %opts = ();getopts('sc:d:rzV:', \%opts);

my $param_remove    = $opts{'r'} || 0;
my $param_select    = $opts{'s'} || 0;
my $param_cutoff    = $opts{'c'} || 0;
my $param_select_dir= $opts{'d'} || "";
my $param_check_zip = $opts{'z'} || 0;
my $param_loglevel  = $opts{'V'} || 0;

my $val = 0;
my $tp  = 0;
my $p   = 0;

while (<STDIN>) {
    next if (/^\#/);
    chop();
    $p++;
    my (@data) = split(/\t/);

    if ($data[0] eq "+")
	$tp += 1;
    
    printf("%f\t%f\n");
