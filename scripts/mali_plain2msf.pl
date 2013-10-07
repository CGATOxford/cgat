
use strict;

my @lines = <STDIN>;

my %nids;

print <<EOF_HEADER
PileUp



   MSF:   82  Type: P    Check:  4661   ..

EOF_HEADER
    ;

my $current_nid = 0;
foreach (@lines) {
    chop();
    last unless ($_);
    my (@data) = split(/\t/);
    my ($from, $ali, $to, $nid) = @data;    
    
    $nid = ++$current_nid unless ($nid);

    if ($nids{$nid}) {
	$nids{$nid} += 1;
	$nid .= chr(ord('a') + $nids{$nid} ); 
    } else {
	$nids{$nid} = 1;
    }

    printf(" Name: %s oo Len: %i  Check: %i  Weight: %4.2f\n", $nid, 82, 1000, 1.0);
}

undef %nids;

print "\n//\n\n\n";
$current_nid = 0;
foreach (@lines) {
    my (@data) = split(/\t/);

    my ($from, $ali, $to, $nid) = @data;    
    
    $ali =~ s/\-/\./g;
    $nid = ++$current_nid unless ($nid);

    if ($nids{$nid}) {
    	$nids{$nid} += 1;
	$nid .= chr(ord('a') + $nids{$nid} ); 
    } else {
	$nids{$nid} = 1;
    }

    printf("%-20s %s\n", $nid, $ali);
}

print "\n";



    
    
