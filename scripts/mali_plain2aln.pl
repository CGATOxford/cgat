
use strict;


my $param_line_width = 60;

my @lines = <STDIN>;

my $lali = 0;

foreach (@lines) {
    my (@data) = split(/\t/);
    my ($from, $ali, $to, $nid) = @data;    
    $lali = length($ali);
    last;
}

print "CLUSTAL";

my $current_pos = 0;
while ($current_pos < $lali) {
    print "\n\n";
    my %nids;

    my $current_nid = 0;

    foreach (@lines) {
	chop();
	last unless ($_);
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
	
	printf("%-20s %s\n", $nid, substr($ali, $current_pos, $param_line_width));
    }
    $current_pos += $param_line_width;
}

print "\n";

