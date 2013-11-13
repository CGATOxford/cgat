use strict;
use Getopt::Std;

my %opts = ();getopts('s:f:cw:', \%opts);

my $param_file_sequences = $opts{'s'} || 0;
my $param_file_fasta = $opts{'f'} || 0;
my $param_offset_correction = $opts{'c'} || 0;
my $param_default_weight = $opts{'w'} || 10;

my %sequences;
my %offsets;

my @lines = <STDIN>;

## collect all pairs and delete from sequences those
## without any links
foreach (@lines) {
    my ($id1, $residue1, $id2, $residue2) = (/^(\S+)-(\d+).*\s(\S+)-(\d+)/);
    $sequences{$id1} = 1;
    $sequences{$id2} = 1;
}

if ($param_file_sequences) {
    open(IN, "<$param_file_sequences") or die "Could not open plain sequences in $param_file_sequences";
    while (<IN>) {
	chop();
	my ($id, $seq) = split(/\t/);
	$sequences{$id} = $seq if ($sequences{$id});
    }
    close(IN);
}

if ($param_file_fasta) {
    open(IN, "<$param_file_fasta") or die "Could not open fasta sequences in $param_file_fasta";
    my $s = "";
    my $id = "";
    while (<IN>) {
	chop();
	if (/^>/) {
	    if ($id) {
		$s =~ s/\s//;
		$sequences{$id} = $s if ($sequences{$id});
	    }
	    ($id) = /^>(\S+)/;
	    $s = "";
	} else {
	    $s .= $_;
	}
    }
    $s =~ s/\s//;
    $sequences{$id} = $s if ($sequences{$id});
}

if ($param_offset_correction) {
    foreach (keys %sequences) {
	my ($offset) = /\d+_(\d+)/;
	$offset--;
	$offsets{$_} = $offset;
    }
}

## map sequences to numbers.
my @x = keys %sequences;
my $nsequences = $#x + 1;

my %identifiers;
printf("%i\n", $nsequences);
my $x = 1;
foreach (keys % sequences) {
    $identifiers{$_} = $x++;
    printf("%s %i %s\n", $_, length($sequences{$_}), $sequences{$_});
}



print "!comment\n";

my $last_id1 = "";
my $last_id2 = "";
foreach (@lines) {
    my ($id1, $residue1, $id2, $residue2) = (/^(\S+)-(\d+).*\s(\S+)-(\d+)/);
    if ($last_id1 ne $id1 || $last_id2 ne $id2) {
	printf("#%i %i\n", $identifiers{$id1}, $identifiers{$id2});
	$last_id1 = $id1;
	$last_id2 = $id2;
    }
    printf("%i %i %i\n", $residue1 - $offsets{$id1}, $residue2  - $offsets{$id2}, $param_default_weight);
}

print "! SEQ_1_TO_N\n";





