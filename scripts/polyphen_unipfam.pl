use warnings;
use strict;

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=head1 NAME

unipfam.pl - a simple script to extract Pfam annotations for UniProt proteins

=head1 SYNOPSIS

unipfam.pl [options] [seq_file]

Where options are:

 -t <name>        output files for this species. The default is <name>=human.

=head1 DESCRIPTION

This script uses UniProt sequence database file (default: 'uniprot_human.seq')
prepared by 'prepare_uniprot.pl' script to get a list of protein accessions
and then parses swisspfam.gz to extract annotations for all proteins.
It creates 'unipfam.annot' file in the current directory.

Modified by Andreas Heger - use for other species than human.

=head1 SUBVERSION

 $LastChangedDate: 2009-11-19 19:03:22 -0500 (Thu, 19 Nov 2009) $
 $LastChangedRevision: 246 $
 $LastChangedBy: ivan $

=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

use LWP::Simple;
use Getopt::Std;
# changed by AH
my %opts;
getopts('t:', \%opts);
my $species = $opts{'t'} || "human"; 


my $uniseq_name = "uniprot_${species}.seq";
$uniseq_name = shift if @ARGV;
my $outfile = "unipfam_${species}.annot";

open(F, $uniseq_name) or die "Can't open input file: $uniseq_name\n";

my %HsAcc;
warn "Reading protein accessions from $uniseq_name...\n";
while (<F>) {
    /^>/ or next;
    my (undef, $acc) = split /\|/;
    $HsAcc{$acc} = 1;
}
close F;
warn "Done.\n";

# Download recent swisspfam.gz file
my $PFAM_FTP = 'ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/swisspfam.gz';
warn "Fetching swisspfam.gz...\n";
unlink 'swisspfam.gz' if -e 'swisspfam.gz';
die "Can't fetch file: swisspfam.gz\n" unless mirror($PFAM_FTP, 'swisspfam.gz') == 200;
warn "Done.\n";

warn "Extracting Pfam annotations...\n";
open(F, 'gunzip -c swisspfam.gz |') or die "Can't open input file: swisspfam.gz\n";
open(FOUT, ">$outfile") or die "Can't create output file: $outfile\n";
while (<F>) {
    /^>/ or next;
    my (undef, undef, $acc) = split;
    $acc =~ s/\.\d+$//;
    exists $HsAcc{$acc} or next;
    while (<F>) {
        /^\s*$/ and last;
        /^Pfam-B/ and next;
        chomp;
        my ($id, $desc, @f) = split;
        my @loc;
        while ($f[$#f] =~ /^\d+-\d+$/) {
            push @loc, pop @f;
        }
        while ($f[$#f] !~ /^PF/) {
            $desc = (pop @f)." $desc";
        }
        print FOUT join("\t", $acc, $f[$#f], $id, $desc, @loc), "\n";
    }
}
close FOUT;
close F;
unlink 'swisspfam.gz' if -e 'swisspfam.gz';
warn "Done.\n";
