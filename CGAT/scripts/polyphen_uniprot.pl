use strict;

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=head1 NAME

prepare_uniprot.pl

=head1 SYNOPSIS

prepare_uniprot.pl [options]

Where options are:

 -s <dat>[,<idm>] specify local source file(s): <dat> with the protein entries
                  in UniProtKB format and optional <idm> with idmapping data;
                  default is to fetch all source files from ftp.uniprot.org
 -k               keep source files; default is to remove source files after
                  they were successfully processed
 -h               print this help
 -t <name>        output files for this species. The default is <name>=human.
 -u <name>        uniprot database to downlead. The default is <name>human.
                  Alternatives are 'mammals', 'rodents', .... See
                  ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/
 -o <org>         The two part species name to extract. The default is <org>"homo sapiens".
                  Uniprot adds abbreviations to the OS lines, and these work as well (mouse = Mus musculus)

=head1 DESCRIPTION

Downloads uniprot_trembl_human.dat.gz, uniprot_sprot_human.dat.gz, and
idmapping_selected.tab.gz from ftp.uniprot.org and prepare protein sequence
database and database ID index files for PolyPhen-2

Modified by Andreas Heger - load feature tables for mouse as well.

=head1 AUTHOR

 Vasily Ramensky
 Steffen Schmidt
 Ivan Adzhubey

=head1 SUBVERSION

 $LastChangedDate: 2010-03-29 01:33:25 -0400 (Mon, 29 Mar 2010) $
 $LastChangedRevision: 288 $
 $LastChangedBy: ivan $

=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

use LWP::Simple;
use Getopt::Std;
use Pod::Usage;
use Storable qw(nstore);

my @fttokens = qw(
  DISULFID THIOLEST THIOETH MOD_RES LIPID CARBOHYD METAL BINDING ACT_SITE SITE
  SIGNAL TRANSIT PROPEP TRANSMEM INIT_MET
  ## PEPTIDE CA_BIND DNA_BIND NP_BIND ZN_FING
);

my $UNIPROT_FTP  = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase';
#-----------------------------------------------------------------------------

my %opts;
# changed by AH
getopts('hks:t:u:o:', \%opts);

exists $opts{'h'} and pod2usage({ -verbose=>99 });

# changed by AH
my $species = $opts{'t'} || "human"; 
my $uniprot = $opts{'u'} || "human"; 
my $organism = $opts{'o'} || "homo sapiens"; 
my $NAME = 'uniprot_' . $species;

my $seqfile     = "$NAME.seq";
my $ftfile      = "$NAME.ft";
my $swfile      = "$NAME.dat";
my $id2accfile  = "$NAME.id2acc";
my $snp2accfile = "$NAME.snp2acc";
my $xdb2accfile = "$NAME.xdb2acc";

#-----------------------------------------------------------------------------

my @names = ( "uniprot_sprot_${uniprot}.dat",
	      "uniprot_trembl_${uniprot}.dat", 
	      "idmapping_selected.tab" );

if (defined $opts{'s'}) {

  my @files = split /,/, $opts{'s'};
  foreach my $file (@files) {
    foreach my $ext (qw{ .seq .ft .dat .id2acc .snp2acc .xdb2acc}) {
      die "Can't simultaneously read from and write to $file, please rename input file\n"
        if $file eq "$NAME$ext";
    }
  }
  @names = @files;

} else {

  my $suf = '.gz';

  foreach my $name (@names) {
    my $infile;
    if ($name =~ /^idmapping/) {
      $infile = "$UNIPROT_FTP/idmapping/${name}${suf}";
    } else {
      $infile = "$UNIPROT_FTP/taxonomic_divisions/${name}${suf}";
    }
    my $outfile = "${name}${suf}";
    unlink $outfile; ###  otherwise will not open them for writing
    warn "Fetching $infile ...\n";
    die "Can't fetch file: $infile\n" unless mirror($infile, $outfile) == 200;
    warn "Decompressing $outfile ...\n";
    `gzip -d $outfile`;
    warn "Done.\n";
  }

}

open(SEQ,     ">$seqfile") || die "Can't open $seqfile";
open(FT,      ">$ftfile" ) || die "Can't open $ftfile";
open(SWALL,   ">$swfile" ) || die "Can't open $swfile";
open(ID2ACC,  ">$id2accfile" ) || die "Can't open $id2accfile";
open(SNP2ACC, ">$snp2accfile" ) || die "Can't open $snp2accfile";

my $dbase_len = 0;

foreach my $file (@names) {

 next if $file =~ /idmapping/;

 warn "Processing $file ...\n";

 my $num_human = my $copy = 0;
 my (@entrylines, @snps, @ftlines, @seqlines, @accs, @emblrefs);
 my ($id, $seq, $de, $len, $varpos);

 open(TMP, $file) || die "Can't open $file\n";

 LINE: while (<TMP>) {

  if(/^ID.*\b(\d+) AA/) {
   @seqlines = @snps = @ftlines = @entrylines = @accs = @emblrefs = ();
   $copy = 0;
   $id = (split)[1];
   $len = $1;
   undef $de;
   undef $varpos;
  } # end if
  elsif( /^OS/ ) {
   $copy++ if m/${organism}/i; ## there exist multiple OS-lines
  } # end if
  elsif( /^AC/ ) {
   s/;//g;
   my @accs_tmp = (split(/\s+/, $_));
   shift @accs_tmp; # remove ^AC
   push @accs, @accs_tmp;
  } # end elsif
  elsif( /^FT   (\w+)\s/ ) { ## important to put exact spaces in regexp
    my $feature = $1;
    if ($feature eq 'VARIANT') {
      # Read and merge all continuation lines for FT VARIANT
      # since the 'dbSNP:' tag(s) can appear on any of them.
      chomp(my $varlines = $_);
      my $inline;
      while ($inline = <TMP>) {
        last unless substr($inline, 2, 32) eq ' ' x 32;
        chomp($varlines .= ' ' . substr($inline, 34));
      }
      my $from = substr($varlines, 14, 6); $from =~ s/\s+//g;
      my $to   = substr($varlines, 21, 6); $to   =~ s/\s+//g;
      # Validate that the line annotates a single residue substitution
      if ($from =~ /^\d+$/ && $to =~ /^\d+$/ && $to - $from == 0) {
        my $description = substr($varlines, 34);
        # There can be multiple SNPs annotated at the same position
        while ($description =~ /\b([A-Z]) -> ([A-Z])\s+\(in (.+?)\)/g) {
          my ($aa1, $aa2, $pointer) = ($1, $2, $3);
          # We have a dbSNP annotation
          if ($pointer =~ /dbSNP:(rs\d+)/) {
            push @snps, [ $1, $from, $aa1, $aa2 ];
          }
        }
      }
      # since we already have next input line in the $inline buffer
      $_ = $inline, redo LINE;
    } else {
      push @ftlines, $_ if grep { $1 eq $_ } @fttokens;
    }
  } # end elsif
  elsif( /^DE\s+(\w.*)\n$/ ) {
   $de .= " $1";
  } # end elsif
  elsif( /^ / ) {
   $seq = $_;
   $seq =~ s/ //g;
   push @seqlines, $seq;
   $dbase_len += length($seq);
  } # end elsif
  elsif( /^\/\// ) {
   # changed by AH
   if($copy && $id !~ /ALU\d_${species}/i) {
    print SWALL @entrylines;
    print SWALL "//\n";

    my $primacc = $accs[0];

    chomp(@seqlines);
    print SEQ ">sp|$primacc|$id$de LENGTH: $len AA\n";
    print SEQ join '', @seqlines, "\n";
    print FT map { "$primacc $_" } @ftlines;
    print ID2ACC map { "$_ $primacc\n" } @accs;
    print ID2ACC "$id $primacc\n";
    print SNP2ACC map { "$_->[0]\t$primacc\t" . join("\t",@$_[1..3]) . "\n" } @snps if @snps;

    $num_human++;

   } # end if
  } # end elsif

  push @entrylines, $_;

 } # end while

 close TMP;

 # changed by AH
 warn "Found $num_human ${species} proteins.\n";

 unlink $file unless defined $opts{'k'};

} # end foreach

warn "Total length of database $NAME: $dbase_len\n";

close SNP2ACC; close ID2ACC; close SWALL; close FT; close SEQ;

warn "Gzipping $NAME.dat ...\n";
`gzip $NAME.dat`;
warn "Done.\n";

if (@names > 1) {
  my $file = $names[-1];
  warn "Processing $file ...\n";
  my @dbtags    = qw( acc pid geneid refseq gi );
  my %db2acc;
  open DBFIN, '<', $file or die "Can't open $file: $!\n";
  while (<DBFIN>) {
    # We only process first 5 columns: acc, pid, geneid, refseq, gi
    # Note: acc and pid are both unique (acc<->pid) and only primary
    # accessions are included. All other IDs can be lists of multiple
    # emtries separated by '; ' characters.
    my @a = split /\t/;
    # changed by AH
    next unless $a[1] =~ /^\w+_${species}$/i;  # skip all but human proteins
    my $acc = $a[0];
    for (my $i=3; $i<=4; $i++) {
      my $dbids = $a[$i];
      next unless length $dbids;
      foreach my $dbid (split /;\s*/, $dbids) {
        $dbid =~ s/^\s+//; $dbid =~ s/\s+$//;
        next unless length $dbid;
        # Skip all except the first xdb ID encountered in case of ambiguous mappings
        next if exists $db2acc{$dbtags[$i]}{$dbid};
        $db2acc{$dbtags[$i]}{$dbid} = $acc;
      }
    }
  }
  close(DBFIN);
  nstore \%db2acc, $xdb2accfile or die "Can't store $xdb2accfile: $!\n";
  warn "File saved: $xdb2accfile\n";
  warn "Warning: file has zero size: $xdb2accfile\n" unless -s $xdb2accfile;
  unlink $file unless defined $opts{'k'};
}

warn "All done.\n";

#-----------------------------------------------------------------------------
