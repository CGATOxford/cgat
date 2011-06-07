
# script to extract repeats from an ENSEMBL database
# in GFF format.
#
# This scripts uses the BioPerl/Ensembl API.

use strict;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw ( throw warning );
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::SeqFeature;
use Getopt::Long;
use Carp;
use FileHandle;

my $host   = '';
my $user   = '';
my $pass   = '';
my $port   = 5306;
my $dbname = '';
my @chromosomes;
my @repeattypes;
my $gtf_file = undef;
my $fix_phases = 0;
my $localcoords = 0;
my $trim_to_slice = 0;
my $include_codons = 1;
my $logic_name;

my $start = undef;
my $end = undef;
my $seqname = undef;
my $coordsystem = 'toplevel';

$| = 1;

&GetOptions(
    'host|h:s'        => \$host,
    'user|u:s'        => \$user,
    'name|d:s'      => \$dbname,
    'pass|p:s'        => \$pass,
    'port:n'        => \$port,
    'start:n'       => \$start,
    'end:n'         => \$end,
    'repeattypes:s'   => \@repeattypes,
    'coordsystem:s' => \$coordsystem,
    );

## connect to ENSEMBL database
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $host,
  -user   => $user,
  -port   => $port,
  -pass   => $pass,
  -dbname => $dbname
);

# build filter for repeat types
my $repeattypes_regex = "";

if (scalar(@repeattypes)) {
  @repeattypes = split (/,/, join (',', @repeattypes));
  $repeattypes_regex = join( "|", @repeattypes);
}

my $slice_adaptor = $db->get_SliceAdaptor();

my @slices = @{ $slice_adaptor->fetch_all($coordsystem) };

my $rfa = $db->get_RepeatFeatureAdaptor();
my $x = 0;
my %found_repeat_types;

# fetch_all_by_Slice returns duplicate entries. In the
# absence of any insight on how to turn this off, duplicate
# locations are removed manually.
my %output_repeats;

foreach my $slice (@slices)
{
    my @repeats = @{$rfa->fetch_all_by_Slice($slice, undef)};
    foreach my $repeat (@repeats) 
    {
	$found_repeat_types{$repeat->repeat_consensus()->repeat_type()} = 1;
	next unless ($repeat->repeat_consensus()->repeat_type() =~ m/${repeattypes_regex}/);
	my $start = $repeat->slice()->start() + $repeat->start() -1;
	my $end = $repeat->slice()->start() + $repeat->end() -1;
	my $key = $repeat->slice()->seq_region_name() . $start . $end;
	next if ($output_repeats{$key});
	$output_repeats{$key} = 1;
	my $strand = ( $repeat->strand == 1 ? '+' : '-');
	$strand = '.' if ($repeat->strand == 0);

	printf("%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\tgene_id \"%06i\"; transcript_id \"%06i\"; type \"%s\"\n",
	       $repeat->slice()->seq_region_name(),
	       "repeat",
	       "exon",
	       $start,
	       $end,
	       $strand,
	       ".",
	       ".",
	       ++$x, 
	       $x,
	       $repeat->repeat_consensus()->repeat_type());
    }
}

printf( "# %s\n", join( ',', keys %found_repeat_types));
