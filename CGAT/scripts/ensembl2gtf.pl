
# This is a script which dumps out a gtf file from an
# Ensemb-schema DB with gene models.

# The option "logic_name" is not compulsory.  When not 
# used, the script assumes genes from all analyses need
# to be dumped out.

# The option "genetypes", however, must be provided. The
# script uses the list of gene biotypes to filter for a
# non-redundant list of chromosomes/seq_region names. 
# The script then fetches genes for each of these seq
# regions.

# Original author: Simon White (Ensembl)
# Modified by: Amy Tan (Sanger center)
# Modified by Andreas Heger: removed unnecessary code

use strict;

use Bio::EnsEMBL::Utils::Exception qw ( throw warning );
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::SeqFeature;
use Getopt::Long;
use Carp;
use FileHandle;

my $host   = '';
my $user   = '';
my $pass   = '';
my $port   = 3306;
my $dbname = '';
my $dnahost;
my $dnauser;
my $dnadbname;
my $dnaport;
my $dnapass = '';
my @chromosomes;
my @genetypes;
my $gtf_file = undef;
my $fix_phases = 0;
my $schema = 20;
my $localcoords = 0;
my $trim_to_slice = 0;
my $include_codons = 1;
my $logic_name;

my $start = undef;
my $end = undef;
my $seqname = undef;
my $coordsystem = 'toplevel';

$| = 1;

# local: coordinates start from the transcript start
# schema: not used

&GetOptions(
            'dbhost|host:s'        => \$host,
            'dbuser|user:s'        => \$user,
            'dbname:s'      => \$dbname,
            'dbpass|pass:s'        => \$pass,
            'dnahost:s'        => \$dnahost,
            'dnauser:s'        => \$dnauser,
            'dnadbname:s'      => \$dnadbname,	
	    'dnaport:s'      => \$dnaport,	
            'dnapass:s'        => \$dnapass,        
            #'path:s'        => \$path,
            'dbport|port:n'        => \$port,
            'start:n'       => \$start,
            'end:n'         => \$end,
            'seqname:s'     => \$seqname,
            'local'         => \$localcoords, 
            'trim_to_slice' => \$trim_to_slice,
            'codons!'       => \$include_codons,
            'chromosomes:s' => \@chromosomes,
	    'logic_name:s'  => \$logic_name,
            'genetypes:s'   => \@genetypes,
            'gtffile:s'     => \$gtf_file,
            'coordsystem:s' => \$coordsystem,
            'schema:n'      => \$schema,
            'fixphases'     => \$fix_phases,
           );



if (scalar(@genetypes)) {
  @genetypes = split (/,/, join (',', @genetypes));
}
 
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $host,
  -user   => $user,
  -port   => $port,
  -pass   => $pass,
  -dbname => $dbname
);

#$db->assembly_type($path);
if($dnadbname){
  my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                                 -host   => $dnahost,
                                                 -user   => $dnauser,
                                                 -port   => $dnaport,
                                                 -pass   => $dnapass,
                                                 -dbname => $dnadbname
                                                );
  $db->dnadb($dnadb);
}
my $sa  = $db->get_SliceAdaptor();
my $aga = $db->get_GeneAdaptor();
my $adx = $db->get_DBEntryAdaptor();
my $ga = $db->get_GeneAdaptor;

my $gtffp = new FileHandle;
$gtffp->open(">$gtf_file") or croak "Unable to open $gtf_file for write";
$gtffp->autoflush(1);


my @gene_ids;

my $sql = "select distinct(name) ".
  "from gene, seq_region ".
  "where seq_region.seq_region_id = gene.seq_region_id ".
  "and biotype = ?";
my $sth = $db->dbc->prepare($sql);

my %chr_filter_hash;

foreach my $type(@genetypes){
  $sth->execute($type);
  while(my ($name) = $sth->fetchrow){
    #push(@chromosomes,$name);    # this is deprecated as it can push redundant seq_region names into @chromosomes
    $chr_filter_hash{$name} = 1;
  }
}
@chromosomes= keys%chr_filter_hash;

my @slices;
if(@chromosomes)
{
    foreach my $name(@chromosomes)
    {
	my $slice = $sa->fetch_by_region($coordsystem, $name);
	throw("Failed to find ".$name) if(!$slice);
	push(@slices, $slice);
    }
}
else 
{
    @slices = @{$sa->fetch_all($coordsystem)};
}

foreach my $slice(@slices)
{
    #throw("Failed to find ".$name) if(!$slice);
    foreach my $type(@genetypes)
    {
	my @genes = @{$slice->get_all_Genes_by_type($type)};
	foreach my $gene(@genes)
	{
	    if ($logic_name)
	    {
		next unless $gene->analysis->logic_name eq $logic_name;
	    }
	    foreach my $trans (@{$gene->get_all_Transcripts}) 
	    {
		write_transcript_gtf($gtffp,$slice,$gene,$trans,$localcoords,
				     $trim_to_slice,$include_codons,$seqname);
	    } 
	}
    }
}

sub make_start_codon_features {
  my ($trans,$id) = @_;


  if (!$trans->translation) {
    return (());
  }

  my @translateable = @{$trans->get_all_translateable_Exons};

  my @pepgencoords = $trans->pep2genomic(1,1);

  if(scalar(@pepgencoords) > 2) {
    die("pep start does not map cleanly\n");
    next;
  } elsif (scalar(@pepgencoords) == 2) {
    print "WOW got a 2 feature start codon for " . $trans->stable_id . " strand " . $translateable[0]->strand . ", most likely a split codon.\n";
  }

  # 17 Mar 2010 at6 modified
  #
  # The script originally only handles cases where first aa residue of translation lies on an exon with start phase = 0
  # (i.e. 1st codon not split between two exons.), or cases where the first aa residue is split between two exons.
  # However, it could not handle cases where the first codon was truncated, e.g. "NNA" or "NNG" (the resulting aa residue is
  # designated "X").  Genomic sequences such as "NNA" and "NNG" would return two Mapper objects; here "NN" will be considered
  # as a Gap, and the "A" or "G" that follows will have a proper genomic coordinate.  In most of the cases, the "NN" bases 
  # did not exist in the actual genomic DNA sequence, but were just inserted to represent the truncated codon, so there is actually
  # no genomic "Gap" in the vast majority these cases.
  #
  # I therefore modified the script so it skips the "Gap" returned from these truncated codons and only look at the 
  # "Mapper::Coordinate" object that's returned.


  if(scalar(@pepgencoords) == 1) { #  start codon not split between two exons, or no truncated first codon, simplest case

    unless($pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
      die("pep start maps to gap\n");
      next;
    }
    unless($pepgencoords[$#pepgencoords]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
      die("pep start (end of) maps to gap\n");
      next;
    }
  }


  if(scalar(@pepgencoords) == 2) {  

    # looking at cases where a start codon is split between two exons, or is truncated.
    # Start codon split is allowed and we do nothing about them.
    # Nonsense cases are caught.  Also handling truncated first codon (skip the Mapper::Gap object).

    if ( ($pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) &&  ($pepgencoords[1]->isa('Bio::EnsEMBL::Mapper:Gap')) ) {
      die("Nonsense! The first codon is truncated, which is OK, but nothing from the first aa residue can be mapped to genomic coordinates!\n");
    } elsif ( $pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Gap') ) {
       print "Transcript " . $trans->stable_id . " has a truncated start codon: ". $trans->slice->seq_region_name . " " . $trans->start . " - " .  $trans->end ."\n";
       shift @pepgencoords;  # the gap is removed (ignored) so the @pepgencoords has only one Mapper::Coordinate object left.
    } 
  }

  @translateable = @{$trans->get_all_translateable_Exons};
  my @startc_feat;
  my $phase = 0;
  foreach my $pepgencoord (@pepgencoords) {
    push @startc_feat, new Bio::EnsEMBL::SeqFeature(
                             -seqname => $id,
                             -source_tag => 'starttrans',
                             -primary_tag => 'similarity',
                             -start => $pepgencoord->start,
                             -end   => $pepgencoord->end,
                             -phase => $phase,
                             -strand => $translateable[0]->strand);
    $phase = 3 - ($pepgencoord->end - $pepgencoord->start + 1);
    # print "New phase = $phase given " . $pepgencoord->start . " " . $pepgencoord->end . " " . ( $pepgencoord->end - $pepgencoord->start+1) . "\n";
  }
  if ($translateable[0]->strand == 1) {
    @startc_feat = sort {$a->start <=> $b->start } @startc_feat;
  } else {
    @startc_feat = sort {$b->start <=> $a->start } @startc_feat;
  }
  return @startc_feat;

}

sub make_stop_codon_features {
  my ($trans,$id) = @_;

  if (!$trans->translation) {
    return (());
  }
  my @translateable = @{$trans->get_all_translateable_Exons};

  my $cdna_endpos = $trans->cdna_coding_end;

  my @pepgencoords = $trans->cdna2genomic($cdna_endpos-2,$cdna_endpos);

  if(scalar(@pepgencoords) > 2) {
    die("pep end does not map cleanly\n");
    next;
  } elsif (scalar(@pepgencoords) == 2) {
    print "WOW got a 2 feature stop codon for " . $trans->stable_id . " strand " . $translateable[0]->strand . "\n";
  }

  unless($pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
    die("pep end maps to gap\n");
    next;
  }
  unless($pepgencoords[$#pepgencoords]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
    die("pep end (end of) maps to gap\n");
    next;
  }

  my @stopc_feat;
  my $phase = 0;
  foreach my $pepgencoord (@pepgencoords) {
    push @stopc_feat, new Bio::EnsEMBL::SeqFeature(
                             -seqname => $id,
                             -source_tag => 'endtrans',
                             -primary_tag => 'similarity',
                             -start => $pepgencoord->start,
                             -end   => $pepgencoord->end,
                             -phase => $phase,
                             -strand => $translateable[0]->strand);
    $phase = 3 - ($pepgencoord->end-$pepgencoord->start+1);
    # print "New phase = $phase given " . $pepgencoord->start . " " . $pepgencoord->end . " " . ($pepgencoord->end-$pepgencoord->start+1) . "\n";
  }

  if ($translateable[0]->strand == 1) {
    @stopc_feat = sort {$a->start <=> $b->start } @stopc_feat;
  } else {
    @stopc_feat = sort {$b->start <=> $a->start } @stopc_feat;
  }
  #print "Ended with " . scalar(@stopc_feat) . " stop codon features\n";
  return @stopc_feat;
}


sub write_transcript_gtf {
  my ($fh,$slice,$gene,$transcript,$localcoords,$trim_to_slice,$include_codons,$seqname) = @_;

  my $sliceoffset = 0;
  if (!$localcoords) {
    $sliceoffset = $slice->start-1;
  }

  my @startcs =  make_start_codon_features($transcript,$transcript->stable_id);
  my @endcs   =  make_stop_codon_features($transcript,$transcript->stable_id);


  my $chrname;
  $chrname = $slice->seq_region_name;

  my $idstr;

  if (defined($seqname)) {
    $idstr = $seqname;
  } else {
    $idstr = $chrname;
  }

  my ($hasstart,$hasend) = check_start_and_stop($slice,$transcript);

  if (!$include_codons) {
    $hasstart = $hasend = 0;
  }

  my @translateable_exons = @{$transcript->get_all_translateable_Exons} if $transcript->translation;


  if ($fix_phases) {
    my $phase = 0;
    foreach my $exon (@translateable_exons) {
      $exon->phase($phase);
      $exon->end_phase(($exon->length + $exon->phase) % 3);
      $phase = $exon->end_phase;
    }
  }

  my $show_trim_attrib = $trim_to_slice && ($transcript->start < 1 || $transcript->end > $slice->length);

  my $count=1;
  my $intrans = 0;
  my $instop = 0;
  
  foreach my $exon (@{$transcript->get_all_Exons}) {
    my $strand = $exon->strand;

    my $trimmed = $trim_to_slice && ($exon->start < 1 || $exon->end > $slice->length);

    if ($exon->strand == -1) {
        $strand = "-";
    } elsif ($exon->strand == 1) {
        $strand = "+";
    } elsif ($exon->strand == 0) {
        $strand = ".";
    }

    if ($transcript->translation && $exon == $transcript->translation->start_Exon) {
      $intrans = 1;
    }

    if (!$trimmed) {
      print $fh $idstr . "\t" . 
                $gene->biotype . "\t" . 
               'exon' . "\t" . 
               ($exon->start+$sliceoffset) . "\t". 
               ($exon->end+$sliceoffset) . "\t". 
               "." . "\t". 
               $strand . "\t". 
               "." . "\t";
      print_attribs($fh,$gene,$transcript,$count,'exon',$show_trim_attrib);

    # remarks don't seem to exist in the ensembl schema as they do in otter
    # gene_description is used for something else, so is not appropriate
    #if ($count == 1) { print_description_attribs($fh,$gene,$transcript,$count); }

      print $fh "\n";
    }

    if ($intrans) {

      my $cdsexon = shift @translateable_exons;
      my $phase = $cdsexon->phase;
      if ($cdsexon->phase == 1) {
        $phase = 2;
      } elsif ($cdsexon->phase == 2) {
        $phase = 1;
      } elsif ($cdsexon->phase == -1) {
        $phase = 0;
      }

      my $exon_start = $cdsexon->start;
      my $exon_end   = $cdsexon->end;
      if ($transcript->translation && 
          # $exon == $transcript->translation->end_Exon &&
          $hasend && 
          overlaps($exon, $endcs[0])) {

        if ($cdsexon->strand == 1) {
          $exon_end = $cdsexon->end - $endcs[0]->length;
        } else {
          $exon_start = $cdsexon->start + $endcs[0]->length;
        }
      }

      if ($exon_start <= $cdsexon->end &&
          $exon_end >= $cdsexon->start &&
          !$instop) {
        if (!$trimmed) {
          print $fh $idstr . "\t" . 
                   $gene->biotype . "\t" . 
                   'CDS' . "\t" . 
                   ($exon_start+$sliceoffset) . "\t". 
                   ($exon_end+$sliceoffset) . "\t". 
                   "." . "\t". 
                   $strand . "\t". 
                   $phase . "\t";
          print_attribs($fh,$gene,$transcript,$count,'CDS',$show_trim_attrib);
          print $fh "\n";
        }
      }
    }
    if ($transcript->translation && 
        $exon == $transcript->translation->start_Exon && $hasstart) {
      my $tmpcnt = $count;
      foreach my $startc (@startcs) {
        if (!$trimmed) {
          print $fh $idstr . "\t" . 
                   $gene->biotype . "\t" . 
                   'start_codon' . "\t" . 
                   ($startc->start+$sliceoffset) . "\t". 
                   ($startc->end+$sliceoffset) . "\t". 
                   "." . "\t". 
                   $strand . "\t". 
                   $startc->phase . "\t";
          print_attribs($fh,$gene,$transcript,$tmpcnt++,'start_codon',$show_trim_attrib);
          print $fh "\n";
        }
      }
    }
    if ($transcript->translation && 
        ($exon == $transcript->translation->end_Exon)) {
      if ($hasend) {
        my $tmpcnt = $count - $#endcs;
        if (!$trimmed) {
          foreach my $endc (@endcs) {
            print $fh $idstr . "\t" . 
                      $gene->biotype . "\t" . 
                      'stop_codon' . "\t" . 
                      ($endc->start+$sliceoffset) . "\t". 
                      ($endc->end+$sliceoffset) . "\t". 
                      "." . "\t". 
                      $strand . "\t". 
                      $endc->phase . "\t";
            print_attribs($fh,$gene,$transcript,$tmpcnt++,'stop_codon',$show_trim_attrib);
            print $fh "\n";
          }
        }
      }
      $intrans = 0;
    }

    if (scalar(@endcs) && overlaps($exon, $endcs[0])) {
      $instop = 1;
    }

    $count++;
  }
}

sub print_attribs {
  my ($fh,$gene,$transcript,$count,$type,$show_trim_attrib) = @_;


  my $gene_name;
  $gene_name = $gene->external_name;

  my $trans_name;
  $trans_name = $transcript->external_name;

  print $fh " gene_id \"" .  get_gene_id($gene) . "\";" .
            " transcript_id \"" . get_transcript_id($transcript) . "\";";
  print $fh " exon_number \"$count\";";
  print $fh " gene_name \"" . $gene_name . "\";" if ($gene_name);
  print $fh " transcript_name \"" . $trans_name . "\";" if ($trans_name);
  if ($type eq 'exon') {
  #  print $fh ' gbkey "mRNA";';
  } elsif ($type eq 'CDS') {
    print $fh ' evidence_id "' . get_translation_id($transcript) . '";';
  }
  if ($show_trim_attrib) {
    print $fh ' clipped_transcript "full_extent_of_transcript_beyond_bounds_of_region";'
  }
}


sub get_gene_id {
  my $gene = shift;

  if (defined($gene->stable_id)) {
    return $gene->stable_id;
  }
  return $gene->dbID;
}

sub get_transcript_id {
  my $transcript = shift;

  if (defined($transcript->stable_id)) {
    return $transcript->stable_id;
  }
  return $transcript->dbID;
}
sub get_translation_id {
  my $transcript = shift;
  my @tsf = @{$transcript->get_all_supporting_features};
  foreach my $tsf(@tsf){
    return $tsf->hseqname if($tsf->hseqname);
  }
}

sub check_start_and_stop {
  my ($slice,$trans) = @_;

  return (0,0) if (!defined($trans->translation));

  my $tln = $trans->translation;

  #$trans->sort;

  my $coding_start = $trans->cdna_coding_start;
  my $coding_end   = $trans->cdna_coding_end;
  my $cdna_seq     = uc($trans->spliced_seq);

  my $startseq     = substr($cdna_seq,$coding_start-1,3);
  my $endseq       = substr($cdna_seq,$coding_end-3,3);

  #print "Codons: " . $startseq . " " . $endseq . "\n";

  my $has_start = 1;
  my $has_end = 1;

  $has_start = 0  if ($startseq ne "ATG");
  $has_end = 0 if ($endseq ne "TAG" && $endseq ne "TGA" && $endseq ne "TAA");

  return ($has_start, $has_end);
}


sub overlaps{
  my ($f, $f2) = @_;
  return ($f->end >= $f2->start and $f->start <= $f2->end);
}
