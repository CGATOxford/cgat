################################################################################
#   Gene prediction pipeline 
#
#   $Id: split_fasta.pl 1313 2007-06-20 13:12:16Z andreas $
#
#   Copyright (C) 2004 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
#
# Splits a fasta file on STDIN into files of $STEP_SIZE
# Note: with ping you can test a range of nodes.
#	Pinging a node is not enough, check also if you can 
#	log in by scp'ing something to a target directory.
#	note: node ranges are given in real numbering (starting from 0)
#	thus, for mosix numbering, an offset of 1 has to be added.

use strict;

use Getopt::Std;

my %opts = ();getopts('n:o:p:s:i:l:t:fm:a:db:', \%opts);

$|= 1;

my $param_use_fasta_header = $opts{'f'} || 0;
my $param_num_components = $opts{'n'} || 0;
my $param_offset = $opts{'o'} || 0;
my $param_suffix = $opts{'s'} || "";
my $param_prefix = $opts{'p'} || "";
my $param_check_ping   = $opts{'i'} || "";
my $param_check_alive    = $opts{'l'} || "";
my $param_node_template = $opts{'t'} || "node%i";
my $param_use_map = $opts{'m'} || "";
my $param_pattern = $opts{'a'} || "";
my $param_pattern_input = $opts{'b'} || "^>(\S+)";
my $param_dry_run = $opts{'d'} || 0;

my $iteration = 0;
# give scp two seconds to react
my $sleep_intervall = 4;

my ($param_step_size) = @ARGV;

print "# param_use_fasta_header=$param_use_fasta_header\n";
print "# param_offset=$param_offset\n";
print "# param_num_components=$param_num_components\n";
print "# param_step_size=$param_step_size\n";
print "# param_suffix=$param_suffix\n";
print "# param_prefix=$param_prefix\n";
print "# param_check_ping=$param_check_ping\n";
print "# param_check_alive=$param_check_alive\n";
print "# param_node_template=$param_node_template\n";
print "# param_use_map=$param_use_map\n";
print "# param_pattern=$param_pattern\n";
print "# param_pattern_input=$param_pattern_input\n";
print "# param_dry_run=$param_dry_run\n";

my $nfiles = 0;

if (!$param_pattern) 
{
    $param_pattern = "${param_prefix}%s${param_suffix}";
}

if ($param_use_fasta_header) 
{
    while (<STDIN>)
    {
	if (/^>/) 
	{
	    close(OUT);
	    my ($id) = /$param_pattern_input/;
	    my $filename = $param_pattern;
	    $filename =~ s/%s/$id/g;
	    open(OUT, ">>$filename") or die "Could not open $filename\n";
	    $nfiles += 1;
	}
	print OUT $_;
    }
    close(OUT);
}
elsif ($param_use_map) 
{
    open(INFILE, "<$param_use_map") or die "Could not open $param_use_map\n";
    my %map_id2file;
    my %map_file2id;
    my $nids;
    print "reading map...";
    while (<INFILE>)
    {
	next if (/^[\#>]/);
	chomp();
	my ($id,$name) = split( /\t/ );
	$map_id2file{$id} = $name;
	$map_file2id{$name} .= "$id;";
	$nids += 1;
    }
    close(INFILE);
    
    print "done\nread $nids identifiers\n";
    
    my %seqs;
    
    my $id;
    my $keep = 0;
    while (<STDIN>)
    {
	if (/^>/) 
	{
	    ($id) = /^>(\S+)/;
	    if (defined $map_id2file{$id}) 
	    { 
		$keep = 1;
	    } 
	    else
	    { 
		$keep = 0;
	    }
	    $seqs{$id} = "" if ($keep);
	}
	$seqs{$id}.=$_ if ($keep);
    }
    
    my $name;
    for $name (keys %map_file2id)
    {
	my $ids = $map_file2id{$name};
	chop($ids);
	my $filename = $param_pattern;
	$filename =~ s/%s/$name/g;
	if ($param_dry_run)
	{
	    print "opening file file=$filename pattern=$param_pattern name=$name\n";
	    open (OUT, ">/dev/null");
	}
	else
	{
	    if (!open(OUT, ">$filename"))
		{
		    print "Could not open $filename\n";
		    next;
		}
	}
	for $id (split(/;/, $ids))
	{
	    if (!defined $seqs{$id}) 
	    {
		print "# Warning: no sequence for $id.\n";
		next;
	    }
	    print OUT $seqs{$id};
	}
	close(OUT);
	$nfiles += 1;
    }
}
else
{ 

    ## read all to find out how much data we have
    print "reading sequences...";
    my @lines = <STDIN>;
    my $nlines = $#lines + 1;
    my @headers = grep(/^>/, @lines);
    my $nheaders = $#headers + 1;
    
    print "done\nread $nlines lines for $nheaders sequences\n";
    my $nfiles = 0;

    # use ping to determine whether nodes are alive
    my @indices;
    my $chunk_size = 0;

    if ($param_check_ping || $param_check_alive) 
    {
	my ($from, $to) = split(/-/, $param_check_ping);
	my $node_id;
	foreach $node_id ($from..$to) 
	{
	    my $node=sprintf($param_node_template, $node_id-$param_offset);
	    my $retval_alive = system("$param_check_alive $node >& /dev/null ");
	    print "# node $node_id: alive=$retval_alive\n";
	    push (@indices, $node_id) if (!$retval_alive);
	}
	print "# taking nodes " . join (' ', @indices) . "\n";
	$param_num_components = $#indices + 1;
	$param_offset = 0;
    } 
    
    if ($param_num_components) 
    {
	$chunk_size = int( ($nheaders / $param_num_components) + 0.5);
	$chunk_size = 1 if ($chunk_size < 1);
    } 
    else 
    {
	$chunk_size = $param_step_size;
    }

    print "# splitting stdin into chunks of size $chunk_size\n";

    my $index = $param_offset;
    my $component = $index;
    
    $component = $indices[$index] if (@indices);

    my $filename;

    if ($param_pattern ne "")
    {
	$filename = sprintf( $param_pattern, $component);
    }
    else
    {
	$filename = sprintf("%s%05i%s", $param_prefix, $component, $param_suffix);
    }
    open (OUT, ">$filename");

    foreach (@lines) 
    {
	# continue, until you reach a different query_nid - sbjct_nid combination
	if (/^>/) 
	{
	    $iteration++;
	    
	    if ($iteration > $chunk_size) 
	    {    
		close (OUT);
	    
		$index++;
		if (@indices) 
		{
		    $component=$indices[$index];
		} 
		else 
		{
		    $component=$index;
		}
		if ($param_pattern ne "")
		{
		    $filename = sprintf( $param_pattern, $component);
		}
		else
		{
		    $filename = sprintf("%s%05i%s", $param_prefix, $component, $param_suffix);
		}
		open (OUT, ">$filename");
		$nfiles += 1;
		$iteration = 1;
	    }
	} 
	
	print OUT $_;
    }
    close (OUT);
}

print "wrote $nfiles files\n";


sub BuildFilename
{
    my ($pattern, $name) = @_;

    my $new_name = $pattern;
    $new_name =~ s/%s/$name/g;

    return $new_name;
}
