################################################################################
#   Gene prediction pipeline 
#
#   $Id: extract_fasta.pl 2233 2008-10-28 10:55:23Z andreas $
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
## select lines from a fasta file according to a list of nids
## 
## Also works with domain sequences.
## This file will remove any white space characters and put sequences on
## a single line unless -e is set (for preserve)
## 
## Set token (-t) to '@' to split fastq files.
##
use strict;

use Getopt::Std;

my %opts = ();getopts('mu:l:s:n:f:rp:et:', \%opts);

my $param_min_length = $opts{'l'} || 0;
my $param_max_length = $opts{'u'} || 100000000000;
my $param_file_nids = $opts{'f'} || "";
my $param_masked = $opts{'m'} || 0;
my $param_max_segment = $opts{'s'} || 5;
my $param_nid = $opts{'n'} || "";
my $param_pattern = $opts{'p'} || "";
my $param_remove = $opts{'r'} || 0;
my $param_preserve = $opts{'e'} || 0;
my $param_token = $opts{'t'} || ">";

my $pattern = "^$param_token(\\S+)";

my %nids;

my $use_all_nids = 0;

if ($#ARGV == 0)
{
    ($param_file_nids) = @ARGV;
}


if ($param_file_nids) 
{
    open(IN1, "<$param_file_nids") or die "Could not open $param_file_nids\n";

    while(<IN1>) { 
	chop(); 
	my @data=split(/\s+/);
	$nids{$data[0]} = 1;
    }
} 
elsif ($param_nid)
{
    $nids{$param_nid} = 1;
} else 
{
    $use_all_nids = 1;
}
    
close(IN1);

my $keep = 0;
my $sequence = "";

sub printSequence 
{
    my ($keep, $sequence) = @_;

    if (!$param_preserve)
    {
	$sequence =~ s/\s//g;
	$sequence .= "\n";
    }
    my $l = length($sequence);
    if ($param_masked) 
    {
	my $s = $sequence;
	$s =~ s/X[A-WYZ]{1,$param_max_segment}X/XX/g;
	$s =~ s/X[A-WYZ]{1,$param_max_segment}X/XX/g;
	$s =~ s/X//gi;
	$l = length($s);
    }
    if ($l >= $param_min_length &&
	$l <= $param_max_length) 
    {
	print $keep;
	print $sequence; 
    }
}

while (<STDIN>) 
{
    next if (/^#/);

    if (/$pattern/)
    {
	if ($keep) 
	{
	    printSequence( $keep, $sequence);
	}
	$sequence = "";
	if ( ( !$param_pattern &&  $use_all_nids || 
	     ($nids{$1} && !$param_remove) ||
	     (!$nids{$1} && $param_remove) ) || 
	     ($param_pattern && /$param_pattern/ ))
	{
	    $keep = $_;
	    next;
	} 
	else 
	{ 
	    $keep = 0;
	};
    }
    
    $sequence .= $_;
}

printSequence( $keep, $sequence) if ($keep);


