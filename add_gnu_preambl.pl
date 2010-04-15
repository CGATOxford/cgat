## add a GNU PREAMBL to a source code file
##

use strict;

use Getopt::Std;

my %opts = ();getopts('m:t:s:', \%opts);

my $mode = $opts{'m'} || 'c++';
my $title = $opts{'t'} || "alignlib - a library for aligning protein sequences";
my $subtitle = $opts{'s'} || "\$Id\$";

my $preamble_prefix = "/*\n";
my $line_prefix = "  ";
my $line_suffix = "\n";
my $preamble_suffix = "*/\n";

my $year = "2004";
my $author = "Andreas Heger";

## use the following for alignlib:

if ($mode =~ /c\+\+/) {
    $preamble_prefix = "/*\n";
    $line_prefix = "  ";
    $line_suffix = "\n";
    $preamble_suffix = "*/\n";
    print "bhere\n";
} elsif ($mode =~ /python/) {
    $preamble_prefix = "";
    $line_prefix = "# ";
    $line_suffix = "\n";
    $preamble_suffix = "#\n";
}

##############################################
my $PREAMBL="Copyright (C) $year $author

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
";

##############################################
print $preamble_prefix;
print $line_prefix.$title.$line_suffix;
print $line_prefix.$line_suffix;
print $line_prefix.$subtitle.$line_suffix;
print $line_prefix.$line_suffix;

my $line;
for $line (split(/\n/, $PREAMBL)) {
    print $line_prefix.$line.$line_suffix;
}
print $preamble_suffix;

# print the rest
while (<STDIN>) {print $_;}



# how to strip away old headers:
# for file in *.{h,cpp}; do perl -i -n -e '{ if (!/^\/\//) { $start = 1}; print if ( $start); };' $file; done
# for file in *.{h,cpp}; do perl ~/group/scripts/tools/add_gnu_preambl.pl < $file > tmp; mv tmp $file; done







