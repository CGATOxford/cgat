$|=1;

# use strict;
# global variables for the parser :
use vars '$query';
use vars '$query_length'; 
use vars '$sbjct_length'; 

use vars '$began';       # boolean : did we get the whole blast report for this query ? (might be truncated) -
use vars '$infinite';    # high int number - to initialize stuff.
use vars '@lines';       # lines of the last query's blast report.
use vars '$iter';        # counting the iterations for PSI-blast.
use vars '$toprint';     # boolean : shall we print that ?
use vars '$bug_warning'; # print warning if parser detects inconsistency
use vars '%iter_inclusion';	# key: sequence-id, value: first inclusion in profile
use vars '$current_iteration';	# current iteration (PSI-BLAST)
use vars '$incremental_offset';	# offset for incremental PSI-BLAST runs
######################################################################
# selection options.
######################################################################
use Getopt::Long;

my %opt = ();

my %legend = (
    "blastx!" =>  "consider input as blastx output.",
    "tblastn!" => "consider input as tblastn output. (not yet implemented)",
    "blastn!" =>  "consider input as blastn output (will show strand information)",
    "ends!" =>    "report end positions as well (default : only start position are reported)",
    "its!" =>	 "report psi-blast iteration",
    "inc!"      => "result is from incremental run",
    "frame!"    => "when -blastx : use frame information of the HSP and report start/end positions relative to the translated protein in the frame, rather than position in the DNA.",
    "dna!"      => "when -blastx -dna report the position in the DNA as well (columns after alignment strings) (and the frame also)",
    "wu!"       => "consider as wash-U blast output (relevant only for the score)",
    "bits!"     => "get the bits score.",
    "len!"      => "print length of query and sbjct.",
    "h!"        => "this help message.",
    "noheader!" => "do not put a header in the output.",
    "table!"    => "table output",
    "log!"      => "convert evalues to logs.",
    "all!"      => "report all iterations.",
    "cov!"      => "print coverage of query and sbjct.",
    "zero!"     => "use zero based/open-closed coordinates.",
);
	      

GetOptions(\%opt, "help!", keys %legend);

if($opt{'h'} || $opt{'help'}) {
    Usage();
    exit(-1);
}

if($opt{'dna'}) { $opt{'frame'} = 1; }

my $o_dna   = $opt{'dna'};
my $o_frame = $opt{'frame'};
my $o_ends  = $opt{'ends'};
my $o_wu    = $opt{'wu'};
my $o_bits  = $opt{'bits'};
my $o_its   = $opt{'its'};
my $o_table = $opt{'table'};
my $o_log   = $opt{'log'};
my $o_all   = $opt{'all'};
my $o_inc   = $opt{'inc'};
my $o_len   = $opt{'len'};
my $o_cov   = $opt{'cov'};
my $o_zero  = $opt{'zero'};
my $query;
my $current_iteration;
my $incremental_offset = 0;

# minimum evalue
my $MINE = "1e-200";
######################################################################

$infinite=1e42;

$bug_warning='!';

$began = -1;

&printformat;

# read blast output from STDIN; parse one chunk at a time
# note: incomplete entries are not echoed.
while(<STDIN>) 
{ 

    if ($o_inc && /Incremental PSIBLAST: starting from iteration (\d+)/) 
    {
	$incremental_offset = $1;
    }

    if(/^Query=\s*(\S+)/) 
    { # start of chunk
	if ($query)
	{
	    ## call rdb_parser with current_chunk
	    rdb_parser();
	    undef(@lines);
	}
	
	## setup a new chunk
	$query = $1;
	# read consequent lines until an empty line - catch the sequence length on the way.
	$query_length = 0;
	while(<STDIN>) 
	{
	    if(/([0-9,]+)\s+letters/) 
	    {
		$query_length = $1 ;
	    }
	    # for blast+
	    elsif(/Length=(\d+)/)
	    {
		$query_length = $1;
	    }
	    ## explicitely test for an empty line
	    ## last if (/^\s+$/);
	    ## changed for blastp+
	    last if (/Sequences producing significant alignments:/);
	    if(/^Database:/) { $began = 1; last; }
	}
	$query_length =~ s/,//;
	print ">$query\t$query_length\n" unless ($o_table); # pseudo-fasta format

	# start with 0, as first call to rdb_parser is empty;
	$current_iteration=$incremental_offset;
	# cue ^Sequences producing significant ... missing if no hits found
    }

    # take last PSI-Blast iteration; only one if Blast
    if(/^Sequences producing significant alignments:/) 
    { 
	# parse lines quickly update inclusion hash
	rdb_parser() if ($o_all);	# output every iteration
	$current_iteration++;
	undef(@lines);
    }

    push(@lines,$_);

    # end of file: 
    if(/Posted date/) { 
	$began = 0;
    }
}

rdb_parser() if ($query);

exit($began);

sub printformat {
    my($extracolumns,$extracomment);
    
    $extracolumns = ($opt{'blastx'} ? 
		     "\tframe" .
		     "\tquerystart" . ($o_ends ? "\tqueryend" : "") .
		     "\tsbjctstart" . ($o_ends ? " \tsbjctend" : "") . 
		     "\tqueryali\tsbjctali" . 
		     ($o_dna ? 
		      "\tdnaquerystart" . ($o_ends ? "\tdnaqueryend" : "") :
		      "") 
		     :
		     ($opt{'blastn'} ? "\tstrand" 
		      : ($opt{'tblastn'} ? "\tstrand" :
		      "")) .
		     "\tquerystart" . ($o_ends ? "\tqueryend" : "") .
		     "\tsbjctstart" . ($o_ends ? " \tsbjctend" : "") . 
		     "\tqueryali\tsbjctali" );
;
    $extracomment = $opt{'blastn'} ? "and strand is + for plus, and - for minus (strand on the hit)" : "" ;
    
    if (!$opt{'noheader'})
    {
	if (${o_table})
	{
	    printf "qId\tsId" .
		($opt{'tblastn'} ? "\tstrand" : "") . 
		"\teValue\tqStart" . ($o_ends ? "\tqEnd" : "") . "\tqAli" .
		"\tsStart" . ($o_ends ? "\tsEnd" : "") . "\tsAli" .
		"\tscore\tpide" . 
		($o_its ? "\titeration" : "") . 
		( ($opt{'tblastn'} || $o_len) ? "\tqLen\tsLen" : "") .
		( ($opt{'tblastn'} || $o_cov) ? "\tqCov\tsCov" : "" ) .
	      ($o_bits ? "\tbits" : "") . 
	      "\n";
	}
	else
	{
	    printf <<EOF;
\# Blast output parsed into pseudofastaformat by 
\# @ARGV
\#
\# >queryhid\tquerylength
\# sbjcthid\tsbjctlength\tpide\tscore\tevalue$extracolumns
\#
\# where ali=+<number of match-or-insert states>-<number of delete states>
\# $extracomment
EOF
	}
    }
}

sub rdb_parser {

    my($pid,$name,$len2,$score,$evalue,$from1,$from2,$to1,$to2);
    my($sign,$n,$oldpid,$oldname,$oldlen2,@ali1,@ali2);
    my($lali,$pide,$frame,$gaps,$strand);
    
    my($invquery,$invsbjct) = (0,0);
    my($blastn);
    my($one,$two,$three,$four);

    $blastn = $opt{'blastn'};
    $toprint = 0; 
    foreach(@lines) { 

	if($_ eq " ***** No hits found ******\n") {
	    print "No-hits-found\n" unless ($o_table);
	    last;
	}
	# print $_;
	if(/^>/) {
	    ($pid,$name)=/^>\s*[\/\:]*(\S+)\s+(\S{0,1}.*)\s*$/;
#	} elsif(/^\s+Length =\s+(\d+)/) {
#\	    $len2=$1;
	} elsif(/^\s*Length\s*=\s*(\d+)/) {
	    $len2=$1;
	    # a bitscore can look like this: 1.058e+04 
	} elsif( /Score =\s+([\d\.e\+\-]+)\s+bits\s+\(([\d\.]+).* Expect\S* =\s+(\S+)/ ) 
	{
	    $one = $1; $two = $2; $three=$3;
	    &printold($oldpid,$oldlen2,$pide,$score,$evalue,$frame,$from1,$to1,$from2,$to2,\@ali1,\@ali2,$strand,$bitscore);
	    $toprint=1;
	    $oldpid=$pid; $oldname=$name; $oldlen2=$len2;
	    $score=$two;
	    $bitscore=$one;
	    $evalue=$three; $evalue =~ s!,$!!;
	    $from1=$infinite; $from2=$infinite; 
	    $to1 = -$infinite; $to2 = -$infinite;
	    $frame = undef;
	    undef(@ali1); undef(@ali2);
	} elsif($blastn && m,Strand\s+=\s+(\S+)\s+/\s+(\S+),) {
	    $strand = ( $1 eq 'Minus' ? "-" : "+") . ($2 eq 'Minus' ? "-" : "+");
	    $invquery = ($1 eq 'Minus'); if($invquery) { $from1 = -$infinite; $to1 = $infinite; }
	    $invsbjct = ($2 eq 'Minus'); if($invsbjct) { $from2 = -$infinite; $to2 = $infinite; }
	} elsif(m,^Query[\:\s]+(\d+)\s*([\w\-\*]+)\s+(\d+)\s*$,) {
	    # Query   123  --AEAA--- 456
	    if((! $invquery && $1<$from1) || ($invquery && $from1<$1)) { $from1=$1; }
	    if((!$invquery && $3>$to1) || ($invquery && $3>$to1)) { $to1=$3; }
	    $n=0; $sign='+';
	    foreach(split(//,$2)) {
		if(/\-/) { 
		    if($sign eq '+') { 
			push(@ali1,$sign.$n);
			$sign='-'; 
			$n=0;
		    }
		} else {
		    if($sign eq '-') { 
			push(@ali1,$sign.$n);
			$sign='+'; 
			$n=0;
		    }
		}
		$n++;
	    }
	    push(@ali1,$sign.$n);
	} elsif(m,^Query[\:\s]+\s*([\-]+)\s*$,) {
	    # completely empty query lines
	    # query    -------   
	    push(@ali1,'-'.length($1));
	} elsif(/^Sbjct[\:\s]+(\d+)\s*([\w\*\-]+)\s+(\d+)\s*$/) {
	    if((!$invsbjct && $1<$from2) || ($invsbjct && $1>$from2)) { $from2=$1; }
	    if( (!$invsbjct && $3>$to2) || ($invsbjct && $3<$to2)) { $to2=$3; }
	    $n=0; $sign='+';
	    foreach(split(//,$2)) {
		if(/\-/) { 
		    if($sign eq '+') { 
			push(@ali2,$sign.$n);
			$sign='-'; 
			$n=0;
		    }
		} else {
		    if($sign eq '-') { 
			push(@ali2,$sign.$n);
			$sign='+'; 
			$n=0;
		    }
		}
		$n++;
	    }
	    push(@ali2,$sign.$n);
	} elsif(m,^Sbjct[\:\s]+\s*([\-]+)\s*$,) {
	    # completely empty sbjct lines
	    # sbjct    ------- 
	    push(@ali2,'-'.length($1));
	}
	if(/Identities =\s+\d+\/(\d+)\s+\((\d+)\%/) {
	    $lali = $1;
	    $pide = $2;
	}
	if(m,^\s+Frame = ([\+-]\d),) {
	    $frame = $1;
	    $invquery = ($frame =~ m/^-/ && $opt{'blastx'}) ;
	    $invsbjct = ($frame =~ m/^-/ && $opt{'tblastn'});
	    if($invquery) { $from1 = -$infinite; $to1 = $infinite; }
	    if($invsbjct) { $from2 = -$infinite; $to2 = $infinite; }
	}
    }
    
    &printold($oldpid,$oldlen2,$pide,$score,$evalue,$frame,$from1,$to1,$from2,$to2,\@ali1,\@ali2,$strand,$bitscore);
    $strand = undef;
}


sub printold {
    my($hitid,$hitlen,$pide,$score,$evalue,$frame,$qstart,$qend,$sstart,$sEnd,$ali1ref,$ali2ref,$strand,$bitscore) = @_;
    return unless $toprint;
    my ($ali1)=&compressali(@{$ali1ref});
    my ($ali2)=&compressali(@{$ali2ref});
    
    $evalue =~ s/^e/1e/; # e-162 => 1e-162
    $evalue = $MINE if ( $evalue < $MINE );             
    $evalue = log( $evalue) if ($o_log);


    my($qdnastart,$qdnaend);
    if($opt{'blastx'}) {
	if($o_frame) {
	    my($x);
	    $qdnastart = $qstart;
	    $qdnaend   = $qend;
	    if($frame >0) {
		$qstart = ($qdnastart - $frame + 3 )/3;
		$qend   = ($qdnaend +1 -$frame)/3;
	    }elsif($frame<0) {
		$qstart = ($query_length + $frame +1  - $qdnastart)/3 + 1;
		$qend   = ($query_length + $frame +1  - $qdnaend +1)/3;
	    }
	}
    }elsif($opt{'tblastn'}) {
	if($frame<0) {
	    $sstart = ($hitlen - $sstart + 1);
	    $sEnd   = ($hitlen - $sEnd + 1);
	    $strand = "-";
	} else {
	    $strand = "+";
	}
    }

    # convert to zero based coords
    if ($o_zero)
    {
	$qstart -= 1;
	$sstart -= 1;
    }

    my $st;
    if (!$o_table) 
    {
	$st= "$hitid\t$hitlen\t$pide\t$score\t$evalue";
	$st .= ($opt{'blastx'} ? 
		("\t$frame" .
		 "\t$qstart" . ($o_ends ? "\t$qend" : "") .
		 "\t$sstart" . ($o_ends ? " \t$sEnd" : "") . 
		 "\t$ali1\t$ali2" . 
		 ($o_dna ? "\t$qdnastart" . ($o_ends ? "\t$qdnaend" : "") : "")
		 ) 
		:
		($opt{'blastn'} ? "\t$strand" 
		 : 
		 ($opt{'tblastn'} ? "\t$strand" 
		 :
		 "")) .
		"\t$qstart" . ($o_ends ? "\t$qend" : "") .
		"\t$sstart" . ($o_ends ? " \t$sEnd" : "") . 
		"\t$ali1\t$ali2" );
	$st .= "\t$current_iteration" if ($o_its);
	
    } 
    else 
    {
	$st = "$query\t$hitid" .
	    ($opt{'tblastn'} ? "\t$strand" : "") . 
	    "\t$evalue\t$qstart" . ($o_ends ? "\t$qend" : "") . "\t$ali1" .
	      "\t$sstart" . ($o_ends ? "\t$sEnd" : "") . "\t$ali2" .
	      "\t$score\t$pide" . 
	      ($o_its ? "\t$current_iteration" : "") . 
	      ( ($opt{'tblastn'} || $o_len) ? "\t${query_length}\t$hitlen" : "") .
	      ( ($opt{'tblastn'} || $o_cov) ? 
		sprintf("\t%5.2f\t%5.2f", 
			($qend - $qstart + 1) * 100 / $query_length, 
			($sEnd - $sstart + 1) * 100 / $hitlen) : "") .
	      ($o_bits ? "\t$bitscore" : "");
    }
    
    my ($ali11,$ali22) = ($ali1,$ali2);
    $ali11 =~ s/-/+/g; $ali22 =~ s/-/+/g;
    my $lali1=eval "0$ali11";
    my $lali2=eval "0$ali22";
    if( $lali1 != $lali2 ) 
    {
	warn "problem with ali1/ali2 for query=$qstart-$qend ; sbjct=$sstart-$ssend ; $lali1 != $lali2; on $query versus $hitid";
	$st = "$bug_warning$st";
    }
   
    print "$st\n";

    $toprint=0; undef(@{$ali1ref}); undef(@{$ali2ref});
}

sub compressali {
    local(@ali1)=@_;
    local($ali1,$n,$sign,$i);
    $ali1='';
    $_=$ali1[0]; $n=$_; ($sign)=/^([+-])/;
    $i=0; 
    while($i<$#ali1) {
	while($ali1[$i+1]=~/[$sign]/ || $ali1[$i+1] eq '+0') { $i++; $n+=$ali1[$i]; }
	$n=~s/[\-\+]//g; $ali1.=$sign.$n; 
	$i++;
	$_=$ali1[$i]; $n=$_;($sign)=/^([+-])/;
    }	
    $n=~s/[\-\+]//g; $ali1.=$sign.$n; 
    return($ali1);
}

sub Usage {
    print <<"EofMsg";
$0 : parse various types of blast output. 
     (possibility a stream of blast outputs)

$0 [-blastx|-tblastn|-blastn] [-ends] [-frame] [-dna] 
    [-wu] [-bits] [-h|-help]

By default : will consider ncbi blast output similar to blastp.
             will detect output is from wu-blast if can
	     recognize wu-blast header (BLAST..WashU first line).

EofMsg

    my $k;
    foreach $k ('blastx','tblastn','blastn',
	       'ends','frame','dna','wu','bits','h') {
	print "-$k : " . $legend{"$k!"} . "\n";
    }

}
