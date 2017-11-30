#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FuhaoPerl5Lib::CmdKit qw/exec_cmd_return CanRun/;
use FuhaoPerl5Lib::FastaKit qw/IndexFasta ExtractFastaSeqtk NumSeq/;
use FuhaoPerl5Lib::FileKit qw /DeletePath AddFilePath MergeFiles CountLines/;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: LUFUHAO20171110

Requirements:
    Programs: MUMmer v3.23, quickmerge, seqtk
    FuhaoPerl5Lib: CmdKit, FastaKit, FileKit
    Modiles: Cwd, Getopt::Long

Descriptions:
    Determine the insert size given pairs of seqing data by
    mapping them to a reference.

Options:
    --help|-h
        Print this help/usage;
    --input | -i
        Mummer coord file
    --query | -q
        Query fasta
    --reference | -r
        Reference fasta file
    --identity | -d
        Min percentage for delta-filter
    --alignlen | -l
        Min alignment length for delta-filter
    --prefix | -p
        Output prefix
    --verbose
        Detailed output for trouble-shooting;
    --version|v!
        Print current SCRIPT version;

Example:
	perl $0 

Author:
	Fu-Hao Lu
	Post-Doctoral Scientist in Micheal Bevan laboratory
	Cell and Developmental Department, John Innes Centre
	Norwich NR4 7UH, United Kingdom
	E-mail: Fu-Hao.Lu\@jic.ac.uk

EOH
###HELP ends#########################################################
die USAGE unless @ARGV;



###Receving parameter################################################
my ($help, $verbose, $debug, $version);
my ($input, $queryfasta, $referencefasta, $percentage, $alignlen, $prefix);

GetOptions(
	"help|h!" => \$help,
	"input|i:s" => \$input,
	"query|q:s" => \$queryfasta,
	"reference|r:s" => \$referencefasta,
	"identity|d:s" => \$percentage,
	"alignlen|l:i" => \$alignlen,
	"prefix|p:s" => \$prefix,
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$version) or die USAGE;
($help or $version) and die USAGE;



### Defaults ########################################################
$debug=0 unless (defined $debug);
$verbose=0 unless (defined $verbose);
$percentage=0 unless (defined $percentage);
$alignlen=0 unless (defined $alignlen);
$prefix='MyOut' unless (defined $prefix);
my $foldername='TEMP0000000000';
my $curdir = getcwd;
my %rfn2qry=();
my %qry2rfn=();
my %final_rfn2qry=();
my %final_qry2rfn=();
my $nucmer_option=" -c 100 ";
my $showcoords_options=" -r -l -T ";
my $deltafilter_option=" -i 99.00 -r -q ";
my $quickmerge_option="  -hco 5.0 -c 1.5 -l 5000 -ml 500 ";
my $out_assembled="$prefix.assembled.fasta";
$out_assembled=AddFilePath($out_assembled);
my $out_not_assembled="$prefix.not.assembled.query.list";
$out_not_assembled=AddFilePath($out_not_assembled);
my $out_rfn2qry="$prefix.assembled.rfn2qry";
$out_rfn2qry=AddFilePath($out_rfn2qry);
my $out_twice="$prefix.assembled.over-assembled";
$out_twice=AddFilePath($out_twice);
my $total_runs=0;
my $successful_runs=0;



### Commands ########################################################
unless (CanRun("nucmer") and CanRun("delta-filter") and CanRun("show-coords") and CanRun("mummerplot")) {
	die "Error: mummerplot not found\n";
}
unless (CanRun("quickmerge")) {
	die "Error: quickmerge not found\n";
}
unless (CanRun("samtools")) {
	die "Error: samtools not found\n";
}



### input and output ################################################

($input)=AddFilePath($input);
unless (defined $input and -s $input) {
	die "Error: invalid nucmer coord input: $input\n";
}
($queryfasta)=AddFilePath($queryfasta);
unless (defined $queryfasta and -s $queryfasta) {
	die "Error: invalid query fasta\n";
}
($referencefasta)=AddFilePath($referencefasta);
unless (defined $referencefasta and -s $referencefasta) {
	die "Error: invalid reference fasta\n";
}
unless (defined $percentage and $percentage=~/^\d+\.*\d*$/) {
	die "Error: invalid percentage\n";
}
unless (defined $alignlen and $alignlen=~/^\d+$/) {
	die "Error: invalid min alignemnt length\n";
}
unlink $out_assembled if (-e $out_assembled);
unlink $out_not_assembled if (-e $out_not_assembled);
unlink $out_rfn2qry if (-e $out_rfn2qry);
unlink $out_twice if (-e $out_twice);



### Summary
print "### INPUT SUMMARY ###\n";
print "Coord:         $input\n";
print "Query:         $queryfasta\n";
print "Reference:     $referencefasta\n";
print "Percentage:    $percentage\n";
print "AlignLen:      $alignlen\n";
print "Prefix:        $prefix\n";
print "Assembled:     $out_assembled\n";
print "NOT assebleed: $out_not_assembled\n";
print "Final rfn2qry: $out_rfn2qry\n";
print "OverAssembled: $out_twice\n";



### Main ############################################################
open (COORDFILE, "< $input") || die "Error: Can not open nucmer coord input\n";
<COORDFILE>;<COORDFILE>;<COORDFILE>;<COORDFILE>;### ignore first 4 header line
my $linenum=4;
while (my $line=<COORDFILE>) {
	chomp $line;
	$linenum++;
	my @arr=split(/\t/, $line);
#0		1		2		3		4		5		6		7		8		9		10		11		12
#[S1]	[E1]	[S2]	[E2]	[LEN 1]	[LEN 2]	[% IDY]	[LEN R]	[LEN Q]	[COV R]	[COV Q]	[ref]	[query]
	unless (defined $arr[11] and $arr[11]=~/^\S+/ and defined $arr[12] and $arr[12]=~/^\S+/) {
		die "Error: invalid coord line ($linenum)\n$line\n";
	}
	$rfn2qry{$arr[11]}{$arr[12]}++;
	$qry2rfn{$arr[12]}{$arr[11]}++;
}
close COORDFILE;



### Check if one query mapped to >1 reference
if($verbose) {
	foreach my $qry1 (sort keys %qry2rfn) {
		if (scalar (keys %{$qry2rfn{$qry1}}) >1) {
			print STDERR "Warnings: Query hit >1 reference\n";
			print STDERR "    $qry1\n";
			foreach (sort keys %{$qry2rfn{$qry1}}) {
				print STDERR "        -$_\n";
			}
		}
	}
}



unless (-s "$queryfasta.fai") {
	IndexFasta("$queryfasta") || die "Error: can not index query fasta\n";
}
unless (-s "$referencefasta.fai") {
	IndexFasta("$referencefasta") || die "Error: can not index reference fasta\n";
}
open (OUTASSEM, " > $out_assembled") || die "Error: can not write assembled output: $out_assembled\n";
open (OUTRFN2QRY, " > $out_rfn2qry") || die "Error: can not write rfn2qry: $out_rfn2qry\n";

foreach my $ref (sort keys %rfn2qry) {
	$foldername++;
	$total_runs++;
	print "\n\n\n### Serial\t$foldername\tREF_SEQIDS: $ref\n";
	print STDERR "\n\n\n### Serial\t$foldername\tREF_SEQIDS: $ref\n";
	chdir $curdir || die "Error: can not change back to running dir: $curdir\n";
	DeletePath "$curdir/$foldername" if (-d "$curdir/$foldername");
	mkdir "$curdir/$foldername",0766 || die "Error: can not create dir: $foldername";
	chdir "$curdir/$foldername" || die "Error: can not change dir: $curdir/$foldername\n";

### Defined out
	my $rawquerylist="$curdir/$foldername/$foldername.query.list";
	my $rawqueryfa="$curdir/$foldername/$foldername.query.fa";
	my $rawreflist="$curdir/$foldername/$foldername.rfn.list";
	my $rawreffa="$curdir/$foldername/$foldername.rfn.fa";
	my $rawmerged="$curdir/$foldername/$foldername.query.rfn.fa";
	my $rawassemble="$curdir/$foldername/merged.fasta";
	my $newquerylist="$curdir/$foldername/$foldername.query2.list";
	my $newqueryfa="$curdir/$foldername/$foldername.query2.fa";
	my $newreflist="$curdir/$foldername/$foldername.rfn2.list";
	my $newreffa="$curdir/$foldername/$foldername.rfn2.fa";

### Extract raw Reference
	close LISTREF if (defined fileno(LISTREF));
	unless (open(LISTREF, " > $rawreflist")) {
		print STDERR "Warnings: Failed to write raw ref list: $rawreflist\n";
		next;
	}
	print LISTREF $ref, "\n";
	close LISTREF;
	unless (&GetSubSeq($referencefasta, $rawreffa, $rawreflist, 1)) {
		print STDERR "Warnings: failed to extract reference fa: $ref\n";
		next;
	}

### Extract raw query
	close LISTQUERY if (defined fileno(LISTQUERY));
	unless (open (LISTQUERY, " > $rawquerylist")) {
		print STDERR "Warnings: failed to write query list: $ref\n";
		next;
	}
	my $numqry1=0;
	foreach (sort keys %{$rfn2qry{$ref}}) {
		print LISTQUERY $_, "\n";
		$numqry1++;
	}
	close LISTQUERY;
	if ($numqry1>0) {
		unless (&GetSubSeq($queryfasta, $rawqueryfa, $rawquerylist, $numqry1)) {
			print STDERR "Warnings: failed to extract query fa: $ref\n";
			next;
		}
	}
	else {
		print STDERR "Warnings: empty query list: $ref\n";
		next;
	}

### Run mummer
	unless (&RunMummerplot($rawreffa, $rawqueryfa, "$curdir/$foldername/$foldername.run1")) {
		print STDERR "Error: RunMummerplot running failed 1: $ref\n";
		next;
	}
	unless (-s "$curdir/$foldername/$foldername.run1.rq.delta") {
		print STDERR "Error: nucmer output failed: $ref\n";
		next;
	}

### Run quickmerge
	unless (exec_cmd_return("quickmerge -d $curdir/$foldername/$foldername.run1.rq.delta -q $rawqueryfa -r $rawreffa $quickmerge_option > quickmerge.log 2>quickmerge.err")) {
		print STDERR "Error: quickmerge running failed: $ref\n";
		next;
	}
	unless (-s "$curdir/$foldername/merged.fasta") {
		print STDERR "Error: quickmerge output failed: $ref\n";
		next;
	}
	unless (rename $rawassemble,"$curdir/$foldername/$foldername.merged.fasta") {
		print STDERR "Warnings: rename failed: $ref\n";
	}
	else {
		$rawassemble="$curdir/$foldername/$foldername.merged.fasta";
	}

### Test
	if ($verbose) {
		unless (MergeFiles($rawmerged, $rawqueryfa, $rawreffa)) {
			print STDERR "Error: merge files failed\n";
			next;
		}
		if (-s $rawmerged) {
			unless (&RunMummerplot($rawassemble, $rawmerged, "$curdir/$foldername/$foldername.run2")) {
				print STDERR "Error: RunMummerplot running failed 2: $ref\n";
				next;
			}
		}
	}

### Analyze assembly seqids
	IndexFasta("$rawassemble") || die "Error: can not index reference fasta: $rawassemble\n";
	unless (-s "$rawassemble.fai") {
		print STDERR "Error: failed to index fasta: $rawassemble\n";
		next;
	}
	my %mergedfalen=();
	close MERGEDFASTA if (defined fileno(MERGEDFASTA));
	unless (open (MERGEDFASTA, " < $rawassemble.fai")) {
		print STDERR "Error: failed to open merged fasta index: $rawassemble.fai\n";
		next;
	}
	while (my $line=<MERGEDFASTA>) {
		chomp $line;
		my @arr=split(/\t/, $line);
		if (exists $mergedfalen{$arr[0]}) {
			die "Error: suplicated fasta ids: $arr[0]\n";
		}
		else {
			$mergedfalen{$arr[0]}=$arr[1];
		}
	}
	close MERGEDFASTA;

### Get raw query seqids and length
	IndexFasta($rawqueryfa) || die "Error: failed to index query: $rawqueryfa\n";
	unless (-s "$rawqueryfa.fai") {
		print STDERR "Error:  not existing index fasta: $rawqueryfa\n";
		next;
	}	
	my %queryfalen=();
	close QUERYFASTA if (defined fileno(QUERYFASTA));
	unless (open (QUERYFASTA, " < $rawqueryfa.fai")) {
		print STDERR "Error: failed to open query fasta index: $rawqueryfa.fai\n";
		next;
	}
	while (my $line=<QUERYFASTA>) {
		chomp $line;
		my @arr=split(/\t/, $line);
		if (exists $queryfalen{$arr[0]}) {
			die "Error: suplicated fasta ids: $arr[0]\n";
		}
		else {
			$queryfalen{$arr[0]}=$arr[1];
		}
	}
	close QUERYFASTA;

### Tell which raw query is assembled and which not
	my %mergedqryseqids=();
	my %notmergedseqids=();
	
	foreach my $idvseq (sort keys %queryfalen) {
		if (exists $mergedfalen{$idvseq}) {
			unless ($mergedfalen{$idvseq}=~/^\d+$/ and $queryfalen{$idvseq}=~/^\d+$/) {
				die "Error: invalid sequence length\n";
			}
			if ($mergedfalen{$idvseq}==$queryfalen{$idvseq}) {
				$notmergedseqids{$idvseq}++;
				delete $mergedfalen{$idvseq};
			}
			else {
				die "Error: same seqid but different length: $idvseq\n";
			}
		}
		else {
			$mergedqryseqids{$idvseq}++;
		}
	}

### OUTPUT summary
	print "##############################\n";
	print "Assembled TO:\n";
	foreach (sort keys %mergedfalen) {
		print "        $_\n";
	}
	print "Assembled FROM:\n";
	foreach (sort keys %mergedqryseqids) {
		print "        $_\n";
	}
	print "NOTAssembled:\n";
	foreach (sort keys %notmergedseqids) {
		print "        $_\n";
	}
	print "##############################\n";

### Check if there is sequence assembled
	unless (scalar(keys %mergedqryseqids)>0) {
		print STDERR "AssembleInfo: NO sequence assebled from: $ref\n";
		next;
	}
	unless (scalar(keys %queryfalen)>0) {
		print STDERR "AssembleInfo: NO sequence assebled into: $ref\n";
		next;
	}

### Extract Assembled sequences
	if (scalar(keys %notmergedseqids)>0) {
		print "AssembleInfo: partial query assembled: $ref\n";
#		unless (open (FASTAEXTQUERY, " > $newquerylist")) {
#			print STDERR "Error: failed to write list1: $newquerylist\n";
#			next;
#		}
#		foreach (sort keys %mergedqryseqids) {
#			print FASTAEXTQUERY "$_\n";
#		}
#		close FASTAEXTQUERY;
#		unless (open (FASTAEXTQUERY, " > $newquerylist")) {
#			print STDERR "Error: failed to write list1: $newquerylist\n";
#			next;
#		}
#		foreach my $idvseq (sort keys %mergedqryseqids) {
#			print FASTAEXTQUERY $idvseq, "\n";
#		}
#		close FASTAEXTQUERY;
#		unless (-s "$newquerylist") {
#			print STDERR "Error: failed to extract new query list\n";
#			next;
#		}
#		unless (&GetSubSeq($rawqueryfa, $newqueryfa, $newquerylist, scalar(keys %mergedqryseqids))) {
#			print STDERR "Error: failed to extract new query fasta\n";
#			next;
#		}
#		unless (-s $newqueryfa) {
#			print STDERR "Error: empty new query fasta\n";
#			next;
#		}

		unless (open (FASTAEXTREF, " > $newreflist")) {
			print STDERR "Error: failed to write list2: $newreflist\n";
			next;
		}
		foreach my $idvseq (sort keys %mergedfalen) {
			print FASTAEXTREF $idvseq, "\n";
		}
		close FASTAEXTREF;
		unless (-s "$newreflist") {
			print STDERR "Error: failed to extract new ref list\n";
			next;
		}
		unless (&GetSubSeq($rawassemble, $newreffa, $newreflist, scalar(keys %mergedfalen))) {
			print STDERR "Error: failed to extract new ref fasta\n";
			next;
		}
		unless (-s $newreffa) {
			print STDERR "Error: empty new ref fasta\n";
			next;
		}
	}
	else {
		print "AssembleInfo: all query assembled: $ref\n";
		$newqueryfa=$rawqueryfa;
		$newreffa=$rawassemble;
	}

### print Assembled into OUTPUT
	close FINALFASTA if (defined fileno(FINALFASTA));
	unless (open(FINALFASTA," < $newreffa")) {
		print "Error: open final assembly error: $newreffa\n";
		next;
	}
	while (my $line=<FINALFASTA>) {
		print OUTASSEM $line;
	}
	close FINALFASTA;
	print OUTRFN2QRY "REF_SEQIDS\t$ref\tQRY_SEQIDS";
	foreach (keys %mergedqryseqids) {
		$final_rfn2qry{$ref}{$_}++;
		$final_qry2rfn{$_}{$ref}++;
		print OUTRFN2QRY "\t$_";
	}
	print OUTRFN2QRY "\n";
	
	chdir $curdir || die "Error: can not change dir2: $curdir\n";
	DeletePath "$curdir/$foldername" unless ($verbose);
	$successful_runs++;
}
close OUTASSEM;
close OUTRFN2QRY;



open (NOTASSEMBLED, " > $out_not_assembled") || die "Error: can not write Not Assembled query seqid list: $out_not_assembled\n";
foreach my $indseq (sort keys %qry2rfn) {
	unless (exists $final_qry2rfn{$indseq}) {
		print NOTASSEMBLED $indseq, "\n";
	}
}
close NOTASSEMBLED;



open (OVERASSEMBLED, " > $out_twice") || die "Error: can not wirte over-assembled list: $out_twice\n";
foreach my $indseq (sort keys %final_qry2rfn) {
	if (scalar(keys %{$final_qry2rfn{$indseq}})>1) {
		print OVERASSEMBLED "QRY_SEQIDS\t$indseq\tREFSEQIDS";
		foreach (sort keys %{$final_qry2rfn{$indseq}}) {
			print OVERASSEMBLED "\t$_";
		}
		print OVERASSEMBLED "\n";
	}
}
close OVERASSEMBLED;



print "\n\n\n########################################\n";
print "Script running successfully finished\n";
print "Total runs:     $total_runs\n";
print "Successful:     $successful_runs\n";
print "Failed:         ", ($total_runs-$successful_runs), "\n";
print "Total refseqs:  ", scalar(keys %rfn2qry), "\n";
print "  -successful   ", scalar(keys %final_rfn2qry), "\n";
print "Total query:    ", scalar(keys %qry2rfn), "\n";
print "  -successful   ", scalar(keys %final_qry2rfn), "\n";



#####################################################################
###                         sub functions                         ###
#####################################################################
### ReadSam
###&ReadSam(sam,ref, 1/2/3)
###Global: $nucmer_option, $deltafilter_options, $deltafilter_option
###Dependency: exec_cmd_return
###Note:
sub RunMummerplot {
	my ($RMrfn, $RMqry, $RMpfx)=@_;
	
	my $RMsubino="SUB(RunMummerplot)";
	
	unless (exec_cmd_return("nucmer $nucmer_option -prefix $RMpfx $RMrfn $RMqry > /dev/null 2>&1")) {
		print STDERR $RMsubino, "Error: nucmer running failed\n";
		return 0;
	}
	unless (-s "$RMpfx.delta") {
		print STDERR $RMsubino, "Error: nucmer output error\n";
		return 0;
	}
	if (0) {### Test ###
		unless (exec_cmd_return("show-coords  $showcoords_options $RMpfx.delta > $RMpfx.delta.coord 2> /dev/null")) {
			print STDERR $RMsubino, "Error: show-coords running failed\n";
		}
	}
	if (1) {### Test ###
		unless (exec_cmd_return("mummerplot --large --png -p $RMpfx $RMpfx.delta > /dev/null 2>&1")) {
			print STDERR $RMsubino, "Error: mummerplot running failed 1\n";
		}
	}
	if (0) {### Test ###
		unless (exec_cmd_return("delta-filter -1 -i 99.00 $RMpfx.delta > $RMpfx.1to1.delta 2> /dev/null")) {
			print STDERR $RMsubino, "Error: delta-filter running failed 2\n";
		}
		unless (exec_cmd_return("mummerplot --layout --large --png -p $RMpfx.1to1 $RMpfx.1to1.delta > /dev/null 2>&1")) {
			print STDERR $RMsubino, "Error: mummerplot running failed 2\n";
		}
		
	}
	unless (exec_cmd_return("delta-filter $deltafilter_option $RMpfx.delta > $RMpfx.rq.delta 2> /dev/null")) {
		print STDERR $RMsubino, "Error: delta-filter running failed 3\n";
	}
	if ($verbose) {
		unless (exec_cmd_return("show-coords  $showcoords_options $RMpfx.rq.delta > $RMpfx.rq.delta.coord 2> /dev/null")) {
			print STDERR $RMsubino, "Error: show-coords running failed\n";
		}
		unless (exec_cmd_return("mummerplot --layout --large --postscript -p $RMpfx.rq $RMpfx.rq.delta > /dev/null 2>&1")) {
			print STDERR $RMsubino, "Error: mummerplot running failed 3\n";
		}
	}
	return 1;
}




sub GetSubSeq {
	my ($GSSin, $GSSout, $GSSlist, $GSSexp_num)=@_;
	
	my $GSSsubinfo='SUB()';
	
	unless (defined $GSSin and -s $GSSin) {
		print STDERR $GSSsubinfo, "Error: invalid fasta input\n";
		return 0;
	}
	unless (defined $GSSlist and -s $GSSlist) {
		print STDERR $GSSsubinfo, "Error: invalid list\n";
		return 0;
	}
	unless (defined $GSSout) {
		print STDERR $GSSsubinfo, "Error: invalid fasta output\n";
		return 0;
	}
	unlink "$GSSout" if (-s "$GSSout");
	unless (defined $GSSexp_num) {
		$GSSexp_num=CountLines($GSSlist);
	}
	
	unless (ExtractFastaSeqtk($GSSin, $GSSout, $GSSlist, 'seqtk')) {
		print STDERR $GSSsubinfo, "Warnings: failed to extract seq list : $GSSlist\n";
		return 0;
	}
	my $GSSnumseq=NumSeq($GSSout);
	unless ($GSSexp_num==$GSSnumseq and $GSSnumseq>0) {
		print STDERR "Warnings: unexpected seq number : $GSSlist\n";
		return 1;
	}
}
