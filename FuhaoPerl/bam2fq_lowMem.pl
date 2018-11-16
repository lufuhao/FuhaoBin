#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: LUFUHAO20150602

Requirements:
	Programs: SamTools
	Modiles: Getopt::Long

Descriptions:reverse
	Determine the insert size given pairs of seqing data by
	mapping them to a reference.

Options:
	--help|-h
		Print this help/usage;
	--mode|-m		<sam/fq>
		[Opt] Extract reads from bam or fastq, default: sam
	--input|-i		<Bam/Sam>
		[Msg] Bam/Sam to extract reads;
	--output|-o		<R1R2.fq>
		[Opt] All the reads = (--forword) + (--reverse)
	--mapped	<0/1/2>
		[Opt] 0=all, default; 1=mapped reads only; 2=unmapped
	--forward|-f	<R1.fq>
		[Opt] Forward reads with SamFlag[64]
	--reverse|-r	<R2.fq>
		[Opt] Reverse reads with SamFlag[128]
	--pair1|-1		<R1.paired.fq>
		[Opt] Paired forward reads from R1.fq vs R2.fq
	--pair2|-2		<R2.paired.fq>
		[Opt] Paired reverse reads from R1.fq vs R2.fq
	--shuffle|-s	<R1R2.paired.fq>
		[Opt] All paired reads = --pair1 + --pair2
	--unpaired1		<R1.unpaired.fq>
		[Opt] Unpaired forward reads = (--forward) - (--pair1)
	--unpaired2		<R2.unpaired.fq>
		[Opt] Unpaired reverse reads = (--reverse) - (--pair2)
	--debug
		[Opt] Keep temporary files to find bugs
	--verbose
		Detailed output for trouble-shooting;
	--version|v!
		Print current SCRIPT version;

Example:
	perl $0 --input

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
our ($help, $input, $output, $verbose, $debug, $ver);
our ($mapped, $mode, $file_forward, $file_reverse, $readlist, $paired_R1, $paired_R2, $unpaired_R1, $unpaired_R2, $file_shuffle);
GetOptions(
	"help|h!" => \$help,
	"mode|m:s" => \$mode,
	"input|i:s" => \$input,
	"output|o:s" => \$output,
	"mapped:i" => \$mapped,
	"forward|f:s" => \$file_forward,
	"reverse|r:s" => \$file_reverse,
	"readlist|l:s" => \$readlist,
	"pair1|1:s" => \$paired_R1,
	"pair2|2:s" => \$paired_R2,
	"shuffle|s:s" => \$file_shuffle,
	"unpaired1:s" => \$unpaired_R1,
	"unpaired2:s" => \$unpaired_R2,
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$ver) or die USAGE;
($help or $ver) and die USAGE;



### Defaults ########################################################
$debug=0 unless (defined $debug);
$mode='sam' unless (defined $mode);
$mapped=0 unless (defined $mapped);
unless ($mode eq 'sam' or $mode eq 'fq') {
	if ($mode =~ /sam/i or $mode =~ /bam/i) {
		$mode='sam';
	}elsif ($mode =~/fastq/i or $mode =~/fastq/i) {
		$mode='fq';
	}
	else {
		die "Error: unknown extraction mode: $mode\n";
	}
}



### input and output ################################################
if ($mode eq 'sam'){
	die "Error: can not find bam/sam input: $input\n" unless (defined $input and -s $input);
	&DeleteIfExist($output) if (defined $output);
}elsif ($mode eq 'fq') {
	die "Error: can not find fastq input specified by --forward and --reverse\n" unless (defined $file_forward and -s $file_forward and -s $file_reverse);
	if (defined $output) {
		print "Warning: --output option is not necessary in 'fq' mode\nIgnoring the output to $output\n";
		undef $output;
	}
	if (defined $unpaired_R1) {
		print "Warning: --unpaired1 option is not necessary in 'fq' mode\nIgnoring the output to $unpaired_R1\n";
		undef $unpaired_R1;
	}
	if (defined $unpaired_R2) {
		print "Warning: --unpaired2 option is not necessary in 'fq' mode\nIgnoring the output to $unpaired_R2\n";
		undef $unpaired_R2;
	}
}
if (defined $paired_R2 or defined $unpaired_R2 or defined $paired_R1 or defined $unpaired_R1 or defined $file_shuffle) {
	$file_forward="temp_R1.fq" if (! defined $file_forward);
	$file_reverse="temp_R2.fq" if (! defined $file_reverse);
}
$paired_R1="temp_R1.paired.fq" if (defined $file_shuffle and ! defined $paired_R1);
$paired_R2="temp_R2.paired.fq" if ((defined $file_shuffle or defined $paired_R1 or defined $unpaired_R1) and ! defined $paired_R2);

if ($mode eq 'sam') {
	&DeleteIfExist($file_forward) if (defined $file_forward);
	&DeleteIfExist($file_reverse) if (defined $file_reverse);
}
&DeleteIfExist($paired_R1) if (defined $paired_R1);
&DeleteIfExist($paired_R2) if (defined $paired_R2);
&DeleteIfExist($unpaired_R1) if (defined $unpaired_R1);
&DeleteIfExist($unpaired_R2) if (defined $unpaired_R2);
&DeleteIfExist($file_shuffle) if (defined $file_shuffle);



### Main ############################################################
if (defined $verbose) {
	print "\n\n#####1. Extracting fastq reads into fastq #####\nSetting:\n";
	print "\t--output:  $output\n" if (defined $output);
	print "\t--forward: $file_forward\n" if (defined $file_forward);
	print "\t--reverse: $file_reverse\n" if (defined $file_reverse);
}
&Bam2Fq($input) if ($mode eq 'sam');
if (defined $paired_R1 or defined $paired_R2 or defined $unpaired_R1 or defined $unpaired_R2) {
	if (-e $file_forward and -e $file_reverse) {
		if (defined $verbose) {
			print "\n\n##### 2. Extract paired reads only #####\nSetting:\n";
			print "\tForward in: $file_forward\n";
			print "\tForward in: $file_reverse\n";
			print "\tPaired F out: $paired_R1\n" if (defined $paired_R1);
			print "\tPaired R out: $paired_R2\n" if (defined $paired_R2);
			print "\tUnpaired F out: $unpaired_R1\n" if (defined $unpaired_R1);
			print "\tUnpaired R out: $unpaired_R2\n" if (defined $unpaired_R2);
		}
		&ExtractPaired($file_forward, $file_reverse);
	}
	else {
		die "Error: can not find raw fastq files to extract paired reads\n";
	}
}
if (defined $file_shuffle) {
	if ($paired_R1 and -e $paired_R2) {
		if (defined $verbose) {
			print "\n\n##### 3. shuffle paired reads into one fastq file #####\nSetting:\n";
			print "Paired F in: $paired_R1\n";
			print "Paired R in: $paired_R2\n";
			print "Shuffle out: $file_shuffle\n";
		}
		&ShuffleFq($paired_R1, $paired_R2);
	}
	else {
		die "Error: Wrong input to call subfunction ShuffleFq\n";
	}
}



#####################################################################
###                         sub functions                         ###
#####################################################################
###Convert bam/sam into fastq files using col0,col9,col10 info
###&Bam2Fq(bamfile);
###Global: $output, $file_forward, $file_reverse, $file_singletons



###Delete that file if exists;
###&DeleteIfExist(file)
sub DeleteIfExist {
	my @DIE_files=@_;
	foreach my $ind_file (@DIE_files) {
		unlink($ind_file) if (defined $ind_file and -e $ind_file);
	}
}



sub Bam2Fq {
	print "##########\n# subfunction Bam2Fq #\n##########\n" if (defined $verbose);
	my $bam_file=shift @_;
	if ($bam_file=~/\S+\.bam/i) {
		print "Input format is BAM\n" if (defined $verbose);
		open (INPUT, "samtools view $bam_file |") || die "Error: can not open BAM input in subfunction Bam2Fq\n";
	}
	elsif ($bam_file=~/\S+\.sam/i) {
		print "Input format is SAM\n" if (defined $verbose);
		open (INPUT, "samtools view -S $bam_file |") || die "Error: can not open SAM input in subfunction Bam2Fq\n";
	}
	else {
		die "Error: can not detect file format\n";
	}
	open (OUTPUT, ">>$output") || die "Error: can not write to file : $output\n" if (defined $output);
	open (FORWARD, ">>$file_forward") || die "Error: can not write to file : $file_forward\n" if (defined $file_forward);
	open (REVERSE, ">>$file_reverse") || die "Error: can not write to file : $file_reverse\n" if (defined $file_reverse);
	while (my $line=<INPUT>) {
		chomp $line;
		next if ($line=~/^@/);
		my @arr=();
		@arr=split(/\t/,$line);
		next if (scalar(@arr)<11);
		if ($mapped==1) {
			next if ($arr[1] & 0x0004);
		}
		elsif ($mapped==2) {
			next unless ($arr[1] & 0x0004);
		}
		my $read_seq='';
		my $read_qual='';
		if ($arr[1] & 0x0010) {
			$read_seq=reverse ($arr[9]);
			$read_seq=~tr/ATCGatcg/TAGCtagc/;
			$read_qual=reverse($arr[10]);
		}else {
			$read_seq=$arr[9];
			$read_qual=$arr[10];
		}		
		if ($arr[1] & 0x0001) {

			if ($arr[1] & 0x0040) {
				print OUTPUT '@'.$arr[0].'/1'."\n".$read_seq."\n+\n".$read_qual."\n" if (defined $output);
				print FORWARD '@'.$arr[0].'/1'."\n".$read_seq."\n+\n".$read_qual."\n" if (defined $file_forward);
				
			}
			elsif ($arr[1] & 0x0080) {
				print OUTPUT '@'.$arr[0].'/2'."\n".$read_seq."\n+\n".$read_qual."\n" if (defined $output);
				print REVERSE '@'.$arr[0].'/2'."\n".$read_seq."\n+\n".$read_qual."\n" if (defined $file_reverse);
			}
			else {
				print "unknown paired SAMFLAG: $arr[1]\n" if (defined $verbose);
			}
		}
		else {
			print OUTPUT '@'.$arr[0]."\n".$read_seq."\n+\n".$read_qual."\n" if (defined $output);
		}
	}
	close OUTPUT if (defined $output);
	close FORWARD if (defined $file_forward);
	close REVERSE if (defined $file_reverse);
}



### extract paired read only
### &ExtractPaired($file_forward, $file_reverse);
### global $paired_R1, $paired_R2, $unpaired_R1, $unpaired_R2, $debug, $readlist
sub ExtractPaired {
    my ($EP_R1, $EP_R2)=@_;
   	print "##########\n# subfunction ExtractPaired #\n##########\n" if (defined $verbose);
    die "Error: subfunction ExtractPaired input error\n" if (! -e $EP_R1 or ! -e $EP_R2);
    my $temp_R1_id='temp_R1.id';
    my $temp_share_id='temp_share.id';
    my $num_lines_R1_id=0;
    my $num_lines_R2_id=0;
    my $num_lines_share_id=0;
    &DeleteIfExist($temp_R1_id, $temp_share_id);
#write R1 ID into file: temp_R1.id
	if (defined $paired_R2 or defined $unpaired_R2 or defined $paired_R1 or $unpaired_R1) {
		if (defined $readlist and -s $readlist) {
			$temp_R1_id=$readlist;
			goto BLOCK2;
		}
		else {
			goto BLOCK1;
		}
	}else {
		goto BLOCK5;
	}
BLOCK1: {
	print "Extract $EP_R1 ID into file: $temp_R1_id\n" if (defined $verbose);
    open (EPR1, "$EP_R1") || die "Error: can not open file in ExtractPaired: \$EP_R1\n";
    open (IDR1, ">>$temp_R1_id") || die "Error: can not write to file in ExtractPaired: \$temp_R1_id\n";
    while (my $line=<EPR1>) {
	if ($line=~ /^@(\S+)\/[12]\s*/) {
	    print IDR1 $1."\n";
	    $num_lines_R1_id++;
	}
    }
    close IDR1;
    close EPR1; }
#parse R2, if ID in temp_R1.id, output to paired R2 and output ID to temp_share.id
#if not, output to unpaired r2
	if (defined $paired_R2 or defined $unpaired_R2 or defined $paired_R1 or $unpaired_R1) {
		goto BLOCK2;
	}else {
		goto BLOCK5;
	}
BLOCK2: {
	if (defined $verbose) {
		print "Extract paired R2 from $EP_R2 to $paired_R2\n" if (defined $paired_R2);
		print "Extract unpaired R2 from $EP_R2 to $unpaired_R2\n" if (defined $unpaired_R2);
	}
    open (EPR2, "<$EP_R2") || die "Error: can not open file in ExtractPaired: \$EP_R2\n";
    if (defined $paired_R2) {
		open (PAIRED2, ">>$paired_R2") || die "Error: can not write to file in ExtractPaired: \$paired_R2\n";
    }
    if (defined $unpaired_R2) {
		open (UNPAIRED2, ">>$unpaired_R2") || die "Error: can not write to file in ExtractPaired: \$unpaired_R2\n";
    }
    open (SHAREID, ">>$temp_share_id") || die "Error: can not write to file in ExtractPaired: \$temp_share_id\n";
    while (my $line2=<EPR2>) {
		my $test_paired=0;
		if ($line2=~/^\@(\S+)\/[12]\s*/) {
			my $cur_id='';
			$cur_id=quotemeta($1);
			$num_lines_R2_id++;
			open (R1IDLIST, "<$temp_R1_id") || die "Error: can not open file in ExtractPaired: \$temp_R1_id\n";
    	    while (my $line3=<R1IDLIST>) {
				chomp $line3;
				($line3=$line3)=~s/^\@//;
				($line3=$line3)=~s/\/[12]\S*\s*//;
				$line3=quotemeta($line3);
				if ($line3 eq $cur_id) {
					$test_paired=1;
#					print $line3." MATCH ".$cur_id."\n";###for test###
					print SHAREID $cur_id, "\n";
					$num_lines_share_id++;
					last;
				}
#				else {###for test###
#					print $line3." NOTMATCH ".$cur_id."\n";###for test###
#				}###for test###
			}
			close R1IDLIST;
		}else{next;}
		if ($test_paired==1 and defined $paired_R2) {
			print PAIRED2 $line2;
			$_=<EPR2>; print PAIRED2 $_;
			$_=<EPR2>; print PAIRED2 $_;
			$_=<EPR2>; print PAIRED2 $_;
		}
		elsif ($test_paired==0 and defined $unpaired_R2) {
			print UNPAIRED2 $line2;
			$_=<EPR2>; print UNPAIRED2 $_;
			$_=<EPR2>; print UNPAIRED2 $_;
			$_=<EPR2>; print UNPAIRED2 $_;
		}
	}
	close SHAREID;
	close UNPAIRED2 if (defined $unpaired_R2);
	close PAIRED2 if (defined $paired_R2);
	close EPR2;}
#parse R1, if ID in temp_share.id, output to paired R1 ($paired_R1.temp); if not, output to unpaired R1
	if (defined $paired_R1 or $unpaired_R1) {
		goto BLOCK3;
	}else {
		goto BLOCK5;
	}
BLOCK3: {
	if (defined $verbose) {
		print "Extract paired R1 from $EP_R1 to dis-ordered $paired_R1.temp\n" if (defined $paired_R1);
		print "Extract unpaired R1 from $EP_R1 to $paired_R1\n" if (defined $paired_R1);
	}
	open (EPR1, "<$EP_R1") || die "Error: can not open file in ExtractPaired: \$EP_R1\n";
	if (defined $paired_R1) {
		open (PAIRED1, ">>$paired_R1.temp") || die "Error: can not write to file in ExtractPaired: \$paired_R1\n";
	}
	if (defined $unpaired_R1) {
		open (UNPAIRED1, ">>$unpaired_R1") || die "Error: can not write to file in ExtractPaired: \$unpaired_R1\n";
	}
	while (my $line4=<EPR1>) {
		my $test_paired2=0;
		if ($line4=~/^\@(\S+)\/[12]\s*/) {
			my $cur_id=quotemeta($1);
#			print $1."\n".$cur_id."\n";###for test###
			open (SHAREIDLIST, "<$temp_share_id") || die "Error: can not open file in ExtractPaired: \$temp_R1_id\n";
			while (my $line5=<SHAREIDLIST>) {
				chomp $line5;
#				$line5=quotemeta($line5);###for test###
				if ($line5 eq $cur_id) {
					$test_paired2=1;
#					print $line5." MATCH ".$cur_id."\n";###for test###
					last;
				}
#				else {###for test###
#					print $line5." NOTMATCH ".$cur_id."\n";###for test###
#				}###for test###
			}
			close SHAREIDLIST;
			if ($test_paired2==1 and defined $paired_R1) {
				print PAIRED1 $line4;
				$_=<EPR1>; print PAIRED1 $_;
				$_=<EPR1>; print PAIRED1 $_;
				$_=<EPR1>; print PAIRED1 $_;
			}
			elsif ($test_paired2==0 and defined $unpaired_R1) {
				print UNPAIRED1 $line4;
				$_=<EPR1>; print UNPAIRED1 $_;
				$_=<EPR1>; print UNPAIRED1 $_;
				$_=<EPR1>; print UNPAIRED1 $_;
			}
		}else{next;}
	}
	close UNPAIRED2 if (defined $unpaired_R1);
	close PAIRED1 if (defined $paired_R1);
	close EPR1;}
#Correct order of R1 to shuffle the paired R1 and R2
	if (defined $paired_R1) {
		goto BLOCK4;
	}else {
		goto BLOCK5;
	}
BLOCK4: {
	if (defined $verbose) {
		print "Re-order $paired_R1.temp into file $paired_R1\n" if (defined $paired_R2);
	}
	open (SHAREIDLIST2, "<$temp_share_id") || die "Error: can not open file in ExtractPaired: \$temp_R1_id\n";
	open (PAIRED1FINAL, ">>$paired_R1") || die "Error: can not write to file in ExtractPaired: \$paired_R1\n";
	while (my $line6=<SHAREIDLIST2>) {
		chomp $line6;
		open (PAIRED1TEMP, "<$paired_R1.temp") || die "Error: can not open file in ExtractPaired: \$paired_R1.temp\n";
		while (my $line7=<PAIRED1TEMP>) {
			if ($line7=~/^\@(\S+)\/[12]\s*/) {
				my $cur_id2=quotemeta($1);
#				print $line6, "\t", $cur_id2,"\n";###for test###
				if ($line6 eq $cur_id2) {
					print PAIRED1FINAL $line7;
					$_=<PAIRED1TEMP>; print PAIRED1FINAL $_;
					$_=<PAIRED1TEMP>; print PAIRED1FINAL $_;
					$_=<PAIRED1TEMP>; print PAIRED1FINAL $_;
					last;
				}
			}else {next;}
		}
		close PAIRED1TEMP;
	}
	close PAIRED1FINAL;
	close SHAREIDLIST2;}
###Delete temporary files
BLOCK5: {
	if ($debug==0){
		print "Deleting temp files: $temp_R1_id, $temp_share_id, $paired_R1.temp\n" if (defined $verbose);
		&DeleteIfExist($temp_R1_id) unless (defined $readlist and -s $readlist);
		&DeleteIfExist($temp_share_id);
		&DeleteIfExist("$paired_R1.temp") if (defined $paired_R1);
		if ($file_forward eq "temp_R1.fq") {
			print "Deleting temp files: $file_forward\n" if (defined $verbose);
			&DeleteIfExist($file_forward);
		}
		if ($file_reverse eq "temp_R2.fq") {
			print "Deleting temp files: $file_forward\n" if (defined $verbose);
			&DeleteIfExist($file_reverse);
		}
		if (defined $paired_R2 and $paired_R2 eq "temp_R2.paired.fq" and ! defined $file_shuffle) {
#		$paired_R1 ne "temp_R1.paired.fq" and 
			print "Deleting temp files: $paired_R2\n" if (defined $verbose);
			&DeleteIfExist($paired_R2);
		}
	}}
}



### shuffle paired fastq sequences into one fastq file
### &ShuffleFq(paired_R1, paired_R2);
###global: $file_shuffle, $debug
sub ShuffleFq{
	my ($SF_r1, $SF_r2)=@_;
	print "##########\n# subfunction ShuffleFq #\n##########\n" if (defined $verbose);
	die "Error: subfunction ShuffleFq input error\n" if (! -e $SF_r1 or ! -e $SF_r2);
	open (SHUFFLE, ">>$file_shuffle") || die "Error failes to write file in ShuffleFq: \$file_shuffle\n";
	open (PAIREDR1SF, "<$SF_r1") || die "Error: failed to open file in ShuffleFq: \$SF_r1\n";
	open (PAIREDR2SF, "<$SF_r2") || die "Error: failed to open file in ShuffleFq: \$SF_r2\n";
	while (my $line=<PAIREDR1SF>){
		print SHUFFLE $line;
		$_=<PAIREDR1SF>; print SHUFFLE $_;
		$_=<PAIREDR1SF>; print SHUFFLE $_;
		$_=<PAIREDR1SF>; print SHUFFLE $_;
		$_=<PAIREDR2SF>; print SHUFFLE $_;
		$_=<PAIREDR2SF>; print SHUFFLE $_;
		$_=<PAIREDR2SF>; print SHUFFLE $_;
		$_=<PAIREDR2SF>; print SHUFFLE $_;
	}
	close PAIREDR2SF;
	close PAIREDR1SF;
	close SHUFFLE;
###Delete temporary files
	if ($debug==0){
		if ($SF_r1 eq "temp_R1.paired.fq") {
			print "Deleting temp_file: $SF_r1\n" if (defined $verbose);
			&DeleteIfExist($SF_r1);
		}
		if ($SF_r2 eq "temp_R2.paired.fq") {
			print "Deleting temp_file: $SF_r1\n" if (defined $verbose);
			&DeleteIfExist($SF_r2);
		}
	}
}
