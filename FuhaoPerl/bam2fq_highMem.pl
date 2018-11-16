#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: LUFUHAO20141118

Requirements:
	Programs: SamTools
	Modiles: Scalar::Util, Cwd, Getopt::Long, FindBin

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
our ($mode, $file_forward, $file_reverse, $readlist, $paired_R1, $paired_R2, $unpaired_R1, $unpaired_R2, $file_shuffle);
GetOptions(
	"help|h!" => \$help,
	"mode|m:s" => \$mode,
	"input|i:s" => \$input,
	"output|o:s" => \$output,
	"forward|f:s" => \$file_forward,
	"reverse|r:s" => \$file_reverse,
	"readlist|l:s" => \$readlist,
	"pair1|1:s" => \$paired_R1,
	"pair2|2:s" => \$paired_R2,
	"shuffle|s:s" => \$file_shuffle,
	"unpaired1|3:s" => \$unpaired_R1,
	"unpaired2|4:s" => \$unpaired_R2,
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$ver) or die USAGE;
($help or $ver) and die USAGE;



### Defaults ########################################################
$debug=0 unless (defined $debug);
$mode='sam' unless (defined $mode);
if (! defined $input) {
	$mode='fq' if (defined $file_forward and -s $file_forward and defined $file_reverse and -s $file_reverse);
}
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
if ($mode eq 'sam') {print "##### Extract mode: Bam2fq #####\n";}
elsif ($mode eq 'fq') {print "##### Extract mode: fq2fq #####\n";}


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
}
if (defined $paired_R2 or defined $unpaired_R2 or defined $paired_R1 or defined $unpaired_R1 or defined $file_shuffle) {
	$file_forward="temp_R1.fq" if (! defined $file_forward);
	$file_reverse="temp_R2.fq" if (! defined $file_reverse);
}
$paired_R1="temp_R1.paired.fq" if (defined $file_shuffle and ! defined $paired_R1);
$paired_R2="temp_R2.paired.fq" if (defined $file_shuffle and ! defined $paired_R2);

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
	open (OUTPUT, ">$output") || die "Error: can not write to file : $output\n" if (defined $output);
	open (FORWARD, ">$file_forward") || die "Error: can not write to file : $file_forward\n" if (defined $file_forward);
	open (REVERSE, ">$file_reverse") || die "Error: can not write to file : $file_reverse\n" if (defined $file_reverse);
	while (my $line=<INPUT>) {
		chomp $line;
		next if ($line=~/^@/);
		my @arr=();
		@arr=split(/\t/,$line);
		next if (scalar(@arr)<12);
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
    my ($num_lines_R1_id, $num_lines_R2_id)=(0, 0);
    my ($num_lines_R1_shared, $num_lines_R2_shared)=(0, 0);
    my ($num_lines_R1_unique, $num_lines_R2_unique)=(0, 0);
    my %share=();
	open (IN1, "<$EP_R1") || die "Error: can not open IN1\n";
	while (my $line1=<IN1>) {
	    my $id1='';
		if ($line1=~/^\@/) {
			if ($line1=~/^\@(\S+)\/[12]/) {
				$id1=$1;
				 
			}elsif ($line1=~/^\@(\S+)\s*/) {
				$id1=$1;
			}else {
				die "Uknown In1 id\n";
			}
			$num_lines_R1_id++;
		}else {next;}
		$share{$id1}=$line1;
		$share{$id1}.=<IN1>;
		$share{$id1}.=<IN1>;
		$share{$id1}.=<IN1>;
	}
	close IN1;
	open (IN2, "<$EP_R2") || die "Error: can not open IN2\n";
	open (PAIR1, ">$paired_R1") || die "Error: can not write PAIR1\n" if (defined $paired_R1);
	open (PAIR2, ">$paired_R2") || die "Error: can not write PAIR2\n" if (defined $paired_R2);
	open (UNIQUE2, ">$unpaired_R2") || die "Error: can not write UNIQUE2\n" if (defined $unpaired_R2);
	while (my $line2=<IN2>) {
		my $id2='';
		if ($line2=~/^\@/) {
			if ($line2=~/^\@(\S+)\/[12]/) {
				$id2=$1;
				 
			}elsif ($line2=~/^\@(\S+)\s*/) {
				$id2=$1;
			}else {
				die "Uknown In1 id\n";
			}
			$num_lines_R2_id++;
		}else {next;}
		if (exists $share{$id2}) {
			if (defined $paired_R2) {
				print PAIR2 $line2;
				$_=<IN2>; print PAIR2 $_;
				$_=<IN2>; print PAIR2 $_;
				$_=<IN2>; print PAIR2 $_;
				$num_lines_R2_shared++;
			}
			if (defined $paired_R1) {
				print PAIR1 $share{$id2};
				$num_lines_R1_shared++;
			}
			delete $share{$id2};
		}elsif (defined $unpaired_R2) {
			print UNIQUE2 $line2;
			$_=<IN2>; print UNIQUE2 $_;
			$_=<IN2>; print UNIQUE2 $_;
			$_=<IN2>; print UNIQUE2 $_;
			$num_lines_R2_unique++;
		}
	}
	close UNIQUE2 if (defined $unpaired_R2);
	close PAIR2 if (defined $paired_R2);
	close PAIR1 if (defined $paired_R1);
	close IN2;
	if (defined $unpaired_R1) {
		open (UNIQUE1, ">$unpaired_R1") || die "Error: can not write UNIQUE2\n";
		foreach (keys %share) {
			print UNIQUE1 $share{$_};
			$num_lines_R1_unique++;
		}
		close UNIQUE1;
	}
	if ($mode eq 'sam' and $debug==0) {
		print "#Clean temp files...\n";
		&DeleteIfExist("temp_R1.fq", "temp_R2.fq")
	}
	print "### Paired reads extraction Summary ###\n";
	print "Total  files1: $num_lines_R1_id\nTotal  files2: $num_lines_R2_id\n";
    print "Total paired1: $num_lines_R1_shared\n" if (defined $paired_R1);
    print "Total paired2: $num_lines_R2_shared\n" if (defined $paired_R2);
    print "Total unique1: $num_lines_R1_unique\n" if (defined $unpaired_R1);
    print "Total unique2: $num_lines_R2_unique\n" if (defined $unpaired_R2);
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
