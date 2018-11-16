#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FindBin qw($Bin);
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa --function Int [Options]
Version: LUFUHAO20150311

Requirements:
	Programs:
	Modiles: Scalar::Util, Cwd, Getopt::Long, FindBin

Descriptions:
	Determine the insert size given pairs of seqing data by
	mapping them to a reference.

Options:
	--help|-h
		Print this help/usage;
	--input|-i	<PASA_GFF3>
		[Msg] 
	--output|-o	<File>
		[Opt] 
	--function|-f <Int>
		[Opt] Only 1 value as follows:
			1 - Check identity
			2 - Check exon numbers
			3 - Check overlap between two genes
			4 - Add gene/mRNA features
	--1identity_file	<File>
		[Opt] Output from 'calc_pc_id_between_seqs.pl',
		Line format: ID1 ID2 Identity_percentage
	--1maxidentity	<Int>
		[Opt] Will report and select between sequences 
		with similarity >= Thisvalue. Default: 70
	--1autoselect
		[Opt] 
	--2exonnum	<Int[-Int]>
		[Opt] 
	--verbose
		Detailed output for trouble-shooting;
	--version|v!
		Print current SCRIPT version;

Example:
	#1identity_file
	perl ~/perlscript/Translate_DNA_6frames_luf.pl -i PASA.cDNA.fa -lo PASA.PP.longestORF.fa
	calc_pc_id_between_seqs.pl PASA.PP.longestORF.fa [1identity_file] . ggsearch
	cat percentid | perl -ne '\$line=<>;chomp \$line;\@arr=split(/\\s+/,\$line);if (\$arr[2]>=70){print \$line."\\n";}'
	
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
my ($input, $output, $function);
my ($check_identity, $file_identity, $maxidentity, $autoselect);
my ($exonfilter, $exonnum);
my $checkoverlap;
my $check_mrna;



GetOptions(
	"help|h!" => \$help,
	"input|i:s" => \$input,
	"output|o:s" => \$output,
	"function|f:i" => \$function,
	"1identity_file:s" => \$file_identity,
	"1maxidentity:i" => \$maxidentity,
	"1autoselect!" => \$autoselect,
	"2exonnum:s" => \$exonnum,
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$version) or die USAGE;
($help or $version) and die USAGE;



### Configure #######################################################
my $feature1='ID';
my $feature2='Parent';
my $index='ID';



### Defaults ########################################################
die "MainError: Invalid input\n" unless (defined $input and -s $input);
$check_identity=0;
$exonfilter=0;
$checkoverlap=0;
$check_mrna=0;
if (defined $function) {
	if ($function==1) {
		$check_identity=1;
	}
	elsif ($function==2) {
		$exonfilter=1;
	}
	elsif ($function==3) {
		$checkoverlap=1;
	}
	elsif ($function==4) {
		$check_mrna=1;
	}
	else {
		die "MainError: unknown function $function\n";
	}
}
else {
	die "Need to specify which function to use by --function\n";
}



$exonnum='1-99999' unless (defined $exonnum);
my ($min_numexon, $max_numexon)=(1, 99999);
if (defined $exonnum) {
	if ($exonnum=~/^(\d+)-(\d+)$/) {
		$min_numexon=$1;
		$max_numexon=$2;
	}
	elsif ($exonnum=~/^(\d+)$/) {
		$min_numexon=$1;
	}
}
$maxidentity=70 unless (defined $maxidentity);
$autoselect=0 unless (defined $autoselect);
if ($debug) {
	print "\n\n\n###Function activated ###\n";
	print "1. identity check\n" if ($check_identity);
	print "2. Exon number filter\n" if ($exonfilter);
	print "3. Overlap checking\n" if ($checkoverlap);
	print "4. add gene/mRNA feature\n" if ($check_mrna);
}
my $func_sum=$check_identity+$exonfilter+$checkoverlap+$check_mrna;
if ($func_sum != 1) {
	die "MainError: choose one function $func_sum only: 1check_identity/2exonfilter/3check_overlap/4add_mrna\n";
}



### input and output ################################################
my $input_prefix=&RetrvNoExt($input);
if ($check_mrna) {
	$output=$input_prefix.'.WithmRNA.gff3' unless (defined $output);
}
my %id70targets=();
my %tocheckIDsByID70=();
if ($check_identity) {
	die "MainError: Invalid --1identity_file: $file_identity\n" unless (defined $file_identity and -s $file_identity);
	$output=$input_prefix.'.NoID70.gff3' unless (defined $output);
	open (GFFIN3, "<$file_identity") || die "MainError: open file_identity: $file_identity error\n";
	my $num_line2=0;
	my $num_line3=0;
	print "\n\n\n### Sequence identity >= $maxidentity ###" if ($verbose);
	while (my $line=<GFFIN3>) {
		$num_line2++;
		next if ($line=~/^#/);
		chomp $line;
		my @arr=split(/\s+/, $line);
		die "MainError file_identity $file_identity at line $num_line2. Seems NumCol !=3:\n$line\n" unless (scalar(@arr)==3);
		if ($arr[2] >= $maxidentity) {
			$num_line3++;
			print $arr[0]."\t".$arr[1]."\t".$arr[2]."\n" if ($verbose);
			@{$id70targets{$num_line3}}=($arr[0], $arr[1]);
			$tocheckIDsByID70{$arr[0]}++;
			$tocheckIDsByID70{$arr[1]}++;
		}
	}
	close GFFIN3;
}
if ($exonfilter) {
	$output=$input_prefix.'.filter.gff3' unless (defined $output);
}
unlink $output if (defined $output and -e $output);



### Main ############################################################
open (GFFIN, "<$input") || die "MainError: $input open error\n";
my %feature=();
my %target2group=();
my %numexons=();
my %lengthexons=();
my %mrna_startend=();
my $last_end=0;
my $num_lines=0;
while (my $line=<GFFIN>) {
	$num_lines++;
	next if ($line=~/^#/);
	chomp $line;
#0		1		2			3start	4.end	5score	6strand	7frame	8atrributes
#chr	source	feature		113		672		.		+		.		ID=align_8652;Target=asmbl_1 1 560 +

	my @arr1=split (/\t/, $line);
	$last_end=0 unless (exists $feature{$arr1[0]});
	my @attributes=split(/;/, $arr1[8]);
	my ($attr_id, $attr_target)=('', '');
	foreach my $ind_attr (@attributes) {
		if ($ind_attr=~/$feature1=(.*)$/) {
			$attr_id=$1;
		}
		elsif ($ind_attr=~/$feature2=(.*)$/) {
			$attr_target=$1;
		}
	}
	$attr_id=$arr1[8] if ($attr_id eq '');
	$attr_target=$arr1[8] if ($attr_target eq '');
	(my $attr_target2=$attr_target)=~s/\s+.*$//;
	my $hashindex=$attr_id;
	if ($index eq $attr_target2) {
		$hashindex=$attr_target2;
	}
#	print "HashIndex: $hashindex\n";### For test ###
	if ($check_mrna) {
		if (exists $mrna_startend{$hashindex}){
			${$mrna_startend{$hashindex}}{'end'}=$arr1[4];
		}
		else {
			${$mrna_startend{$hashindex}}{'start'}=$arr1[3];
			${$mrna_startend{$hashindex}}{'end'}=$arr1[4];
		}
	}
	if ($check_identity) {
		if (exists $tocheckIDsByID70{$hashindex}) {
			push (@{$target2group{$hashindex}}, $line);
		}
	}
	$numexons{$hashindex}++;
	$lengthexons{$hashindex}+=abs($arr1[4]-$arr1[3]);
	next if ($checkoverlap==0);
#	print $num_lines."\t".$arr1[3]."\t".$last_end."\n";### For test ###
	if ($arr1[3] <= $last_end) {
		my $overlap_started=0;
#		print "\n\n\noverlap\n";### For test ###
		foreach my $ind_numline (sort {$a<=>$b} (keys %{$feature{$arr1[0]}})) {
#			print $arr1[3]."\t".${${$feature{$arr1[0]}}{$ind_numline}}{'start'}."\t".${${$feature{$arr1[0]}}{$ind_numline}}{'end'}."\n";### For test ###
			if ($overlap_started==0 and $arr1[3] <=${${$feature{$arr1[0]}}{$ind_numline}}{'end'}) {
#				print "overlap: yes\n";### For test ###
				print "OverlapStart\t".$num_lines."\t".$arr1[3]."\t".$arr1[4]."\t".$attr_id."\t".$attr_target2."\tvs\t".$ind_numline."\t".${${$feature{$arr1[0]}}{$ind_numline}}{'start'}."\t".${${$feature{$arr1[0]}}{$ind_numline}}{'end'}."\t".${${$feature{$arr1[0]}}{$ind_numline}}{'ID'}."\t".${${$feature{$arr1[0]}}{$ind_numline}}{'target'}."\n";
				$overlap_started=1;
				next;
			}
#			else {
#				print "Ignore1...\n";### For test ###
#			}
			if ($overlap_started==1 and $arr1[4]>=${${$feature{$arr1[0]}}{$ind_numline}}{'start'}) {
				print "OverlapEnd\t".$num_lines."\t".$arr1[3]."\t".$arr1[4]."\t".$attr_id."\t".$attr_target2."\tvs\t".$ind_numline."\t".${${$feature{$arr1[0]}}{$ind_numline}}{'start'}."\t".${${$feature{$arr1[0]}}{$ind_numline}}{'end'}."\t".${${$feature{$arr1[0]}}{$ind_numline}}{'ID'}."\t".${${$feature{$arr1[0]}}{$ind_numline}}{'target'}."\n";
				next;
			}
#			else {
#				print "Ignore2...\n";### For test ###
#			}
		}
#		print "\n\n\n";### For test ###
	}
	else {
		$last_end=$arr1[4];
		${${$feature{$arr1[0]}}{$num_lines}}{'start'}=$arr1[3];
		${${$feature{$arr1[0]}}{$num_lines}}{'end'}=$arr1[4];
		${${$feature{$arr1[0]}}{$num_lines}}{'ID'}=$attr_id;
		${${$feature{$arr1[0]}}{$num_lines}}{'target'}=$attr_target2;
#		print "New\t".$arr1[0]."\t".${${$feature{$arr1[0]}}{$num_lines}}{'start'}."\t".${${$feature{$arr1[0]}}{$num_lines}}{'end'}."\t".${${$feature{$arr1[0]}}{$num_lines}}{'ID'}."\t".${${$feature{$arr1[0]}}{$num_lines}}{'target'}."\n";### For test ###
	}
}
print "\n\n\n### Raw Summmary ###\n";
print "Number of GFF lines: $num_lines\n";
print "Number of chromosomes: ".scalar(keys %feature)."\n";
close GFFIN;
my %numexon_dist=();
foreach (keys %numexons) {
	my $idv_numexon=$numexons{$_};
	$numexon_dist{$idv_numexon}++;
}
print "\n\n\n### SUMMARY: numberExons vs Count ###\nnExons\tcount\n";
foreach (sort {$a<=>$b} (keys %numexon_dist)) {
	print $_."\t".$numexon_dist{$_}."\n";
}



if ($check_identity) {
	my %excludedIDs=();
	my $total_number=scalar(keys %id70targets);
	foreach my $iden_number (sort {$a<=>$b} (keys %id70targets)) {
		my ($id1, $id2)=@{$id70targets{$iden_number}};
		next if (exists $excludedIDs{$id1} or exists $excludedIDs{$id2});
		next if ($id1 eq $id2);
		if ($autoselect) {
			if (exists $numexons{$id1} and exists $numexons{$id2}) {
				if ($numexons{$id1}==$numexons{$id2}) {
					if (exists $lengthexons{$id1} and exists $lengthexons{$id2}) {
						if ($lengthexons{$id1}==$lengthexons{$id2}) {
							$excludedIDs{$id2}++;
						}
						elsif ($lengthexons{$id1}>$lengthexons{$id2}) {
							$excludedIDs{$id2}++;
						}
						elsif ($lengthexons{$id1}<$lengthexons{$id2}) {
							$excludedIDs{$id1}++;
						}
					}
					else {
						$excludedIDs{$id2}++;
					}
				}
				elsif ($numexons{$id1}>$numexons{$id2}) {
					$excludedIDs{$id2}++;
				}
				elsif ($numexons{$id1}<$numexons{$id2}) {
					$excludedIDs{$id1}++;
				}
			}
		}
		else {
			print "\n\n\n#######################\n\n\n";
			if (exists $target2group{$id1} and exists $target2group{$id2}) {
				my $no1_exons=(exists $numexons{$id1}) ? $numexons{$id1} : 'Unknown';
				my $len1_exons=(exists $lengthexons{$id1}) ? $lengthexons{$id1} : 'Unknown';
				print "\n\n### ID1: $id1\tNum_Exons: $no1_exons\tLength_Exons: $len1_exons\tGroup: $iden_number/$total_number ###\n\n";
				if (exists $target2group{$id1}) {
					map {print $_."\n"} @{$target2group{$id1}};
				}
				else {
					print "NaN\n";
				}
				my $no2_exons=(exists $numexons{$id2}) ? $numexons{$id2} : 'Unknown';
				my $len2_exons=(exists $lengthexons{$id2}) ? $lengthexons{$id2} : 'Unknown';
				print "\n\n### ID2: $id2\tNum_Exons: $no2_exons\tLength_Exons: $len2_exons\tGroup: $iden_number/$total_number ###\n\n";
				if (exists $target2group{$id2}) {
					map {print $_."\n"} @{$target2group{$id2}};
				}
				else {
					print "NaN\n";
				}
				print "\n???Please select which one to keep [default 1]: 1 [ID1] or 2 [ID2] or 0 [None] or 3 [Both]";
				my $choice=&AskData();
				chomp $choice;
				if ($choice=~m/^1$/) {
					print "Keep ID1, discard ID2\n";
					$excludedIDs{$id2}++;
				}
				elsif ($choice=~m/^2$/) {
					print "Discard ID1, keep ID2\n";
					$excludedIDs{$id1}++;
				}
				elsif ($choice=~m/^0$/) {
					print "Discard ID1 and ID2\n";
					$excludedIDs{$id1}++;
					$excludedIDs{$id2}++;
				}
				elsif ($choice=~m/^3$/) {
					print "Keep ID1 and ID2\n";
				}
				else {
					print "Default: Keep ID1, discard ID2\n";
					$excludedIDs{$id2}++;
				}
			}
		}
	}
	print "Expected NumIds to discard: $total_number\n";
	print "Selected NumIds to discard: ".scalar(keys %excludedIDs)."\n";
	open (GFFIN3, "<$input") || die "MainError: $input open error\n";
	open (NOIDOUT, ">$output") || die "MainError: $output write error\n";
	while (my $line=<GFFIN3>) {
		$num_lines++;
		next if ($line=~/^#/);
		chomp $line;
		my @arr1=split (/\t/, $line);
		$last_end=0 unless (exists $feature{$arr1[0]});
		my @attributes=split(/;/, $arr1[8]);
		my ($attr_id, $attr_target)=('', '');
		foreach my $ind_attr (@attributes) {
			if ($ind_attr=~/$feature1=(.*)$/) {
				$attr_id=$1;
			}
			elsif ($ind_attr=~/$feature2=(.*)$/) {
				$attr_target=$1;
			}
		}
		$attr_id=$arr1[8] if ($attr_id eq '');
		$attr_target=$arr1[8] if ($attr_target eq '');
		(my $attr_target2=$attr_target)=~s/\s+.*$//;
		my $hashindex=$attr_id;
		if ($index eq $attr_target2) {
			$hashindex=$attr_target2;
		}
		print NOIDOUT $line."\n" unless (exists $excludedIDs{$hashindex});
	}
	close NOIDOUT;
	close GFFIN3;
}



###Filter exon number
if ($exonfilter) {
	print "\n\n\nFilter number of exons: $min_numexon ~ $max_numexon\n";
	my %filter_numexons=();
	open (GFFIN4, "<$input") || die "MainError: $input open error\n";
	open (FILTER, ">$output") || die "MainError: $output write error\n";
	while (my $line=<GFFIN4>) {
		$num_lines++;
		next if ($line=~/^#/);
		chomp $line;
		my @arr1=split (/\t/, $line);
		$last_end=0 unless (exists $feature{$arr1[0]});
		my @attributes=split(/;/, $arr1[8]);
		my ($attr_id, $attr_target)=('', '');
		foreach my $ind_attr (@attributes) {
			if ($ind_attr=~/$feature1=(.*)$/) {
				$attr_id=$1;
			}
			elsif ($ind_attr=~/$feature2=(.*)$/) {
				$attr_target=$1;
			}
		}
		$attr_id=$arr1[8] if ($attr_id eq '');
		$attr_target=$arr1[8] if ($attr_target eq '');
		(my $attr_target2=$attr_target)=~s/\s+.*$//;
		my $hashindex=$attr_id;
		if ($index eq $attr_target2) {
			$hashindex=$attr_target2;
		}
		if (defined $numexons{$hashindex} and $numexons{$hashindex}>=$min_numexon and $numexons{$hashindex}<=$max_numexon) {
			print FILTER $line."\n";
			$filter_numexons{$hashindex}++;
		}
	}
	close FILTER;
	close GFFIN4;
	my %temp_count=();
	foreach (keys %filter_numexons) {
		$temp_count{$filter_numexons{$_}}++;
	}
	print "\n\n\n###After filtering, Number of exons: ###\n";
	map {print $_."\t".$temp_count{$_}."\n"} (sort {$a<=>$b} (keys %temp_count));
}



###Add mRNA and gene feature
if ($check_mrna){
	print "\n\n\n### 4. Add gene and mRNA features ###\n";
	open (GFFIN2, "<$input") || die "MainError: $input open error2\n";
	open (GFFOUT, ">$output") || die "MainError: $output write error\n";
	my %writeIDs=();
	my @print_exon=();
	my $numth_exon=1;
	while (my $line=<GFFIN2>) {
		next if ($line=~/^#/);
		chomp $line;
		my @arr1=split (/\t/, $line);
		my @attributes=split(/;/, $arr1[8]);
		my ($attr_id, $attr_target)=('', '');
		foreach my $ind_attr (@attributes) {
			if ($ind_attr=~/$feature1=(.*)$/) {
				$attr_id=$1;
			}
			elsif ($ind_attr=~/$feature2=(.*)$/) {
				$attr_target=$1;
			}
		}
		$attr_id=$arr1[8] if ($attr_id eq '');
		$attr_target=$arr1[8] if ($attr_target eq '');
		(my $attr_target2=$attr_target)=~s/\s+.*$//;
		my $hashindex=$attr_id;
		if ($index eq $attr_target2) {
			$hashindex=$attr_target2;
		}
		if (exists $writeIDs{$hashindex}) {
#			print GFFOUT $arr1[0]."\tPASA\tCDS\t".$arr1[3]."\t".$arr1[4]."\t".$arr1[5]."\t".$arr1[6]."\t".$arr1[7]."\tID=$hashindex.CDS;Parent=$hashindex.mRNA\n";
			my $exon_print=$arr1[0]."\tPASA\texon\t".$arr1[3]."\t".$arr1[4]."\t".$arr1[5]."\t".$arr1[6]."\t".$arr1[7]."\tID=$hashindex.mRNA.exon$numth_exon;Parent=$hashindex.mRNA\n";
			push (@print_exon, $exon_print);
			$numth_exon++;
		}
		else {
			$writeIDs{$hashindex}=1;
			map {print GFFOUT $_} @print_exon if (scalar(@print_exon)>0);
#			print ${$mrna_startend{$hashindex}}{'start'}."\t".${$mrna_startend{$hashindex}}{'end'}."\n";### For test ###
			print GFFOUT $arr1[0]."\tPASA\tgene\t".${$mrna_startend{$hashindex}}{'start'}."\t".${$mrna_startend{$hashindex}}{'end'}."\t".$arr1[5]."\t".$arr1[6]."\t".$arr1[7]."\tID=$hashindex\n";
			print GFFOUT $arr1[0]."\tPASA\tmRNA\t".${$mrna_startend{$hashindex}}{'start'}."\t".${$mrna_startend{$hashindex}}{'end'}."\t".$arr1[5]."\t".$arr1[6]."\t".$arr1[7]."\tID=$hashindex.mRNA;Parent=$hashindex\n";
#			print GFFOUT $arr1[0]."\tPASA\tCDS\t".$arr1[3]."\t".$arr1[4]."\t".$arr1[5]."\t".$arr1[6]."\t".$arr1[7]."\tID=$hashindex.CDS;Parent=$hashindex.mRNA\n";
			$numth_exon=1;
			my $exon_print=$arr1[0]."\tPASA\texon\t".$arr1[3]."\t".$arr1[4]."\t".$arr1[5]."\t".$arr1[6]."\t".$arr1[7]."\tID=$hashindex.mRNA.exon$numth_exon;Parent=$hashindex.mRNA\n";
			push (@print_exon, $exon_print);
			$numth_exon++;
		}
	}
	map {print GFFOUT $_} @print_exon if (scalar(@print_exon)>0);
	close GFFIN2;
	close GFFOUT;
	print "finished\n";
}
#####################################################################
###                         sub functions                         ###
#####################################################################
### Retrieve filebasename without extension
### &RetrvNoExt(file)
### Global:
### Dependency:
### Note:
sub RetrvNoExt {
	my $RNE_ori=shift @_;
	chomp $RNE_ori;
	my $RNE_new='';
	my $RNE_base='';
	($RNE_base=$RNE_ori)=~ s/.*\///s;
	($RNE_new=$RNE_base)=~s/^(\S+)\.\w+$/$1/;
	return $RNE_new;
}



### Read STDIN from keyboard
###
sub AskData {
	my $ADanswer;
	print "Enter the data before 5 minutes: ";
	eval {
		local $SIG{ALRM} = sub { die "timeout reading from keyboardn" };
		alarm 300;
		$ADanswer = <STDIN>;
		alarm 0;
		chomp $ADanswer;
    };
    if ($@) {
        die $@ if $@ ne "timeout reading from keyboardn";
        $ADanswer = 'No answer given';
    }
    return $ADanswer;
}
