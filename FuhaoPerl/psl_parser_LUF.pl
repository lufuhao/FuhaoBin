#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;



###HELP################################
use constant USAGE=><<END;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: LUFUHAO20140321

		
Options:
	--help/-h
		Print this help/usage;
	--input|-in <FILE>
		[necessary] Input fasta file to extract from 
	--align_length/-l <INT>
		[Optional] The minimum alignment length you want to 
		output (default:100);
	--align_retio/-r
	--output/-o <FILE>
		[Optional] New file for the extracted sequences
	--verbose
		[Optional] Detailed output for trouble-shooting;
	--version/-v
		[Optional] Print current SCRIPT version;

Example:
	perl $0 --input my.fasta --min_length 200 \
		-s id01,id02,id03 \
		--seqid_file myID.file
		--output Extracted.my.fa
		

Author:
	Fu-Hao Lu
	Post-Doctoral Scientist in Micheal Bevan laboratory
	Cell and Developmental Department, John Innes Centre
	Norwich NR4 7UH, United Kingdom
	E-mail: Fu-Hao.Lu\@jic.ac.uk
END
###HELP ends##########################
die USAGE unless @ARGV;


###Receving parameter#################
our ($help, $input, $aln_length, $aln_ratio, $output, $verbose, $ver);

GetOptions(
	"help|h!" => \$help,
	"input|in=s" => \$input,
	"align_length|l:i" => \$aln_length,
	"align_ratio|r:f" => \$aln_ratio,
	"output|o:s" => \$output,
	"verbose!" => \$verbose,
	"version|v!" => \$ver) or die USAGE;
	
($help or $ver) and die USAGE;



###Default and length
$aln_length=0 unless (defined $aln_length);
$aln_ratio=0 unless (defined $aln_ratio);
our ($input_basename, $input_base_no_ext, $sum_output);
($input_basename=$input)=~ s/.*\///s;
($input_base_no_ext=$input_basename)=~s/^(.*)\.\w+$/$1/;
our $output_table=$input_base_no_ext.".detailed";
unlink("$output_table");
our $output_list=$input_base_no_ext.".IDlist";
unlink("$output_list");
print "Input: $input\nAlign_length>=$aln_length\nAlign_ratio>=$aln_ratio\nOutput: $output_table\n$output_list\n" if defined $verbose;


open (INPUT, "$input") || die "Can not input\n";
open (OUTDETAIL, ">>$output_table") || die "Can not output to Detail\n";
print OUTDETAIL "Query\tSubject\tQratio_aligned\tBlock_count\nQalign_length\n";
open (OUTLIST, ">>$output_list") || die "Can not out put to list\n";

our %non_redundant=();
<INPUT>;<INPUT>;<INPUT>;<INPUT>;<INPUT>;
while (my $input_line=<INPUT>) {
	chomp $input_line;
	my @input_arr=();
	@input_arr=split(/\t/, $input_line);
	my $total_ratio=0;
	$total_ratio=($input_arr[12]-$input_arr[11])/$input_arr[10];
	my $align_line='';
	$align_line=$input_arr[18];
	chomp $align_line;
	chop $align_line;
	my @align_arr=();
	@align_arr=split(/,/, $align_line);
	my $align_length=0;
	foreach (@align_arr) {
		$align_length+=$_;
	}
	my $align_ratio=0;
	$align_ratio=$align_length/$input_arr[10];
	print OUTDETAIL "$input_arr[9]\t$input_arr[13]\t$total_ratio\t$input_arr[17]\t$align_length\t$align_ratio\t$align_line\n";
	
	if ($align_length>=$aln_length and $align_ratio>=$aln_ratio and !(exists $non_redundant{$input_arr[9]})) {
		print OUTLIST $input_arr[9]."\n";
		$non_redundant{$input_arr[9]}++;
		print "ID:$input_arr[9]\tTotal_length: $input_arr[10]\tTotal_ratio: $total_ratio\tAln_length: $align_length\tAln_ratio: $align_ratio\n" if defined $verbose;
	}
}
close INPUT;
close OUTDETAIL;
close OUTLIST;



