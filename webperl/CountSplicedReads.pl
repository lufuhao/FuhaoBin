#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;


## SET PATH HERE ##
my $BAMUTIL_BIN = "/usr/bin/";
my $PICARD_BIN = "/usr/bin/";



$BAMUTIL_BIN =~ s/(\/)$//;
$PICARD_BIN =~ s/(\/)$//;

	if(!(scalar(@ARGV))){
	print_help();
	exit;
	}
	
my($inBam,$removeTemp);
$removeTemp= "T";
my $result = GetOptions ("in-bam=s" => \$inBam,
			 "remove-temp=s" => \$removeTemp,
			); ## input bam file

my $fileTag = $inBam;
$fileTag =~ s/(\.bam)$//;

## temporary output files ##
my $splicedAlignBam = $fileTag.".splicedOnly.bam";  ## out
my $splicedAlignSortedBam = $fileTag.".splicedOnly.sorted.bam";  ## out
my $splicedAlignStats = $fileTag.".PicardsSplicedAlign.stats";  ## out
my $logFile = $fileTag.".output.log";


## required executables ##
my $bamutil = $BAMUTIL_BIN."/bam";
my $picardSort  = $PICARD_BIN."/SortSam.jar";
my $picardCollectStats = $PICARD_BIN."/CollectAlignmentSummaryMetrics.jar";


my($spliced);
my %alignStatsHash = get_alignment_stats_hash();  ## alignment hash

## checking input files ##
	if( !( -e $inBam) ){
	die "\n input bam file $inBam does not exist\n";
	}

print "Running bamUtil to collect spliced read alignments\n";

	
`$bamutil findCigars --in $inBam --out $splicedAlignBam --cskip`;

	if(!( -e $splicedAlignBam)){
	die "\n bamutil run failed. Spliced alignment file does not exist\n";
	}

print "\nrunning picard to sort the extracted bam file\n";
`java -Xmx2g -jar $picardSort INPUT=$splicedAlignBam OUTPUT=$splicedAlignSortedBam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate`;

	if(!( -e $splicedAlignSortedBam)){
	die "\n Picard SortSam run failed. Spliced alignment sorted file does not exist\n";
	}
	
print "\nrunning picard to collect alignment metrics\n";
`java -Xmx2g -jar $picardCollectStats INPUT=$splicedAlignSortedBam OUTPUT=$splicedAlignStats ASSUME_SORTED=true`;


get_picard_alignmet_stats($splicedAlignStats,\%alignStatsHash);
$spliced = $alignStatsHash{"total-mapped-all"};  ## spliced reads

print "\nTotal number of spliced reads: $spliced\n\n";

open FH2,">$logFile";  ## output file
print FH2 "Total number of spliced reads: $spliced\n";
close FH2;


	if($removeTemp =~ m/t/i){  ## if removeTemp is true
	`rm $splicedAlignBam`;
	`rm $splicedAlignSortedBam`;
	`rm $splicedAlignStats`;
	}


	###########################
	sub get_alignment_stats_hash{
	my %alignStatsHash = (
		"total-mapped-all" => 0,
		"total-mapped-pf" => 0,
		"first-pair-all" => 0,
		"first-pair-pf" => 0,
		"second-pair-all" => 0,
		"second-pair-pf" => 0,
		);
	return(%alignStatsHash);
	}  ## function ends
	###########################
	sub get_picard_alignmet_stats{
	my($inFile,$hashRef) = @_;
	open FHGetStats,$inFile or die "\n can not open file $inFile\n";
	my($str,@a1);

		while($str = <FHGetStats>){
		$str =~ s/\n//;
		$str =~ s/\r//;
			if( ($str =~ m/^(PAIR)/) || ($str =~ m/^(UNPAIRED)/) ){
			@a1 = split(/\t+/,$str);
			$hashRef->{"total-mapped-all"} = $a1[1];
			$hashRef->{"total-mapped-pf"} = $a1[2];
			}
	
			elsif($str =~ m/^(FIRST)/){
			@a1 = split(/\t+/,$str);
			$hashRef->{"first-pair-all"} = $a1[1];
			$hashRef->{"first-pair-pf"} = $a1[2];
			}
	
			elsif($str =~ m/^(SECOND)/){
			@a1 = split(/\t+/,$str);
			$hashRef->{"second-pair-all"} = $a1[1];
			$hashRef->{"second-pair-pf"} = $a1[2];
			}
		}  ## while(<FHGetStats>) ends
	close FHGetStats;
	}  ## function ends
	###########################
	sub print_help{
	print "\nCalculates the number of spliced reads in given RNA-Seq read mapping bam file
	Usage:
		CountSplicedReads --in-bam <input bam file> --remove-temp <T/F>";
	print "\n";
	}  ## function ends
	###########################
