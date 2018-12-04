#!/usr/bin/perl
use strict;
use warnings;

my $QUALITY_FORMAT="Standard";  ## Could be "Solexa", "Illumina"
my $SAMPLE_NAME = "fqToSam";


## SET PICARD DIRECTORY PATH HERE ##
my $PICARD_BIN="/usr/bin/";





$PICARD_BIN =~ s/(\/)$//;


## executable files ##
my $fqToSam = $PICARD_BIN."/FastqToSam.jar";
my $mergeSam = $PICARD_BIN."/MergeSamFiles.jar";



my $argLen = 4;
	if(scalar(@ARGV) != $argLen){
	print_usage();
	die "\n not enough arguments passed, $argLen required\n";
	}

my $fqIn1 = $ARGV[0];
my $fqIn2 = $ARGV[1];
my $outTag = $ARGV[2];
$QUALITY_FORMAT = $ARGV[3];

	if( ($QUALITY_FORMAT ne "Standard") && ($QUALITY_FORMAT ne "Solexa") && ($QUALITY_FORMAT ne "Illumina") ){
	die "\n quality format not recognized\n";
	}


my($bam1,$bam2,$bam3);
$bam1 = $outTag."_1.fastq.sam";
$bam2 = $outTag."_2.fastq.sam";
$bam3 = $outTag."_merged.sam";

print "Creating first bam file:\n";
`java -Xmx2g -jar $fqToSam FASTQ=$fqIn1 OUTPUT=$bam1 QUALITY_FORMAT=$QUALITY_FORMAT SAMPLE_NAME=$SAMPLE_NAME SORT_ORDER=queryname`;


print "Creating second bam file:\n";
`java -Xmx2g -jar $fqToSam FASTQ=$fqIn2 OUTPUT=$bam2 QUALITY_FORMAT=$QUALITY_FORMAT SAMPLE_NAME=$SAMPLE_NAME SORT_ORDER=queryname`;

print "Merging sam files:\n";
`java -Xmx2g -jar $mergeSam INPUT=$bam1 INPUT=$bam2 OUTPUT=$bam3 SORT_ORDER=queryname`;

print "Removing individual sam files (keeping the merged one):\n";
`rm $bam1`;
`rm $bam2`;



	sub print_usage{
	print "\nUsage:\tFqToSamPicard.pl <fq1> <fq2> <out_tag> <quality_format: Standard/Solexa/Illumina >\n";
	}  ## function ends
	###########################
	
	