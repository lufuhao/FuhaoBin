#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: LUFUHAO20150424

Requirements:
	Programs:
	Modiles: Getopt::Long

Descriptions:
	Given read id infomation, annotate each reads in a Bam
	Output BAM.annot.bam and BAM.log

Options:
	--help|-h
		Print this help/usage;
	--input|-i	<File>
		[Msg] 3 columns of input file (read_ID\tGenome\tChromosome)
		Support gzip file (*.gz)
	--bam|-b	<Bam/Bams>
		[Msg] Bam file to be annotated (comma delimited)
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
our ($help, $verbose, $debug, $ver);
our ($input, $bam);

GetOptions(
	"help|h!" => \$help,
	"input|i:s" => \$input,
	"bam|b:s" => \$bam,
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$ver) or die USAGE;
($help or $ver) and die USAGE;



### Defaults ########################################################




### input and output ################################################
die "Error: annotation input: $input\n" unless (-s $input);
our @bamfiles=split(/,/, $bam);
map {die "Error: invalid bam file: $_\n" unless (-s $_)} @bamfiles;



### Main ############################################################
if ($input=~/\.gz/) {
	open (ANNO, "zcat $input | ") || die "Error: can not open $input\n";
}
else {
	open (ANNO, "cat $input | ") || die "Error: can not open $input\n";
}

my $num_reads=0;
our %read_anno=();
our %readcount_per_chrom=();
while (my $line=<ANNO>) {
	chomp $line;
	my @arr=split(/\t/, $line);
	if (scalar(@arr) != 3) {
		print STDERR "---Unknown annotation (!=3): $line\n";
		next;
	}
	@{$read_anno{$arr[0]}}=($arr[1], $arr[2]);
	my @arr2=&SplitChrom($arr[2]);
	map {$readcount_per_chrom{$_}++} @arr2;
}
close ANNO;


foreach my $bamfile (@bamfiles) {
	my $num_annotalign=0;
	my $total_annotalign=0;
	my $output=&RetrvNoExt($bamfile);
	print "---Input: $bamfile\n";
	print "---Output Bam: $output.annot.bam\n";
	open (BAM, "samtools view -h $bamfile |") || die "Error: can not open BAM: $bam\n";
	open (OUT, " | samtools view -bS - > $output.annot.bam") || die "Error: can not write output $output.annot.bam\n";
	while (my $line2=<BAM>) {
		if ($line2=~/^\@/) {
			print OUT $line2;
			next;
		}
		chomp $line2;
		$total_annotalign++;
		my @arr2=split(/\t/, $line2);
		my $tag_zg='N';
		my $tag_zc='N';
		if (exists $read_anno{$arr2[0]}) {
			$num_annotalign++;
			$tag_zg=${$read_anno{$arr2[0]}}[0];
			$tag_zc=${$read_anno{$arr2[0]}}[1];
		}
		if ($line2=~/\s+zg:Z/) {
			$line2=~s/zg:Z\S+/zg:Z:$tag_zg/g;
		}
		else {
			$line2.="\tzg:Z:".$tag_zg;
		}
		if ($line2=~/\s+zc:Z/) {
			$line2=~s/zc:Z\S+/zc:Z:$tag_zc/g;
		}
		else {
			$line2.="\tzc:Z:".$tag_zc;
		}
#		$line2.="\tzg:Z:".$tag_zg."\tzc:Z:".$tag_zc;
		print OUT $line2."\n";
	}
	close OUT;
	close BAM;


	print "---Output log: $output.log\n";
	open (LOG, ">$output.log") || die "Error: Can not write $output.log\n";
	print LOG "Total alignments: $total_annotalign\nAnnotated alignments: $num_annotalign\n";
	foreach (sort keys %readcount_per_chrom) {
		print LOG $_."\t".$readcount_per_chrom{$_}."\n";
	}
	close LOG;
	print "Total alignments: $total_annotalign\nAnnotated alignments: $num_annotalign\n";
}



#####################################################################
###                         sub functions                         ###
#####################################################################
### SplitChrom
###&SplitChrom($STR)
###Global:
###Dependency:
###Note:
sub SplitChrom {
	my $SCstr=shift;
	my @SCarr=();
	while ($SCstr=~/(\d{1}[A-Z]{1,2})/g) {
		push @SCarr, $1;
	}
	return @SCarr;
}



###Retrieve filebasename without extension
###& RetrvNoExt(file)
###Global:

sub RetrvNoExt {
	my $RNE_ori=shift @_;
	chomp $RNE_ori;
	my $RNE_new='';
	my $RNE_base='';
	($RNE_base=$RNE_ori)=~ s/.*\///s;
	($RNE_new=$RNE_base)=~s/^(\S+)\.\w+$/$1/;
	return $RNE_new;
}
