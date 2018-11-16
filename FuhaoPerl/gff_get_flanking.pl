#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 input.gene.gff genome.fai output.gff

Add 2000 to 5end and 500 to 3end (require col7 +/-) OR add 2000 to both end
Only operate records with col3==gene

Change \$terminus5 or \$terminus3 for different value

v20150727

EOH
die USAGE if (scalar(@ARGV) !=3 or [0] eq '-h' or [0] eq '--help');



my $gffinput=$ARGV[0];
my $lengthfile=$ARGV[1];
my $gffoutput=$ARGV[2];
my $terminus5=2000;
my $terminus3=500;
die "Error: invalid GFF input\n" unless (defined $gffinput and -s $gffinput);
die "Error: invalid Genome length file\n" unless (defined $lengthfile and -s $lengthfile);
unlink $gffoutput if (defined $gffoutput and $gffoutput=~/^\S+$/);



my $linenum=0;
my %genlen=();
open (LENGTH, "< $lengthfile") || die "Error: can not open Genome Length file\n";
while (my $line1=<LENGTH>) {
	$linenum++;
	chomp $line1;
	my @arr1=split(/\t/, $line1);
	if (exists $genlen{$arr1[0]}) {
		print "Error: existing chrom:length $arr1[0]:$arr1[1]\n";
	}
	else {
		$genlen{$arr1[0]}=$arr1[1];
	}
}
close LENGTH;
print STDERR "### SUMMARY ###\n\tTotalGenomeNum: $linenum\n\tRead GenomeNum: ".scalar(keys %genlen)."\n";


open (GFFOUT, "> $gffoutput") || die "Error: can not write GFF output\n";
open (GFFIN, "< $gffinput") || die "Error: can not open GFF input\n";
$linenum=0;
my $startn=0;
my $endn=0;
my $lineout=0;
while (my $line2=<GFFIN>) {
	$linenum++;
	chomp $line2;
	if ($line2=~/^#/) {
		print GFFOUT $line2."\n";
		$lineout++;
		next;
	}
	my @arr2=split(/\t/, $line2);
	unless ($arr2[3] ne 'gene') {
		print STDERR "Warnings: not GENE at col3\n";
		next;
	}
	my ($start, $end)=($arr2[3], $arr2[4]);
	print "Test: $arr2[6] Start:End $start:$end\n";
	if ($arr2[6] eq '+') {
		$start=$start-$terminus5;
		$end=$end+$terminus3;
	}
	elsif ($arr2[6] eq '-') {
		$start=$start-$terminus3;
		$end=$end+$terminus5;
	}
	else {
		print "Warnings: no +/-\n";
		$start=$start-$terminus5;
		$end=$end+$terminus5;
	}
	if ($start<1) {
#		print STDERR "Warnings: line ($linenum): use start 1\n";
		$startn++;
		$start=1;
	}
	print "Test: $arr2[6] Start:End $start:$end\n";
	if ($end=~/^\d+$/ and $end > $genlen{$arr2[0]}) {
#		print STDERR "Warnings: line ($linenum): use end genome length\n";
		$end=$genlen{$arr2[0]};
		$endn++;
	}
	if ($end<=$start) {
		print STDERR "Warnings: line ($linenum): (end ($end) <= start ($start), Igoring\n";
		next;
	}
	$arr2[3]=$start;
	$arr2[4]=$end;
	print GFFOUT join ("\t", @arr2), "\n";
	$lineout++;
}
close GFFIN;
close GFFOUT;
print STDERR "\n### SUMMARY ###\n\tTotal GFF input lines: $linenum\n\tOutput lines: $lineout\n\tUse start1: $startn\n\tUse end length: $endn\n";
