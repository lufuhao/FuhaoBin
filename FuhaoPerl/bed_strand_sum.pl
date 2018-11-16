#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage:  out.sum bedfiles

bed_strand_parser.pl ../final.flanking.sort.merge.bed evm.protein.out evm.protein.sum
bed_strand_parser.pl ../final.flanking.sort.merge.bed MIPS.protein.out protein.sum
bed_strand_parser.pl ../final.flanking.sort.merge.bed MIPS.out MIPS.sum

$0 out.bed evm.protein.sum protein.sum MIPS.sum
v20150710

EOH
die USAGE if (scalar(@ARGV) <2 or [0] eq '-h' or [0] eq '--help');

my $output=shift @ARGV;

my $filenum=0;
my %strandsum=();
foreach my $file (@ARGV) {
	$filenum++;
	close BEDIN if (defined fileno(BEDIN));
	open (BEDIN, "< $file") || die "Error: can not open BED file : $file\n";
	while (my $line=<BEDIN>) {
		chomp $line;
		my @arr=split(/\t/, $line);
		$strandsum{"$arr[0]-$arr[1]-$arr[2]"}{$filenum}=$arr[5];
	}
	close BEDIN,
}
print "Total files: $filenum\n";

open (SUMOUT, "> $output") || die "Error: can not write BED sum\n";
#print SUMOUT "#Ref\tStart\tEnd\t". join ("\t", @ARGV) ."\n";
my $stranded=0;
my $totalindex=0;
foreach my $index (sort keys %strandsum) {
	$totalindex++;
	my %indstrand=();
	my ($ref, $start, $end, $strands);
	my @arr=();
	if ($index=~/^(\S+)-(\d+)-(\d+)$/) {
		$ref=$1; $start=$2; $end=$3;
	}
	else {
		print STDERR "Error: unknown bed line in bedhash: $index\n";
		next;
	}
#	print "Test: ", $ref, "\t", $start, "\t", $end, "\n"; ### For test ###
	for (my $i=1; $i<=$filenum; $i++) {
		my $j=$strandsum{$index}{$i};
		push (@arr, $j);
		$indstrand{$j}++ if ($j ne '.');
	}
	my @arr2=keys %indstrand;
	if (scalar(@arr2)==1) {
		$strands=shift @arr2;
		print SUMOUT $ref, "\t", $start, "\t", $end, "\t", join("\t", @arr), "\t", $strands, "\n";
		$stranded++;
	}
	else {
		print SUMOUT $ref, "\t", $start, "\t", $end, "\t", join("\t", @arr), "\t?", "\n";
	}
}
close SUMOUT;
print "\n### SUMMARY ###\n\tTotal hash: ".scalar(keys %strandsum). "\n\tStranded: $stranded\n";
