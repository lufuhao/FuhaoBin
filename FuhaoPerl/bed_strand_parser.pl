#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

Usage:
  $0 input.bed 12col.intersectbed out.bed

v20150805

Descriptions:
  bedtools intersect -wa -wb -f 1.00 -a xx.gff -b yy.bed > out.intersect.bed
  $0 yybed out.intersect.bed zz.out


EOH
die USAGE if (scalar(@ARGV) !=3 or [0] eq '-h' or [0] eq '--help');


my $inputbed=$ARGV[0];
die "Error: invalid input bed format\n" unless (defined $inputbed and -s $inputbed);
my $intersectbed=$ARGV[1];
die "Error: invalid intersectbed\n" unless (defined $intersectbed and -s $intersectbed);
my $outbed=$ARGV[2];
unlink $outbed if (-s $outbed);

my $linenum=0;
my %bedhash=();
open (BEDIN, "< $inputbed") || die "Error: can not open input bed file\n";
while (my $line=<BEDIN>) {
	$linenum++;
	chomp $line;
	$line=~s/\s+/-/g;
#	print $line."\n"; ### For test ###
	$bedhash{$line}++;
}
close BEDIN;
print "\n### SUMMARY ###\n\tTotal bedline: $linenum\n\tTotal hash: ". scalar(keys %bedhash)."\n";


my %strand=();
$linenum=0;
open (INTERSECT, "< $intersectbed") || die "Error: \n";
while (my $line=<INTERSECT>) {
	$linenum++;
	chomp $line;
	next if ($line=~/^#/);
	my @arr=split(/\t/, $line);
	my $index=$arr[9].'-'.$arr[10].'-'.$arr[11];
#	print $index."\n"; ### For test ###
	unless ($arr[6]=~/^[+-]{1,1}$/) {
		print STDERR "Error: line $linenum: unknown strand $arr[6]\n";
		next;
	}
	unless (exists $bedhash{$index}) {
		print STDERR "Error: line $linenum: unknown bedline $index\n";
		next;
	}
	$strand{$index}{$arr[6]}++;
}
close INTERSECT;
print "\n### SUMMARY (intersectbed) ###\n\tTotal bedline: $linenum\n\tTotal hash: ". scalar(keys %strand)."\n";

my $numstranded=0;
my $numunknown=0;
my $nummixed=0;

open (BEDOUT, "> $outbed") || die "Error: can not write BED out\n";
foreach my $index (sort keys %bedhash) {
	my ($ref, $start, $end, $strand);
	if ($index=~/^(\S+)-(\d+)-(\d+)$/) {
		$ref=$1; $start=$2; $end=$3;
	}
	else {
		print STDERR "Error: unknown bed line in bedhash: $index\n";
		next;
	}
	if (exists $strand{$index}) {
		my @arr2=keys %{$strand{$index}};
		if (scalar(@arr2)==1) {
			$strand=shift @arr2;
			print BEDOUT $ref, "\t", $start, "\t", $end, "\t.\t.\t", $strand, "\n";
			$numstranded++;
		}
		else {
			print BEDOUT $ref, "\t", $start, "\t", $end, "\t.\t.\t?", "\n";
			$nummixed++;
		}
	}
	else {
		print BEDOUT $ref, "\t", $start, "\t", $end, "\t.\t.\t.", "\n";
		$numunknown++;
	}
}
close BEDOUT;
print "\n### SUMMARY (intersectbed) ###\n\tStranded: $numstranded\n\tUnknown: $numunknown\n\tMixed: $nummixed\n";
