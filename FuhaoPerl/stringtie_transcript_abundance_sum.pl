#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 output.sum transcript.abundance.gtf [transcript.abundance.gtf transcript.abundance.gtf ...]

v20171012

EOH
die USAGE if (scalar(@ARGV) <3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



my $output=shift @ARGV;
die "Error: Output exists\n" if (-e $output);


my %sumhash=();
my $colnum="TPM";
### 8 = TPM
### 7 = FPKM
my %header=();

for (my $i=0; $i<scalar(@ARGV); $i++) {
	my $idvfile=$ARGV[$i];
	close INPUT if (defined fileno(INPUT));
	open (INPUT, "< $idvfile") || die "Error: can not open GENE ABUNDANCE file: $idvfile\n";
	while (my $line=<INPUT>) {
		next if ($line=~/^#/);
		chomp $line;
		my @arr=split(/\t/, $line);
		next unless ($arr[2] eq 'transcript');
		my $mrnaid=$arr[8];
		my $exprvalue=$arr[8];
		$mrnaid=~s/^.*transcript_id "//; $mrnaid=~s/".*$//;
		$exprvalue=~s/^.*$colnum "//; $exprvalue=~s/".*$//;
		if (exists $sumhash{$mrnaid} and defined $sumhash{$mrnaid}[$i]) {
			die "Error: defined values: FILE $idvfile mRNA $mrnaid\n";
		}
		$sumhash{$mrnaid}[$i]=$exprvalue;
	}
	close INPUT;
}

open (OUTPUT, " > $output") || die "Error: can not write output\n";

print OUTPUT "###GENEID\t", join ("\t", @ARGV), "\n";
foreach my $geneid (sort keys %sumhash) {
	for (my $i=0; $i<scalar(@ARGV); $i++) {
		unless (defined $sumhash{$geneid}[$i]) {
			$sumhash{$geneid}[$i]='NaN';
		}
	}
	print OUTPUT $geneid, "\t", join("\t", @{$sumhash{$geneid}}), "\n";
}
close OUTPUT;
