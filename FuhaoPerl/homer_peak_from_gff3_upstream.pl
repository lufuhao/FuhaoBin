#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::FastaKit qw/ReadFastaLength/;
use List::Util qw/min max/;
use constant USAGE =><<EOH;

usage: $0 in.gff3 in.fa out.homer.peaks

v20181101

Prepare upstream 2000 bp peak file for HOMER motif analysis

EOH
die USAGE if (scalar(@ARGV) !=4 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


my $ingff3=$ARGV[0];
my $infasta=$ARGV[1];
my $upstream_length=$ARGV[2];
my $outpeak=$ARGV[3];


my ($test, $length)=ReadFastaLength($infasta);
unless ($test) {
	die "Error: read fasta length\n";
}

open (INGFF3FILE, "<", $ingff3) || die "Error: can not open GFF3 input\n";
open (OUTPEAKFILE, ">", $outpeak) || die "Error: can not write peak output\n";
my $numline=0;
my $geneline=0;
my $outlines=0;
while (my $line=<INGFF3FILE>) {
	chomp $line;
	$numline++;
	next if ($line=~/^#/);
	my @arr=split(/\t/, $line);
	next unless ($arr[2] eq 'gene');
	$geneline++;
	$arr[8]=~s/^.*ID=//;$arr[8]=~s/;.*$//;
	if ($arr[6] eq "+") {
		my $start=max (1, $arr[3]-$upstream_length);
		if ($start<($arr[3]-1)) {
			print OUTPEAKFILE $arr[8],"\t",$arr[0],"\t",$start,"\t", $arr[3]-1,"\t+\n";
			$outlines++;
		}
		else {
			print STDERR "Warnings: no upstream region: $line\n";
		}
	}
	elsif ($arr[6] eq "-") {
		die "Error: invalid seq length: $arr[0]\n" unless (exists ${$length}{$arr[0]} and ${$length}{$arr[0]}=~/^\d+$/);
		my $end=min ($arr[4]+$upstream_length, ${$length}{$arr[0]});
		if (($arr[4]+1)<$end) {
			print OUTPEAKFILE $arr[8],"\t",$arr[0],"\t",$arr[4]+1,"\t", $end,"\t-\n";
			$outlines++;
		}
		else {
			print STDERR "Warnings: no upstream region: $line\n";
		}
	}
	else {
		die "invalid line: $line\n";
	}
	
}
print "Info: read GFF lines: $numline\n";
print "Info: read Gene lines: $geneline\n";
print "Info: out peak lines: $outlines\n";

close INGFF3FILE;
close OUTPEAKFILE;
