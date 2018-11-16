#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 input.SortedByReadname.bam Output.stats

Descriptions
	Extract reads pairs information: Rname + loci
	
Requirements
	samtools

v20160627

EOH
die USAGE if (scalar(@ARGV) !=2 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



my $path_samtools='samtools';
my %referlength=();



open (BAMIN, "$path_samtools view -h $ARGV[0] |") || die "Error: can not open BAM file: $ARGV[0]\n";
open (OUTSTATS, "> $ARGV[1]") || die "Error: can not write to $ARGV[1]\n";
while (my $line=<BAMIN>) {
	chomp $line;
	if ($line=~/^\@/) {
		if ($line=~/^\@SQ\s+SN:(\S+)\s+LN:(\d+)$/) {
			if (exists $referlength{$1}) {
				die "Error: reference duplicates in BAM header: $1\n";
			}
			else {
				$referlength{$1}=$2;
			}
		}
		next;
	}
	my @arr=split(/\t/, $line);
	if (scalar(@arr)<11) {
		print STDERR "Warnings: unknown line: $line\n";
		next;
	}
	if ($arr[1] & 0x0004) {#unmapped
		next;
	}
	my $strand='NaN';
	if ($arr[1] & 0x0010) {
		$strand='-';
	}
	else {
		$strand='+';
	}
	my $matenum=0;
	if ($arr[1] & 0x0040) {
		$matenum=1;
	}
	elsif ($arr[1] & 0x0080) {
		$matenum=2;
	}
	print OUTSTATS $arr[0], "\t", $matenum, "\t", $arr[2], "\t";
	if (exists $referlength{$arr[2]}) {
		print OUTSTATS $referlength{$arr[2]};
	}
	else {
		print OUTSTATS "NaN";
	}
	print OUTSTATS "\t", $strand, "\t", $arr[3], "\n";
}
close BAMIN;
close OUTSTATS;
