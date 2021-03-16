#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage:  masked.fa out.gff

Annotate Ns from repeatmasked fasta file into GFF3
Recognise the Ns coordinates

v20171205

EOH
die USAGE if (scalar(@ARGV) !=2 or [0] eq '-h' or [0] eq '--help');

my $mymaskedfasta=$ARGV[0];
my $outgff3=$ARGV[1];

open (MASKED, "< $mymaskedfasta") || die "Error: can not open fasta file\n";
open (GFF3, "> $outgff3") || die "Error: can not write gff file\n";
my @arr=();
$arr[1]='HomeMade';
$arr[2]='repeat';
$arr[5]='.';
$arr[6]='+';
$arr[7]='.';
my $test_start=0;
my $test_end=1;
my $start;
my $end;
my $seqid;
my $basenum=0;
while (my $line=<MASKED>) {
	chomp $line;
	if ($line=~/^>(\S+)/) {
		if ($test_start==1) {
			$arr[3]=$start;
			if ($test_end==0) {
				$arr[4]=$basenum;
				$arr[8]="Target=$seqid";
				if ($arr[4]>=$arr[3]) {
					print GFF3 join("\t", @arr), "\n";
				}
				else {
					print STDERR "Warnings Ignore as end $arr[4] < start $arr[3]\n";
				}
			}		
			$test_start=0;
		}
		$seqid=$1;
		$arr[0]=$seqid;
		$basenum=0;
	}
	else {
		my @arr1=split(//, $line);
		for (my $i=0; $i<scalar(@arr1); $i++) {
			$basenum++;
			if ($arr1[$i] =~/^[nN]$/) {
				if ($test_start==0) {
					$start=$basenum;
					$test_start=1;
					$test_end=0;
				}
			}
			else {
				if ($test_start==1) {
					$end=$basenum-1;
					$arr[3]=$start;
					$arr[4]=$end;
					$arr[8]="Target=$seqid";
					if ($arr[4]>=$arr[3]) {
						print GFF3 join("\t", @arr), "\n";
					}
					else {
						print STDERR "Warnings Ignore as end $arr[4] < start $arr[3]\n";
					}
					$test_start=0;
					$test_end=1;
				}
			}
		}
	}
}
close MASKED;
close GFF3;
