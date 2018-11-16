#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage:  hapcompass.fragment haptree.fragment

v20150710

EOH
die USAGE if (scalar(@ARGV) !=2 or [0] eq '-h' or [0] eq '--help');
my $hapcompassfrag=$ARGV[0];
my $haptreefrag=$ARGV[1];

my $linenum=0;
open (FRAGMENT, "< $hapcompassfrag") || die "Error: can not open hapcompass fragments\n";
open (OUTPUT, "> $haptreefrag") || die "Error: can not write output\n";
while (my $line=<FRAGMENT> ) {
	$linenum++;
	chomp $line;
	my @arr1=split (/\t/, $line);
	unless (scalar(@arr1)>=3) {
		print STDERR "Warnings: no fragment at line($linenum): $line\n";
		next;
	}
	my @allele=();
	for (my $i=2; $i<scalar(@arr1); $i++) {
		my @arr2=split(/\s+/, $arr1[$i]);
		shift @arr2;
		unless (scalar(@arr2)%3 == 0) {
			print STDERR "Warnings: fragment != 3* : $arr1[$i]\n";
			next;
		}
		while (@arr2) {
			my $lineref=shift @arr2;
			my $allelecode=shift @arr2;
			shift @arr2;
			push (@allele, "$lineref: $allelecode");
		}
	}
	print OUTPUT '{', join (', ', @allele), '}', "\n";
	
}
close FRAGMENT;
close OUTPUT;

