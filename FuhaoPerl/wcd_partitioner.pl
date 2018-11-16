#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage:  source.in IDfile out

v20150727

EOH
die USAGE if (scalar(@ARGV) !=3 or [0] eq '-h' or [0] eq '--help');
my $sourcefile=$ARGV[0];
my $idlistfile=$ARGV[1];
my $fileoutput=$ARGV[2];
die "Error: invalid source.in file\n" unless (defined $sourcefile and -s $sourcefile);
die "Error: invalid ID list file\n" unless (defined $idlistfile and -s $idlistfile);
unlink $fileoutput if (-e $fileoutput);
my $linenum=0;
my %idlist=();



open (IDLIST, "< $idlistfile") || die "Error: can not open ID list file\n";
while (my $line1=<IDLIST>) {
	$linenum++;
	chomp $line1;
	my @arr1=split(/\s+/, $line1);
	foreach (@arr1) {
		if (exists $idlist{$_}) {
			print STDERR "Error: repeated ID: $_\n";
			next;
		}
		else {
			$idlist{$_}++;
		}
	}
}
close IDLIST;
print "\n### SUMMARY ###\n\tTotal_line: $linenum\n\tTotalID: ".scalar(keys %idlist)."\n";


$linenum=0;
my $outputlines=0;
open (SOURCEIN, "<$sourcefile") || die "Error: can not open source.in file\n";
open (OUTPUT, "> $fileoutput") || die "Error: can not write output file\n";
while (my $line2=<SOURCEIN>) {
	$linenum++;
	chomp $line2;
	my @arr2=split(/\s+/, $line2);
	if (exists $idlist{$arr2[0]}) {
		$outputlines++;
		print OUTPUT $line2."\n";
	}
}
close SOURCEIN;
close OUTPUT;
print "\n### SUMMARY ###\n\tTital_line: $linenum\n\tOutput_line: $outputlines\n";
