#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 R1 R2

Version: 20150508

Descriptions:
  Quick find the difference LINE number between two lists.

EOH

die USAGE unless (scalar(@ARGV) == 2);

our $in1=$ARGV[0];
our $in2=$ARGV[1];
our $line_num=1;
open (IN1, "<$in1") || die "Error: open in1\n";
open (IN2, "<$in2") || die "Error: open in2\n";
while (my $line1=<IN1>) {
	chomp $line1;
	my $line2=<IN2>;
	chomp $line2;
	if ($line1 ne $line2) {
		print $line1."\tvs\t".$line2."at $line_num\n" ;
		die "##### Not Same#####\n";
	}
	else {
		$line_num++;
	}
}
close IN1;
close IN2;
print "##### All the same: $line_num lines #####\n";
