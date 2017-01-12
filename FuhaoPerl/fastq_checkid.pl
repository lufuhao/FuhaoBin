#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

  usage: $0 R1 R2 PATTERN

  Version: 20160104
  
  Pattern:
    Pattern1: '\\@(\\S+)\\/[12]\\s*\\S*'
    Pattern2: '\\@(\\S+)\\s*\\S*'

  Quick find the difference LINE number between two Fastq or fq.gzip.
  exit code 0 if R1 and R2 are paired, otherwise exit 1

EOH

unless (scalar(@ARGV) == 3) {
	print USAGE; 
	exit 1;
}

our $in1=$ARGV[0];
our $in2=$ARGV[1];
our $pattern=$ARGV[2];
our $line_num=0;
if ($in1=~/\.fastq$|\.fq$/) {
	print "#Read1 in fastq format\n";
	open (IN1, "<$in1") || die "Error: open in1\n";
}
elsif ($in1=~/\.gz$|\.gunzip$/) {
	print "#Read1 in gz format\n";
	open (IN1, "zcat $in1 |") || die "Error: open in1\n";
}
else {
	die "Error: Can not guess input1 format\n";
}
if ($in2=~/\.fastq$|\.fq$|\.gzip$/) {
	print "#Read2 in fastq format\n";
	open (IN2, "<$in2") || die "Error: open in2\n";
}
elsif ($in2=~/\.gz$|\.gunzip$/) {
	print "#Read2 in gz format\n";
	open (IN2, "zcat $in2 |") || die "Error: open in2\n";
}
else {
	die "Error: Can not guess input2 format\n";
}
while (my $line1=<IN1>) {
	$line_num++;
	my $id1='';
	if ($line1=~m/$pattern/) {
		$id1=$1;
		die "Error: invalid pattern\n" if ($id1 eq '');
		<IN1> && <IN1> && <IN1>;
	}
	else {
		die "fastq id1 error: $line1 at line: $line_num\n";
  	}
  	my $line2=<IN2>;
  	my $id2='';
  	if ($line2=~m/$pattern/) {
		$id2=$1;
		die "Error: invalid pattern2\n" if ($id2 eq '');
		<IN2> && <IN2> && <IN2>;
	}
	else {
		die "fastq id1 error: $line1 at line: $line_num\n";
  	}
	if ($id1 ne $id2) {
		print $line1."\tvs\t".$line2."at $line_num\n" ;
		die "##### Not Same#####\n";
	}
}
close IN1;
close IN2;
print "##### All the same: $line_num lines #####\n";
exit 0;
