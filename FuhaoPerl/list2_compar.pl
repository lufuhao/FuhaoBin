#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 List1 List2

		OR

usage: $0 List1 List2 Share.list Unique1.list Unique2.list

v20161115

Description:
	Compare two list and out share and specific list

EOH
die USAGE if ([0] eq '-h' or [0] eq '--help');
die USAGE unless (scalar(@ARGV) ==2 or scalar(@ARGV) ==5);

my $list1=$ARGV[0];
my $list2=$ARGV[1];
my $testout_share=0;
my $testout_unique1=0;
my $testout_unique2=0;
my $share;
if (defined $ARGV[2]) {
	$share=$ARGV[2];
	$testout_share=1;
	unlink $share if (defined $share and -e $share);
}
my $unique1;
if (defined $ARGV[3]) {
	$unique1=$ARGV[3];
	$testout_unique1=1;
	unlink $unique1 if (defined $unique1 and -e $unique1);
}
my $unique2;
if (defined $ARGV[4]) {
	$unique2=$ARGV[4];
	$testout_unique2=1;
	unlink $unique2 if (defined $unique2 and -e $unique2);
}

my %seqidlist=();
my $linenum=0;
my $validlines=0;



open (LIST1, "< $list1") || die "Error: can not open List1\n";
while (my $line1=<LIST1>) {
	$linenum++;
	chomp $line1;
	if ($line1=~/^\S+$/) {
		$seqidlist{$line1}{1}++;
		$validlines++;
	}
	else {
		print STDERR "Warnings: invalid line at line $linenum of list1 : $line1 \n";
	}
}
close LIST1;
print "\n\n\n### Summary LIST1 ###\n\tTotallines: $linenum\n\tValid lines: $validlines\n";


$linenum=0;
$validlines=0;
open (LIST2, "< $list2") || die "Error: can not open List2\n";
while (my $line2=<LIST2>) {
	$linenum++;
	chomp $line2;
	if ($line2=~/^\S+$/) {
		$seqidlist{$line2}{2}++;
		$validlines++;
	}
	else {
		print STDERR "Warnings: invalid line at line $linenum of list2 : $line2 \n";
	}
}
close LIST2;
print "\n\n\n### Summary LIST2 ###\n\tTotallines: $linenum\n\tValid lines: $validlines\n";


my $sharenum=0;
my $unique1num=0;
my $unique2num=0;


if ($testout_share) {
	open (SHARE, "> $share") || die "Error: can not write share\n";
}
if ($testout_unique1) {
	open (UNIQUE1, "> $unique1") || die "Error: can not write share\n";
}
if ($testout_unique2) {
	open (UNIQUE2, "> $unique2") || die "Error: can not write share\n";
}

foreach (sort keys %seqidlist) {
	my $exists1=1 if (exists $seqidlist{$_}{1});
	my $exists2=1 if (exists $seqidlist{$_}{2});
	
	if ($exists1 and $exists2) {
		$sharenum++;
#		print "S:\t", $_, "\n";
		print SHARE $_."\n" if ($testout_share);
	}
	elsif ($exists1) {
		$unique1num++;
#		print "1\t", $_, "\n";
		print UNIQUE1 $_."\n" if ($testout_unique1);
	}
	elsif ($exists2) {
		$unique2num++;
#		print "2:\t", $_, "\n";
		print UNIQUE2 $_."\n" if ($testout_unique2);
	}
}
close SHARE if ($testout_share);
close UNIQUE1 if ($testout_unique1);
close UNIQUE2 if ($testout_unique2);

print "\n\n\n### Summary ###\n\tShared: $sharenum\n\tUnique1: $unique1num\n\tUnique2: $unique2num\n";
