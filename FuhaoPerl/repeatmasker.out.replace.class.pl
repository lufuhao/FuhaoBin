#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 repeatmasker.out convert.class repeatmasker.tab.out

v20181029


### convert.class
### two column: repeatmasker.out[col10]	NewClass
seq000002509	LINE
seq000002515	LINE
seq000000002	LTR
seq000000004	LTR/Gypsy
seq000000012	LTR/Gypsy
seq000000018	LTR/Gypsy
seq000000034	LTR/Copia
seq000000038	LTR/Gypsy

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


my $repeatmasker_out_in=$ARGV[0];
my $convertfile=$ARGV[1];
my $repeatmasker_tab_out=$ARGV[2];
my $numline=0;
my %idhash=();

open (IDS, "<", $convertfile) || die "Error: can not open convert.class\n";
while (my $line=<IDS>) {
	chomp $line;
	$numline++;
	my @arr=split(/\t/, $line);
	unless (scalar(@arr)==2) {
		die "Error: invalid line: $line\n";
	}
	if (exists $idhash{$arr[0]}) {
		if ($idhash{$arr[0]} eq $arr[1]) {
			print STDERR "Info: repeat IDs: $arr[0]\n";
		}
		else {
			die "Warnings: repeat IDs with different repeat class: $arr[0]\n";
		}
	}
	$idhash{$arr[0]} = $arr[1];
}
close IDS;
print "Info: Read ", $numline, " lines\n";
print "Info: Valid ", scalar(keys %idhash), " IDs\n\n";


$numline=0;
my $numvalid=0;
my $number_replace=0;
open (REPEATMASKERIN, "<", $repeatmasker_out_in) || die "Error: invalid repeatmasker.out\n";
open (REPEATMASKEROUT, ">", $repeatmasker_tab_out) || die "Error: invalid repeatmasker.tab.out\n";
while (my $line=<REPEATMASKERIN>) {
	chomp $line;
	$numline++;
	$line=~s/^\s+//;
	$line=~s/\s+/\t/g;
	
	my @arr=split(/\t/, $line);
	unless ($line=~/^\d+\t/) {
		print STDERR "Info: Ignored Line($numline): $line\n";
		next;
	}
	unless (scalar(@arr)>=15 and scalar(@arr)<=16) {
		print STDERR "Warnings: Ignored Line($numline): $line\n";
		next;
	}
	$numvalid++;
	if (exists $idhash{$arr[9]}) {
		$number_replace++;
		$arr[10]=$idhash{$arr[9]};
	}
	print REPEATMASKEROUT join("\t", @arr), "\n";
}
print "Info: Read ", $numline, " lines\n";
print "Info: Valid ", $numvalid, " lines\n";
print "Info: Relaced ", $number_replace, " lines\n";

close REPEATMASKERIN;
close REPEATMASKEROUT;
