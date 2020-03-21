#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 jcvi.bed1 convertfile jcvi.bed2

This script is used to convert JCVI jcvi.formats.gff bed col4 to a new ID

### jcvi.bed1
chr[tab]start[tab]end[tab]ID_bed1[tab]0[tab]strand

### convertfile [2 columns]
ID_bed1[tab]new_ID_bed2

### ### jcvi.bed2
chr[tab]start[tab]end[tab]new_ID_bed2[tab]0[tab]strand

EOH
die USAGE if (scalar(@ARGV)!=3);

my %idhash=();


open (CONVERTID, "<", $ARGV[1]) || die "Error: can not open convertfile\n";

while (my $line=<CONVERTID>) {
	chomp $line;
	my @arr=();
	@arr=split("\t", $line);
	if (exists $idhash{$arr[0]}) {
		die "Error: duplicated ID: $arr[0]\n";
	}
	else {
		$idhash{$arr[0]}=$arr[1];
	}
}
close CONVERTID;


my $linenum_bed1=0;
my $linenum_bed2=0;

open(BED1, "<", $ARGV[0]) || die "Error: can not open bed1\n";
open(BED2, ">", $ARGV[2]) || die "Error: can not open bed2\n";
while (my $line=<BED1>) {
	chomp $line;
	$linenum_bed1++;
	my @arr=();
	@arr=split("\t", $line);
	if (exists $idhash{$arr[3]}) {
		$arr[3]=$idhash{$arr[3]};
		$linenum_bed2++;
	}
	print BED2 join("\t", @arr), "\n";
}
close BED1;
close BED2;
print "Info: Total input lines:     $linenum_bed1\n";
print "Info: Total convert IDs:     ", scalar(keys %idhash), "\n";
print "Info: Total converted lines: $linenum_bed2\n";
