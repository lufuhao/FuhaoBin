#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 in.blast in.pairs out.ids out.ids2

DESC
    This script is used to retrieve those scaffolding evidence for a pair of seq
    by blast to a evidence database
    
    For example:
        seq01 and seq02 was expected to be scaffolded
        Blast to a database
        this script is going to get the shared subject for both
        and this can simplify and synteny plot making

in.blast:    Tabular BLAST result
in.pairs:    Scaffolding links, seqname in pairs
out.ids:     Pairs and potential links
out.ids2:    All links seqname

v20170213

EOH
die USAGE if (scalar(@ARGV) !=4 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



my $blastinput=$ARGV[0];
my $pairinput=$ARGV[1];
my $output1=$ARGV[2];
my $output2=$ARGV[3];


my $linenum=0;
my %idallhash=();
my %pairhash=();
my %finalhash=();

open (PAIRS, "< $pairinput") || die "Error: can not open Pairs\n";
while (my $line=<PAIRS>) {
	
	chomp $line;
	next if ($line=~/^#/);
	$linenum++;
	my @arr=split(/\t/, $line);
	unless (scalar(@arr)==2) {
		die "Error: invalid line($linenum): $line\n";
	}
	unless (exists $pairhash{$arr[0]}{$arr[1]} or exists $pairhash{$arr[1]}{$arr[0]}) {
		$pairhash{$arr[0]}{$arr[1]}++;
	}
	else {
		die "Error: repeated pairs at line($linenum): $line\n";
	}
}
close PAIRS;
print "SUM: read total $linenum pairs\n";

open (BLASTOUT, "< $blastinput") || die "Error: can not open blast input\n";
$linenum=0;
while (my $line=<BLASTOUT>) {
	$linenum++;
	chomp $line;
	my @arr=split(/\t/, $line);
	$idallhash{$arr[0]}{$arr[1]}++;
}
close BLASTOUT;
print "SUM: read total $linenum blast\n";



open (OUTPUT1, " > $output1") || die "Error: can not write $output1\n";

foreach my $first (sort keys %pairhash) {
	foreach my $second (sort keys %{$pairhash{$first}}) {
		my %sharedids=();
		my @arr=();
		unless (exists $idallhash{$first} and exists $idallhash{$second}) {
			next;
		}
		foreach (keys %{$idallhash{$first}}) {
			$sharedids{$_}++;
		}
		foreach (keys %{$idallhash{$second}}) {
			$sharedids{$_}++;
		}
		foreach my $idname (sort keys %sharedids) {
			if (exists $idallhash{$first}{$idname} and exists $idallhash{$second}{$idname}) {
				$finalhash{$first}{$idname}++;
				$finalhash{$second}{$idname}++;
				push (@arr, $idname);
			}
			elsif ($sharedids{$idname}==1) {
				next;
			}
			else {
				die "Error: invalid FIRST $first SECOND $second";
			}
		}
		if (scalar(@arr)>0) {
			print OUTPUT1 "$first\t$second\t", join ("\t", @arr), "\n";
		}
	}
}
close OUTPUT1;

open (OUTPUT2, " > $output2") || die "Error: can not write $output2\n";
foreach (sort keys %finalhash) {
	my @arr=keys %{$finalhash{$_}};
	print OUTPUT2 $_, "\t", join ("\t", @arr), "\n";
}
close OUTPUT2;



