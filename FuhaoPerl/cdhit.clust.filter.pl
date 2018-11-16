#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 in.cdhit.clst out.filtered.clst

Desc: filter cdhit cluster file with those clustered clusters

v20151109

EOH
die USAGE if (scalar(@ARGV) !=2 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



my $cdhitclstin=$ARGV[0];
chomp $cdhitclstin;
die "Error: invalid input cdhit cluster file\n" unless (defined $cdhitclstin and -s $cdhitclstin);
my $CDhitclstout=$ARGV[1];
die "Error: invalid output cdhit cluster file\n" unless (defined $CDhitclstout);
unlink $CDhitclstout if (-e $CDhitclstout);


my %hash=();
my $clustnum;

open (CLUSTIN, "< $cdhitclstin") || die "Error: can not open clustin file: $cdhitclstin\n";
while (my $line=<CLUSTIN>) {
	chomp $line;
	if ($line=~/^>Cluster\s+(\d+)/) {
		$clustnum=$1;
	}
	else {
		$hash{$clustnum}++ if (defined $clustnum);
	}
}
close CLUSTIN;

print "total cluster: ", scalar(keys %hash), "\n";
my $filtered=0;
my $testprint=0;

open (CLUSTIN2, "< $cdhitclstin") || die "Error: can not open clustin file: $cdhitclstin\n";
open (CLUSTOUT, "> $CDhitclstout") || die "Error: can not write clustout file: $CDhitclstout\n";
while (my $line=<CLUSTIN2>) {
	chomp $line;
	if ($line=~/^>Cluster\s+(\d+)/) {
		$clustnum=$1;
		if (exists $hash{$clustnum} and $hash{$clustnum}>1) {
			$testprint=1;
			$filtered++;
			print CLUSTOUT $line, "\n";
		}
		else {
			$testprint=0;
		}
	}
	else {
		print CLUSTOUT $line, "\n" if ($testprint);
	}
}
close CLUSTIN2;
close CLUSTOUT;

print "total cluster filter: $filtered\n";
