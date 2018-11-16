#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 gff cds.out

v20160613

EOH
die USAGE if (scalar(@ARGV) !=2 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


my %trans=();
my $linenum=0;

open (GFFIN, "< $ARGV[0]") || die "Error: can not open GFF input\n";
while (my $line=<GFFIN>) {
	$linenum++;
	chomp $line;
	my @arr=split (/\t/, $line);
	my $transid=$arr[8];
	if ($line=~/\tCDS\t/) {
		$transid=~s/^.*Parent=//; $transid=~s/;.*$//;
		if (exists $trans{$transid} and exists $trans{$transid}{'cds'} and exists $trans{$transid}{'cds'}{$arr[3]}) {
			print STDERR "Error CDS: line ($linenum): $line\n";
			exit 1;
		}
		$trans{$transid}{'cds'}{$arr[3]}=$arr[4];
		$trans{$transid}{'strand'}{$arr[6]}++;
	}
	if ($line=~/\texon\t/) {
		$transid=~s/^.*Parent=//; $transid=~s/;.*$//;
		if (exists $trans{$transid} and exists $trans{$transid}{'exon'} and exists $trans{$transid}{'exon'}{$arr[3]}) {
			print STDERR "Error exon: line ($linenum): $line\n";
			exit 1;
		}
		$trans{$transid}{'exon'}{$arr[3]}=$arr[4];
		$trans{$transid}{'strand'}{$arr[6]}++;
	}
}
close GFFIN;

open (CDSOUT, "> $ARGV[1]") || die "Error: can not write CDS structure\n";
foreach my $transid (sort keys %trans) {
	my @strands=keys %{$trans{$transid}{'strand'}};
	die "Error: $transid strand @strands\n" if (scalar(@strands)!=1);
	my $strand=$strands[0];
	my @length=();
	my $test=0;
	my @introns=();
	my $lastend;
	foreach my $start (sort {$a<=>$b} keys %{$trans{$transid}{'exon'}}) {
		my $end =$trans{$transid}{'exon'}{$start};
		if ($test==0) {
			$lastend=$end;
			$test++;
		}
		else {
			push (@introns, $start-$lastend-1);
			$lastend=$end;
		}
		my $cdslength=$end-$start+1;
		push (@length, $cdslength);
		
	}
	if ($strand eq '-') {
		@length=reverse @length;
		@introns=reverse @introns;
	}
	print CDSOUT $transid, "\tEXON\t", join ("\t", @length), "\tIntrons\t", join ("\t", @introns);
	if (exists $trans{$transid}{'cds'}) {
		$test=0;
		@length=();
		@introns=();
		foreach my $start (sort {$a<=>$b} keys %{$trans{$transid}{'cds'}}) {
			my $end =$trans{$transid}{'cds'}{$start};
			if ($test==0) {
				$lastend=$end;
				$test++;
			}
			else {
				push (@introns, $start-$lastend-1);
				$lastend=$end;
			}
			my $cdslength=$end-$start+1;
			push (@length, $cdslength);
		
		}
		if ($strand eq '-') {
			@length=reverse @length;
			@introns=reverse @introns;
		}
		print CDSOUT "\tCDS\t", join ("\t", @length), "\tIntrons\t", join ("\t", @introns);
	}
	print CDSOUT "\n";
}
close CDSOUT;
