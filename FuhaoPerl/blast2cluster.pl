#!/usr/bin/env perl
use warnings;
use strict;
use constant USAGE =><<EOH;

usage: $0 perl filter.blast.m6/7

Cluster col1 and col2 of blastout into groups

version: 20150630

EOH
my $blastout=$ARGV[0];
my $groupnumber=0;
my %group=();
my %seqid=();


open (BLASTIN, "< $blastout") || die "Error: Can not open blast result\n";
while (my $line=<BLASTIN>) {
	chomp $line;
	next if ($line=~/^#/);
	my @arr=split(/\t/, $line);
	if (exists $seqid{$arr[0]}) {
		$group{$seqid{$arr[0]}}{$arr[0]}++;
		$group{$seqid{$arr[0]}}{$arr[1]}++;
	}
	elsif (exists $seqid{$arr[1]}) {
		$group{$seqid{$arr[1]}}{$arr[0]}++;
		$group{$seqid{$arr[1]}}{$arr[1]}++;
	}
	else {
		$groupnumber++;
		$seqid{$arr[0]}=$groupnumber;
		$seqid{$arr[1]}=$groupnumber;
		$group{$groupnumber}{$arr[0]}++;
		$group{$groupnumber}{$arr[1]}++;
	}
}
close BLASTIN;
my $count=scalar(keys %group);
my $total=0;
my $total_wheat=0;
my $total_tauschii=0;
foreach my $indgroup (sort {$a <=> $b} keys %group) {
	$total+=scalar(keys %{$group{$indgroup}});
	foreach (keys %{$group{$indgroup}}) {
		$total_wheat++ if ($_=~/^wheat/);
		$total_tauschii++ if ($_=~/^tauschii/);
	}
	my $cluster=join ("\t", keys %{$group{$indgroup}});
	print $indgroup."\t".$cluster."\n";
}
print STDERR "Total: $count\n";
print STDERR "Sum: $total\n";
print STDERR "Average: ". $total/$count ."\n";
print STDERR "Wheat: $total_wheat\n";
print STDERR "tauschii: $total_tauschii\n";
exit 0
