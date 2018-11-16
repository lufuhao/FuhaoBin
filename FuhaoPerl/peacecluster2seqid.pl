#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 fasta_index cluster

samtools faidx xx.fa
$0 xx.fa.fai cluster

EOH
die USAGE if (scalar(@ARGV)!=2);
my $fastaindex=$ARGV[0];
my $clusterfile=$ARGV[1];
die "Invalid fasta index file\n" unless ( -s $fastaindex);
die "invalid pease cluster file\n" unless ( -s $clusterfile);

open (INDEX, "<", $fastaindex) || die "Error: can not open fastaindex file\n";
open (CLUSTER, "<", $clusterfile) || die "Error: can not open cluster file\n";

my $count=0;
my %seqindex=();
while (my $line1=<INDEX>) {
	chomp $line1;
	next if ($line1=~/^#/ or $line1 eq '');
	my @arr1=split(/\s+/, $line1);
	die "Error: Inproper fasta index (Num.Col !=5)\n" unless (scalar(@arr1) ==5);
	$seqindex{$count}=$arr1[0];
#	print STDERR $count."\t".$arr1[0]."\n";### For test ###
	$count++;
}
close INDEX;
print STDERR "Info: Load total of $count sequences\n";

my $estcount=0;
my $num_clusters=0;
my @arr2=();
while (my $line2=<CLUSTER>) {
	chomp $line2;
	if ($line2=~/^Cluster/ or $line2 eq '') {
		if (scalar(@arr2)>0) {
			print join("\t", @arr2)."\n";
			$num_clusters++;
			@arr2=();
		}
		next;
	}
	if ($line2=~/^parentIdx.*estIdx.(\d+).*metric.*$/) {
		print STDERR "EST id: $1\n";### For test ###
		if (exists $seqindex{$1}) {
			push (@arr2, $seqindex{$1});
			$estcount++;
		}
		else {
			die "Error: can not find index $1\n";
		}
	}
}
close CLUSTER;
print STDERR "Print total of $estcount sequence ids\n";
print STDERR "Print total if $num_clusters clusters\n";

0;

