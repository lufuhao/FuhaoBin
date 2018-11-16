#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 in.gff3 out

Prepare Exon hint for Augustus

version: 20150709

EOH
die USAGE if (scalar(@ARGV) !=2 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

my $gffin=$ARGV[0];
my $gffout=$ARGV[1];
my $src='E';
my $pri=4;
my $mini_count=2;
open (GFFIN, "< $gffin") || die "Error: can not open input GFF";
open (GFFOUT, "> $gffout") || die "Error: can not write outputGFF";
my %count=();
while (my $line=<GFFIN>) {
	chomp $line;
	next if ($line =~/^#/ or $line =~ /^\s+$/);
	my @arr=split(/\t/, $line);
	$count{$arr[0]}{'start'}{$arr[3]}++;
	$count{$arr[0]}{'end'}{$arr[4]}++;
}
close GFFIN;
my %printed=();
open (GFFIN, "< $gffin") || die "Error: can not open input GFF";
while (my $line2=<GFFIN>) {
	chomp $line2;
	next if ($line2 =~/^#/ or $line2 =~ /^\s+$/);
	my @arr=split(/\t/, $line2);
#	print "Test: $arr[0]\t$count{$arr[0]}{'start'}{$arr[3]}\t$count{$arr[0]}{'end'}{$arr[4]}\n";
	if ($count{$arr[0]}{'start'}{$arr[3]} >=$mini_count or $count{$arr[0]}{'end'}{$arr[4]} >=$mini_count){
		if (exists $printed{$arr[0]} and exists $printed{$arr[0]}{"$arr[3]-$arr[4]"}) {
			next;
		}
		else {
			my $max=($count{$arr[0]}{'start'}{$arr[3]}>=$count{$arr[0]}{'end'}{$arr[4]}) ? $count{$arr[0]}{'start'}{$arr[3]} : $count{$arr[0]}{'end'}{$arr[4]};
			$arr[2]='exonpart';
			$arr[8]="src=$src;mult=$max;pri=$pri";
			print GFFOUT join("\t", @arr), "\n";
		}
		$printed{$arr[0]}{"$arr[3]-$arr[4]"}++;
	}
}

close GFFIN;
close GFFOUT;
