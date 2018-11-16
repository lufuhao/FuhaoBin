#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

$0 bam1.wig1 bam1.wig2

Version: 20150622

EOH

die USAGE if (scalar(@ARGV)!=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

my $wig1=$ARGV[0];
die "Error: invalid WIG1\n" unless (defined $wig1 and -s $wig1);
my $wig2=$ARGV[1];
die "Error: invalid WIG2\n" unless (defined $wig2 and -s $wig2);
my $output=$ARGV[2];
unlink $output if (defined $output and -e $output);

my ($test1, $wig1hash)=&ReadWig($wig1);
if (! $test1) {
	die "Error: ReadWig $wig1\n";
}
my ($test2, $wig2hash)=&ReadWig($wig2);
if (! $test2) {
	die "Error: ReadWig $wig2\n";
}
close WIGOUT if (defined fileno(WIGOUT));
unless (open (WIGOUT, ">$output")) {
	die "Error: write WIG output: $output\n";
}
print WIGOUT "track name=Assigned type=wiggle_0\n";
foreach my $chrom (keys %{$wig1hash}) {
	print "Chrom: $chrom\n";
	die "Error: empty chrom\n" unless ($chrom=~/^\S+$/);
	print WIGOUT "variableStep chrom=$chrom\n";
	foreach my $pos (sort {$a <=>$b} keys %{${$wig1hash}{$chrom}}) {
		if (exists ${${$wig2hash}{$chrom}}{$pos}) {
			if (${${$wig1hash}{$chrom}}{$pos}=~/^\d+$/ and ${${$wig2hash}{$chrom}}{$pos}=~/^\d+$/) {
				my $difference=${${$wig1hash}{$chrom}}{$pos}-${${$wig2hash}{$chrom}}{$pos};
				print WIGOUT "$pos $difference\n"
			}
			else {
				die "Error: none numer detected at chrom:pos $chrom:$pos\n";
			}
		}
		else {
			print "Warnings: wig2 do not have chrom:pos $chrom:$pos, set to 0\n";
			print WIGOUT "$pos ${${$wig1hash}{$chrom}}{$pos}\n"
		}
	}
}
close WIGOUT;






sub ReadWig {
	my $RWwig=shift;
	
	my $RWsubinfo='SUB(ReadWig)';
	my %RWwighash=();
	my $RWline_num=0;
	my $RWchrom='';
	
	return 0 unless (defined $RWwig and -s $RWwig);
	close WIG if (defined fileno(WIG));
	unless (open(WIG, "<$RWwig")) {
		print STDERR "${RWsubinfo}Error: can not open WIG file\n";
		return 0;
	}
	
	while (my $line=<WIG>) {
		chomp $line;
		$RWline_num++;
		next if ($line=~/^#/);
		next if ($line=~/^track/);
		if ($line=~/^variableStep\s+chrom=(\S+)$/) {
			$RWchrom=$1;
			print "${RWsubinfo}Test: $RWchrom\n";### For test ###
			unless (defined $RWchrom and $RWchrom=~/^\S+$/) {
				print STDERR "${RWsubinfo}Error: wrong chrom ID at line $RWline_num of WIG $RWwig\n";
				close WIG;
				return 0;
			}
			next;
		}
		if ($line=~/^(\d+)\s+(\d+)$/) {
			my @RWarr=split(/\s+/, $line);
			${$RWwighash{$RWchrom}}{$1}=$2;
		}
	}
	close WIG;
	return (1, \%RWwighash);
}
