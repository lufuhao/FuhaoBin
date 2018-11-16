#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

Usage: bamreadcount input.bam

[Memory depends on your read number]

EOH
my $bamin=$ARGV[0];
die USAGE if (scalar(@ARGV)!=1 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

die "Error: invalid bam file\n" unless ( -s $bamin);

if ($bamin=~/\.bam$/i) {
	open (BAMIN, "samtools view $bamin |") || die "Error: can not open SAM: $bamin\n";
}
elsif ($bamin=~/\.sam$/i) {
	open (BAMIN, "samtools view -S $bamin |") || die "Error: can not open SAM: $bamin\n";
}
else {
	die "Error: Can not guess format: $bamin\n";
}
my %readcount=();
my $total_alignments=0;
while (my $line=<BAMIN>) {
	$total_alignments++;
	chomp $line;
	my @arr1=split(/\t/, $line);
	next if ($arr1[1] & 0x0004);
	my $readpair=0;
	if ($arr1[1] & 0x0040) {
		$readpair=1;
	}
	elsif ($arr1[1] & 0x0080) {
		$readpair=2;
	}
	if ($readpair==0) {
		die "Error: Unknow reads pair number\n"
	}
	my $annotated=0;
	if ($line=~/zg:Z:N/) {
		$annotated='N';
	}
	else {
		$annotated='Y';
	}
	${${$readcount{$arr1[0]}}{$readpair}}{$annotated}++;
}

my $total_paired_anno=0;
my $total_paired_no=0;
my $total_1_anno=0;
my $total_1_no=0;
my $total_2_anno=0;
my $total_2_no=0;
foreach my $indread (keys %readcount) {
	if (exists ${$readcount{$indread}}{1} and exists ${$readcount{$indread}}{2}) {
		if (exists ${${$readcount{$indread}}{1}}{'Y'} and exists ${${$readcount{$indread}}{2}}{'Y'}) {
			$total_paired_anno++;
		}
		elsif (exists ${${$readcount{$indread}}{1}}{'N'} and exists ${${$readcount{$indread}}{2}}{'N'}) {
			$total_paired_no++;
		}
		else {
			die "Error: deferent annotation of pairs: $indread in $bamin\n";
		}
	}
	elsif (exists ${$readcount{$indread}}{1}) {
		if (exists ${${$readcount{$indread}}{1}}{'Y'}) {
			$total_1_anno++;
		}
		elsif (exists ${${$readcount{$indread}}{1}}{'N'}) {
			$total_1_no++;
		}
		else {
			die "Error: unknown read1 annotation: $indread in $bamin\n";
		}
	}
	elsif (exists ${$readcount{$indread}}{2}) {
		if (exists ${${$readcount{$indread}}{2}}{'Y'}) {
			$total_2_anno++;
		}
		elsif (exists ${${$readcount{$indread}}{2}}{'N'}) {
			$total_2_no++;
		}
		else {
			die "Error: unknown read1 annotation: $indread in $bamin\n";
		}
	}
	else {
		die "Error: unknown read pairs: $indread in $bamin\n";
	}
}

print "Summary: $bamin\n";
print "Total alignment:\t$total_alignments\n";
print "Total pairs annotated:\t$total_paired_anno\n";
print "Total pairs un-annotated: $total_paired_no\n";
print "Total read1 annotated: $total_1_anno\n";
print "Total read1 un-annotated: $total_1_no\n";
print "Total read2 annotated: $total_2_anno\n";
print "Total read2 un-annotated: $total_2_no\n\n\n";
exit 0;
