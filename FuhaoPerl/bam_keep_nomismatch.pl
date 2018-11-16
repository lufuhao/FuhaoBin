#!/usr/bin/env perl
use warnings;
use strict;
use constant USAGE =><<EOH;

usage: $0 bam min_mismatch sam/read

STDOUT: sam/reads
Version: 20150409

EOH

die USAGE if (scalar(@ARGV) !=3 or ! -s $ARGV[0]);



###Default
my $bamfile=$ARGV[0];
my $minmatch= (defined $ARGV[1] and $ARGV[1]=~/^\d+$/) ? $ARGV[1] : 0;
print STDERR "Info: minimum mismatch: $minmatch\n";



###Input and output
my $outformat=$ARGV[2];
if ($outformat=~/^sam$/i) {
	print STDERR "Output Format: SAM\n";
	$outformat='sam';
}
elsif ($outformat=~/^read$/i) {
	print STDERR "Output Format: read ID\n";
	$outformat='read';
}
else {
	die "Error: unknown output format\n";
}



###Main
if ($bamfile=~/\.sam$/i) {
	open (INPUT, "samtools view -S $bamfile |") || die "Error: can not open sam $bamfile\n";
}
elsif ($bamfile=~/\.bam$/i) {
	open (INPUT, "samtools view $bamfile |") || die "Error: can not open bam $bamfile\n";
}
else {
	die "Error: can not guess format SAM or BAM for $bamfile\n";
}

while (my $line=<INPUT>) {
	if ($line=~m/^\@/) {
		print $line if ($outformat eq'sam');
		next;
	}
	if ($line=~/MD:Z:\d+\s+/) {
		my @arr=split(/\t/, $line); 
		next if ($arr[5] =~ /D/);
		unless ($arr[1] & 2048) {
			my $count=0;
			while ($arr[5]=~/(\d+)M/g) {
				$count+=$1;
			}
			next if ($count<$minmatch);
		}
		next if ($arr[1] & 256);
		if ($outformat eq 'sam') {
			print $line;
		}
		else {
			print $arr[0]."\n";
		}
	}
}
close INPUT;
exit 0;
