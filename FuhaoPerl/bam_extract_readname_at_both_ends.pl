#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 input.bam genome.da.fai outreads out.bam

Descriptions:
	Extract read names mapped to each reference ends
	Default: end length = 60KB

Requirements:
	samtools

v20160627

EOH
die USAGE if (scalar(@ARGV) !=4 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


if (defined $ARGV[2] and -s $ARGV[2]) {
	print STDERR "Warnings: output read file $ARGV[2] exists, deleting\n";
	unlink $ARGV[2];
}
if (defined $ARGV[3] and -s $ARGV[3]) {
	print STDERR "Warnings: output BAM file $ARGV[3] exists, deleting\n";
	unlink $ARGV[3];
}

my $maxinsert=60000;
my $path_samtools='samtools';

my %faidx=();
my $linenum=0;
my %readhash=();
open (FASTAINDEX, "< $ARGV[1]") || die "Error: Can not open index file of genome fasta\n";
while (my $line=<FASTAINDEX>) {
	chomp $line;
	$linenum++;
	next if ($line=~/^#/);
	my @arr=split(/\t/, $line);
	unless (exists $faidx{$arr[0]}) {
		$faidx{$arr[0]}=$arr[1];
	}
	else {
		die "Error: duplicated fasta ID: $arr[0]\n";
	}
}
close FASTAINDEX;

print STDERR "Info: read $linenum lines from fasta index: $ARGV[1]\n";
print STDERR "Info: Load ".scalar(keys %faidx)." sequences from fasta index: $ARGV[1]\n\n";

$linenum=0;
my $linecom=0;
my $linealign=0;
open (BAMIN, "$path_samtools view -h $ARGV[0] | ") || die "Error: can not open BAM file: $ARGV[0]\n";
while (my $line=<BAMIN>) {
	chomp $line;
	$linenum++;
	if ($line=~/^\@/) { ### comments lines
		$linecom++;
		next;
	}
	$linealign++;
	my @arr=split(/\t/, $line);
	if (scalar(@arr)<12) {
		print STDERR "Warnings: invalid alignment (line$linenum): $line\n";
		next;
	}
	if ($arr[1] & 0x0004) {
		next;
	}
#	if ($arr[3] >0 and $arr[3] < $maxinsert){
#		$readhash{$arr[0]}++;
#	}
#	elsif ($arr[3] > ($faidx{$arr[2]}-$maxinsert)) {
#		$readhash{$arr[0]}++;
#	}
	if ($faidx{$arr[2]}> 2 * $maxinsert) {
		if ($arr[3]<=$maxinsert) {
			$readhash{$arr[0]}++ if ($arr[1] & 0x0010);
		}
		elsif ($arr[3]>=($faidx{$arr[2]} -$maxinsert)) {
			$readhash{$arr[0]}++ unless ($arr[1] & 0x0010);
		}
	}
	elsif ($faidx{$arr[2]}>= 2 * $maxinsert and $faidx{$arr[2]}< 2 * $maxinsert) {
		if ($arr[3]<=($faidx{$arr[2]} -$maxinsert)) {
			$readhash{$arr[0]}++ if ($arr[1] & 0x0010);
		}
		elsif ($arr[3]>=$maxinsert) {
			$readhash{$arr[0]}++ unless ($arr[1] & 0x0010);
		}
		else {
			$readhash{$arr[0]}++;
		}
	}
	else {
		$readhash{$arr[0]}++;
	}
}
close BAMIN;


print STDERR "Info: read BAM $linenum lines from $ARGV[0]\n";
print STDERR "Info: read BAM header $linecom lines from $ARGV[0]\n";
print STDERR "Info: read BAM alignment $linealign lines from $ARGV[0]\n";
print STDERR "Info: read hash ".scalar(keys %readhash)." from $ARGV[0]\n\n";


%faidx=();
open (READOUT, "> $ARGV[2]") || die "Error: can not write to $ARGV[2]\n";
foreach my $readid (sort keys %readhash) {
	print READOUT $readid, "\n";
}
close READOUT;

$linecom=0;
$linenum=0;
$linealign=0;
open (BAMIN2, "$path_samtools view -h $ARGV[0] | ") || die "Error: can not open BAM file: $ARGV[0]\n";
open (BAMOUT, " | $path_samtools view -h -S -b - > $ARGV[3]") || die "Error: can not write BAM file: $ARGV[3]\n";
while (my $line=<BAMIN2>) {
	chomp $line;
	if ($line=~/^\@/) { ### comments lines
		$linecom++;
		$linenum++;
		print BAMOUT $line, "\n";
		next;
	}
	my @arr=split(/\t/, $line);
	if (scalar(@arr)<12) {
		print STDERR "Warnings: invalid alignment (line$linenum): $line\n";
		next;
	}
	if ($arr[1] & 0x0004) {
		next;
	}
	if (exists $readhash{$arr[0]}){
		print BAMOUT $line, "\n";
		$linealign++;
		$linenum++;
	}
}

close BAMIN2;
close BAMOUT;
print STDERR "Info: write BAM $linenum lines to $ARGV[3]\n";
print STDERR "Info: write BAM header $linecom lines to $ARGV[3]\n";
print STDERR "Info: write BAM alignment $linealign lines to $ARGV[3]\n";
