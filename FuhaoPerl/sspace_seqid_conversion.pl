#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 sspace_formatted.fa name.convertion

convert seqid formatted by SSPACE to it's original

input format in fasta
>contigDDDD|sizeDDDD|seed:previousname
sequence
>...

Output
currentname	previous


v20160927

EOH
die USAGE if (scalar(@ARGV) !=2 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



die "Error: invalid fasta input\n" unless (defined "$ARGV[0]" and -s "$ARGV[0]");
my $inputfasta=$ARGV[0];
die "Error: invalid output\n" unless (defined "$ARGV[1]");
my $outputsum=$ARGV[1];
my $linenum=0;
my $seqnum=0;


open (INPUTFASTA, "< $inputfasta") || die "Error: invalid input fasta\n";
open (OUTPUT, "> $outputsum") || die "Error: can not write output\n";
while (my $line=<INPUTFASTA>) {
	chomp $line;
	$linenum++;
	next unless ($line=~/^>/);
	$seqnum++;
	if ($line=~/^>(contig\d+)\|size\d+\|seed:(\S+)$/) {
		print OUTPUT $1, "\t", $2, "\n";
	}
	else {
		die "Error: invalid fasta header: $line\n";
	}
}
close INPUTFASTA;
close OUTPUT;


print "\n\n\n### Summary: Total lines: $linenum\nTotal seq num: $seqnum\n\n\n";
