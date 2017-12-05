#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::FastqKit qw/ShuffleFastq/;
use constant USAGE =><<EOH;

usage: $0 inR1.fastq[.gz] inR2.fastq[.gz] out_shuffled.fastq[.gz]

Descriotion:

    - Shuffle separate R1+R2 fastq into 1 fastq
    - Support gz (Needs Linux: cat, gzip)

v20171205

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

my $filenameA = $ARGV[0];
my $filenameB = $ARGV[1];
my $filenameOut = $ARGV[2];

unless (ShuffleFastq($filenameA, $filenameB, $filenameOut)) {
	
}
