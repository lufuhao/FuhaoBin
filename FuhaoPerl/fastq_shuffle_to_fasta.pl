#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::FastqKit qw/ShuffleFastq2Fasta/;

use constant USAGE =><<EOH;

usage: $0 R1.fq[.gz] R2.fq[.gz] output.fa[.gz]

v20171205

Description:
    Shuffle fastq to fasta
    Support gzipped input fastq and gzipped output fasta
    Needs Linux: zcat, gzip

Note: fastq head: 
	readname 1:N:0

	fastq support fq[.gz], fastq[.gz]
	fasta support fa, fas, fasta and gz format
	gz format need 'zcat or gzip' command

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

die "Error: failed to shuffle fastq to fasta\n" unless (ShuffleFastq2Fasta($ARGV[0], $ARGV[1], $ARGV[2]));
