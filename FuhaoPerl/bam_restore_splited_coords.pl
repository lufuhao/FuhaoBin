#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::BamKit qw/BamRestoreSplit/;
use constant USAGE =><<EOH;

usage: $0 in.bam in.bed  out.bam

Requirements
	samtools
	
BED file (6columns)
        chr1A_part1     0       471304005       chr1A   0       471304005
        chr1A_part2     0       122798051       chr1A   471304005       594102056
        chr1B_part1     0       438720154       chr1B   0       438720154
        chr1B_part2     0       251131716       chr1B   438720154       689851870
        chr1D_part1     0       452179604       chr1D   0       452179604
        chr1D_part2     0       43273582        chr1D   452179604       495453186

v20180809

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


unless (BamRestoreSplit($ARGV[0], $ARGV[1], $ARGV[2], 'samtools')) {
	die "Main.Error: BamRestoreSplit running error\n";
}
