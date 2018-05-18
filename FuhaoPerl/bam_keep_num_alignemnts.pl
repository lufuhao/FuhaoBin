#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::BamKit qw/BamKeepNumAlignemnts/;
use constant USAGE =><<EOH;

usage: $0 input.bam numK out.bam

v20180517

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


unless (BamKeepNumAlignemnts($ARGV[0], $ARGV[1], $ARGV[2])) {
	die "Error: BamKeepNumAlignemnts running failed\n";
}

