#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::BamKit qw/ExtractMultimapping/;
use constant USAGE =><<EOH;
usage: $0 in.bam out.readID outrefID

Desc
  Evaluate multi-mapping by counting mapping times
  output multi-mapped reads and their refs

Requirements
  samtools

v20180911

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

unless (ExtractMultimapping($ARGV[0], $ARGV[1], $ARGV[2])) {
	die "Error: ExtractMultimapping running error\n";
}
