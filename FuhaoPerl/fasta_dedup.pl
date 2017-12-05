#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::FastaKit qw /FastaDedup/;
use constant USAGE =><<EOH;

usage: $0 input.fasta output.fasta

Deduplicate fasta based on 100% sequence similarity
And will warn if different sequence have the same seqID

v20170106

EOH
die USAGE if (scalar(@ARGV) !=2 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


unless (FastaDedup ($ARGV[0], $ARGV[1])) {
	die "Error: failed to run\n";
}
