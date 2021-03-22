#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::FastaKit qw/FastaSpliterByBed/;
use constant USAGE =><<EOH;

usage: \$0 in.fas in.bed out.fas[.gz] seq_line_len

Version: 20210322

Descriptions
    Split Fasta long sequence into shorter parts
    Especially for mapping (chr > 530M)

EOH
die USAGE if (scalar(@ARGV) != 4);

my $retcode=FastaSpliterByBed($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3]);
unless ($retcode) {
	print STDERR "Error: failed to split\n";
	exit(100)
}

exit(0)
