#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::FastaKit qw/FastaSpliterByBed/;
use constant USAGE =><<EOH;

usage: \$0 in.fas in.bed out.fas[.gz] seq_line_len

Version: 20210323

Descriptions
    Split Fasta long sequence into shorter parts
    Especially for mapping (chr > 530M)

#in.fas
    in FLAT text

#in.bed: 6 columns as follows
1A_part1	0	293171175	1A	0	293171175
1A_part2	0	292095547	1A	293171175	585266722
1B_part1	0	341316541	1B	0	341316541
1B_part2	0	339795971	1B	341316541	681112512
2A_part1	0	387881130	2A	0	387881130
2A_part2	0	387567656	2A	387881130	775448786
2B_part1	0	395347580	2B	0	395347580

#out.fas[.gz]
    Out fasta, with gzip support [.gz]

#seq_line_len
    Fasta sequence line length, default 70

EOH
die USAGE if (scalar(@ARGV) != 4);

my $retcode=FastaSpliterByBed($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3]);
unless ($retcode) {
	print STDERR "Error: failed to split\n";
	exit(100)
}

exit(0)
