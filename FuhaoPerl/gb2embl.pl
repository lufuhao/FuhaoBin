#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use constant USAGE =><<EOH;

usage: $0 infile.gb outfile.embl

Convert seq from GenBank to EMBL using Bio::SeqIO

v20170309

EOH
die USAGE if (scalar(@ARGV) !=2 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

my $seqio = Bio::SeqIO->new('-format' => 'genbank', '-file' => "$ARGV[0]");
my $seqout = new Bio::SeqIO('-format' => 'embl', '-file' => ">$ARGV[1]");
while( my $seq = $seqio->next_seq) {
	$seqout->write_seq($seq)
}
