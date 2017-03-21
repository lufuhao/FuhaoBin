#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use constant USAGE =><<EOH;

usage: $0 infile.gb ourdir

Convert seq from GenBank to split EMBL using Bio::SeqIO

v20170309

EOH
die USAGE if (scalar(@ARGV) !=2 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

my $genbankfile=$ARGV[0];
my $outputdir=$ARGV[1];
$outputdir=~s/\/$//;


my $seqio = Bio::SeqIO->new('-format' => 'genbank', '-file' => "$genbankfile");

while( my $seq = $seqio->next_seq) {
	my $seqid=$seq->id;
	my $seqout = new Bio::SeqIO('-format' => 'embl', '-file' => ">$outputdir/$seqid.embl");
	$seqout->write_seq($seq)
}
