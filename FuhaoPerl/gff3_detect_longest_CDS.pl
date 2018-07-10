#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::GffKit qw/ReadGff3 GuessLongestCDS WriteGff3/;
use constant USAGE =><<EOH;

usage: $0 in.gff3 in.fasta INT out.gff3

Guess the longest CDS given a GFF3 file, especially for pseudogenes

in.gff3   input GFF3 in flat format
in.fasta  multi-fasta file in input.gff3
INT       top INT longest CDS
out.gff3  output.gff3 with predicted CDS

v20180704

EOH
die USAGE if (scalar(@ARGV) !=4 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


my $inputgff3=shift @ARGV;
my $fastafile=shift @ARGV;
my $topnum=shift @ARGV;
my $outgff3=shift @ARGV;


my ($success1, $referenceids, $gene, $gene2mrna, $mrnas, $exons, $cds)=ReadGff3($inputgff3, $fastafile);

unless ($success1) {
	die "Error: ReadGff3 failed\n";
}

my ($success2, $return_gene2mrna, $return_mrna, $return_exon, $return_cdss)=GuessLongestCDS ($fastafile, $mrnas,  $exons, $topnum);

unless ($success2) {
	die "Error: GuessLongestCDS failed\n";
}

my $success3=WriteGff3($outgff3, $referenceids, $return_gene2mrna, $gene, $return_mrna, $return_exon, $return_cdss);
unless ($success3) {
	die "Error: WriteGff3 failed\n";
}
