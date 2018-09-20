#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::GffKit qw/ReadGff3 WriteGff3 GffAddUTR/;
use Data::Dumper qw/Dumper/;
use constant USAGE =><<EOH;

usage: $0 in.gff3 in.fa out.gff3

Desc:
    Add UTR to GFF3

Requirements:
    FuhaoPerl5Lib::GffKit qw/ReadGff3 WriteGff3 GffAddUTR/

v20180920

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


my ($test1, $referenceids, $gene, $gene2mrna, $mrna, $exon, $cds)=ReadGff3($ARGV[0], $ARGV[1]);
unless ($test1) {
	die "Error: ReadGff3 error\n";
}
my ($test2, $utr) = GffAddUTR($mrna, $exon, $cds);
unless ($test2) {
	die "Error: GffAddUTR error\n";
}
#print Dumper $utr;
unless (WriteGff3 ($ARGV[2], $referenceids, $gene2mrna, $gene, $mrna, $exon, $cds, $utr)) {
	die "Error: WriteGff3 error\n";
}
