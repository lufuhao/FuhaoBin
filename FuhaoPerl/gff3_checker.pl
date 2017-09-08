#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::GffKit qw/ReadGff3 WriteGff3/;
use constant USAGE =><<EOH;

usage: $0 input.gff3 [MyOutput.gff3] [genome.fasta ]

 * check GFF3 for errors
 * Correct CDS phased
 * Autocomplete CDS phased
 * Fill in CDS phase if fasta exists

v20170628

EOH
die USAGE if (scalar(@ARGV) <1 or scalar(@ARGV)>3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

my $inputgff3=shift @ARGV;
unless (defined $inputgff3 and -s $inputgff3) {
	die "Error: invalid GFF3 input\n";
}
my $outputgff3=shift @ARGV;
$outputgff3="MyOutput.gff3" unless (defined $outputgff3);
unlink $outputgff3 if (-e $outputgff3);

my $dbfasta=shift @ARGV;
unless (defined $dbfasta and -s $dbfasta) {
	print "Info: No fasta file detect, would not disable phase check\n";
}

my ($success, $referenceids, $gene, $gene2mrna, $mrnas, $exons, $cds)=ReadGff3($inputgff3, $dbfasta);
unless ($success) {
	die "Error: ReadGff3 failed\n";
}

unless (WriteGff3($outputgff3, $referenceids, $gene2mrna, $gene, $mrnas, $exons, $cds)) {
	die "Error: WriteGff3 failed\n";
}
