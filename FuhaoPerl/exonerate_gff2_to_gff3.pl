#!/usr/bin/perl
use strict;
use warnings;
use constant USAGE =><<EOH;

 usage: $0 INFILE genestartID exonstartID > OUTFILE

 Description
   * print out the 'gene' line in GFF3
   * prints out the 'exon' lines in GFF3, link to parent gene

 INPUT FILE:   exonerate should have been run with (at least) the following options:
  --model coding2genome  OR  --model est2genome
  --showtargetgff yes

It is OK to have other output formats in the exonerate output

 v20170619

EOH

die USAGE if (scalar(@ARGV) < 1 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');
my $gffinput=$ARGV[0];
my $geneid = (defined $ARGV[1]) ? $ARGV[1]: "gene000000000";
my $exon_id = (defined $ARGV[2]) ? $ARGV[2]: "exon000000000";
my $est_acc;
print STDERR "\n###Info:\n#\tGFFinput: $gffinput\n#\tGeneIDstart: $geneid\n#\tExonIDstart: $exon_id\n";


# my $alignment_serial_number = 0;
open (GFFIN, "< $gffinput") || die "Error: can not open GFF input\n";
while (<GFFIN>) {
	chomp;
	if (/exonerate:\w+2genome/ && /\tgene\t/ && /sequence (\S+)/) {
		$est_acc = $1;
		$geneid++;
		my @cols = split /\t/, $_;
		my $identity=$cols[8];
		$cols[8] = "ID=$geneid;Name=$est_acc";
		if ($identity=~/;\s+identity\s+(\d+\.*\d*)/) {
			$cols[8]=$cols[8].";Identity=".$1;
		}
		if ($identity=~/\s+similarity\s+(\d+\.*\d*)/) {
			$cols[8]=$cols[8].";Similarity=".$1;
		}
		print join "\t", @cols;
		print "\n";
	}
	elsif (/exonerate:\w+2genome/ && /\texon\t/) {
		$exon_id++;

		my @cols = split /\t/, $_;
		$cols[8] = "ID=$exon_id;Parent=$geneid";

		print join "\t", @cols;
		print "\n";
	}
}
close GFFIN;
