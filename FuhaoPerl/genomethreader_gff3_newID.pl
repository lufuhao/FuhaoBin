#!/usr/bin/perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 gth.gff3 > gth.out.gff3

Version: 20150714

Descriptions:
  Prints out the 'exon' lines in GFF3, ID=geneXXXXXXXXX;Target=source_ID
  It is OK to have other output formats in the genomethreader output

EOH

die USAGE if (scalar(@ARGV) != 1 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

my $est_acc;

my $geneid="gene000000000";
my $exon_id="exon000000000";
my $cds_id="cds000000000";
my $mrnaid="mrna000000000";

open (GTH, "< $ARGV[0]") || die "Error: can not open gth out\n";

while (<GTH>) {
	chomp;
	if (/\tgth\t/ && /\tgene\t/ && /ID=(\S+)/) {
		$est_acc = $1;
		$geneid++;
		my @cols = split (/\t/, $_);
		$cols[8] = "ID=$geneid;Name=$est_acc";
		print join "\t", @cols;
		print "\n";
	}
	elsif (/\tgth\t/ && /\tmRNA\t/) {
		$mrnaid++;
		my @cols = split (/\t/, $_);
		my $target='.';
		if (/Target=(.*)$/) {
			$target=$1;
		}
		$cols[8] = "ID=$mrnaid;Parent=$geneid;Target=$target";
		print join "\t", @cols;
		print "\n";
	}
	elsif (/\tgth\t/ && /\tCDS\t/) {
		$cds_id++;
		my @cols = split (/\t/, $_);
		$cols[8] = "ID=$cds_id;Parent=$mrnaid";
		print join "\t", @cols;
		print "\n";
	}
	elsif (/\tgth\t/ && /\texon\t/) {
		$exon_id++;
		my @cols = split (/\t/, $_);
		$cols[8] = "ID=$exon_id;Parent=$mrnaid";
		print join "\t", @cols;
		print "\n";
	}
}
close GTH;
