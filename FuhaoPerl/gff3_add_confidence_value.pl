#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::GffKit qw/ReadGff3 WriteGff3/;
use constant USAGE =><<EOH;

usage: $0 in.gff3  in.fasta in.confidence  out.gff3

Descriptions:
    Read CONFIDENCE_VALUE from in.confidence
    add ONFIDENCE_VALUE to col6 of GFF3 file

in.gff3           input GFF3 in flat format
in.fasta          Multi-fasta file
in.confidence     Confidence file: 
                  col1=geneID	col2=confidence_value
out.gff3          output.gff3 with predicted CDS

v20180709

EOH
die USAGE if (scalar(@ARGV) !=4 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


my $inputgff3=shift @ARGV;
my $fastafile=shift @ARGV;
my $confid_file=shift @ARGV;
my $outgff3=shift @ARGV;

my %confidence_hash=();
my $linenum=0;


unless (open CONFIDENCEFILE, "<", $confid_file) {
	die "Error: can not open confidence value\n";
}
while (my $line=<CONFIDENCEFILE>) {
	$linenum++;
	chomp $line;
	next if ($line=~/^#/);
	my @conf_arr=split(/\t/, $line);
	unless (scalar(@conf_arr)>=2 and defined $conf_arr[0] and $conf_arr[0]=~/\S+/) {
		print STDERR "Warnings: invalid line($linenum): $line\n";
		next;
	}
	if (exists $confidence_hash{$conf_arr[0]}) {
		die "Error: duplicated ID: $conf_arr[0]\n" unless ($confidence_hash{$conf_arr[0]} eq $conf_arr[1]);
	}
	$confidence_hash{$conf_arr[0]} = $conf_arr[1];
}
close CONFIDENCEFILE;
print "### CONFIDENCE FILE SUM ###\n";
print "Total lines :       $linenum\n";
print "Total valid lines : ", scalar(keys %confidence_hash), "\n\n";


my ($success1, $referenceids, $gene, $gene2mrna, $mrnas, $exons, $cds)=ReadGff3($inputgff3, $fastafile);

unless ($success1) {
	die "Error: ReadGff3 failed\n";
}
foreach my $ind_geneid (sort keys %{$gene}) {
	if (exists $confidence_hash{$ind_geneid}) {
		${$gene}{$ind_geneid}{'score'}=$confidence_hash{$ind_geneid};
		if (exists ${$gene2mrna}{$ind_geneid}) {
			foreach my $ind_mRNA (sort keys %{${$gene2mrna}{$ind_geneid}}) {
				if (exists ${$mrnas}{$ind_mRNA}) {
					${$mrnas}{$ind_mRNA}{'score'}=$confidence_hash{$ind_geneid};
				}
				if (exists ${$exons}{$ind_mRNA}) {
					${$exons}{$ind_mRNA}{'score'}=$confidence_hash{$ind_geneid};
				}
				if (exists ${$cds}{$ind_mRNA}) {
					${$cds}{$ind_mRNA}{'score'}=$confidence_hash{$ind_geneid};
				}
			}
		}
	}
	else {
		print STDERR "Warnings: no value for gene ID: $ind_geneid\n";
	}
}


my $success3=WriteGff3($outgff3, $referenceids, $gene2mrna, $gene, $mrnas, $exons, $cds);
unless ($success3) {
	die "Error: WriteGff3 failed\n";
}
