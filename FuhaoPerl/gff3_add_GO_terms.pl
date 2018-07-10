#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::GffKit qw/ReadGff3 WriteGff3/;
use constant USAGE =><<EOH;

usage: $0 in.gff3  in.fasta in.confidence  out.gff3

Guess the longest CDS given a GFF3 file, especially for pseudogenes

in.gff3           input GFF3 in flat format
in.fasta          Multi-fasta file
in.confidence     Confidence file: 
                  col1=geneID	col2=GO:0000000	col3=GO:000001 ....
out.gff3          output.gff3 with predicted CDS

v20180709

EOH
die USAGE if (scalar(@ARGV) !=4 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


my $inputgff3=shift @ARGV;
my $fastafile=shift @ARGV;
my $confid_file=shift @ARGV;
my $outgff3=shift @ARGV;

my %go_hash=();
my $linenum=0;


unless (open GOFILE, "<", $confid_file) {
	die "Error: can not open GO file\n";
}
while (my $line=<GOFILE>) {
	$linenum++;
	chomp $line;
	next if ($line=~/^#/);
	my @conf_arr=split(/\t/, $line);
	unless (scalar(@conf_arr)>=2 and defined $conf_arr[0] and $conf_arr[0]=~/\S+/) {
		print STDERR "Warnings: invalid line($linenum): $line\n";
		next;
	}
	my $geneID = shift @conf_arr;
	foreach my $x (@conf_arr) {
		$go_hash{$geneID}{$x}++;
	}
}
close GOFILE;
print "### GO FILE SUM ###\n";
print "Total lines :       $linenum\n";
print "Total valid lines : ", scalar(keys %go_hash), "\n\n";


my ($success1, $referenceids, $gene, $gene2mrna, $mrnas, $exons, $cds)=ReadGff3($inputgff3, $fastafile);

unless ($success1) {
	die "Error: ReadGff3 failed\n";
}
foreach my $ind_geneid (sort keys %{$gene}) {
	if (exists $go_hash{$ind_geneid}) {
		if (exists ${$gene}{$ind_geneid}{'Ontology_term'}) { ### merge if there is GO in input GFF files
			my @arr2=split(/,/, ${$gene}{$ind_geneid}{'Ontology_term'});
			foreach my $x (@arr2) {
				$go_hash{$ind_geneid}{$x}++;
			}
		}
		my @arr3=sort keys %{$go_hash{$ind_geneid}};
		${$gene}{$ind_geneid}{'Ontology_term'}=join (",",@arr3);
#		if (exists ${$gene2mrna}{$ind_geneid}) {
#			foreach my $ind_mRNA (sort keys %{${$gene2mrna}{$ind_geneid}}) {
#				if (exists ${$mrnas}{$ind_mRNA}) {
#					${$mrnas}{$ind_mRNA}{'Ontology_term'}=join (",",@arr3);
#				}
#				if (exists ${$exons}{$ind_mRNA}) {
#					${$exons}{$ind_mRNA}{'Ontology_term'}=join (",",@arr3);
#				}
#				if (exists ${$cds}{$ind_mRNA}) {
#					${$cds}{$ind_mRNA}{'Ontology_term'}=join (",",@arr3);
#				}
#			}
#		}
	}
}


my $success3=WriteGff3($outgff3, $referenceids, $gene2mrna, $gene, $mrnas, $exons, $cds);
unless ($success3) {
	die "Error: WriteGff3 failed\n";
}
