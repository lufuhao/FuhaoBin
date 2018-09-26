#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::GffKit qw/ReadGff3/;
use FuhaoPerl5Lib::FastaKit qw/ReadFastaLength2/;
use constant USAGE =><<EOH;

usage: $0 in.gff3 in.fasta out.pfx

NOTE: output might be overlapable
    use bedtools intersect or substract to make unique

v20180926

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

my ($success, $referenceids, $gene, $gene2mrna, $mrnas, $exons, $cds)=ReadGff3($ARGV[0]);


my $downstream_size=2000;
my $promoter_size=2000;
my %cds=();
my %promoter=();
my %downstream=();
my %utr5=();
my %utr3=();


my $seqfulllen = ReadFastaLength2($ARGV[1]);

foreach my $indgene (sort keys %{$gene}) {
	
	next unless (exists ${$gene2mrna}{$indgene});
	my @cdsarr=();
	my $gene_start=${$gene}{$indgene}{'start'};
	my $gene_end=${$gene}{$indgene}{'end'};
	my $strand=${$gene}{$indgene}{'strand'};
	my $ref=${$gene}{$indgene}{'reference'};
	unless (exists ${$seqfulllen}{$ref}) {
		die "Error: can not detect seq length: $ref\n";
	}
	foreach my $indmrna (sort keys ${$gene2mrna}{$indgene}) {
		next unless (exists ${$cds}{$indmrna} and exists ${$cds}{$indmrna}{'cds'});

		foreach my $left (keys %{${$cds}{$indmrna}{'cds'}}) {
			push (@cdsarr, $left);
			foreach my $right (keys %{${$cds}{$indmrna}{'cds'}{$left}}) {
				push (@cdsarr, $right);
			}
		}
	}
	if (scalar(@cdsarr)==0) {### pseudogene
		push (@cdsarr, $gene_start);
		push (@cdsarr, $gene_end);
	}
	@cdsarr=sort {$a<=>$b} @cdsarr;
	
	if ($strand eq '+') {
		$cds{$ref}{$cdsarr[0]}{$cdsarr[-1]}++;
		if ($gene_start<$cdsarr[0]) {
			$utr5{$ref}{$gene_start}{$cdsarr[0]-1}++;
		}
		if ($gene_end>$cdsarr[-1]) {
			$utr3{$ref}{$cdsarr[-1]+1}{$gene_end}++;
		}
		
		my $downright= $gene_end+$downstream_size;
		if ($downright > ${$seqfulllen}{$ref}) {
			$downright=${$seqfulllen}{$ref};
		}
		if ($downright>$gene_end) {
			$downstream{$ref}{$gene_end+1}{$downright}++;
		}
		my $upleft=$gene_start-$promoter_size;
		unless ($upleft>0) {
			$upleft=1;
		}
		if ($upleft<$gene_start) {
			$promoter{$ref}{$upleft}{$gene_start-1}++;
		}
	}
	elsif ($strand eq '-') {
		$cds{$ref}{$cdsarr[0]}{$cdsarr[-1]}++;
		if ($gene_start<$cdsarr[0]) {
			$utr3{$ref}{$gene_start}{$cdsarr[0]-1}++;
		}
		if ($gene_end>$cdsarr[-1]) {
			$utr5{$ref}{$cdsarr[-1]+1}{$gene_end}++;
		}
		
		my $upright= $gene_end+$promoter_size;
		if ($upright > ${$seqfulllen}{$ref}) {
			$upright=${$seqfulllen}{$ref};
		}
		if ($upright>$gene_end) {
			$promoter{$ref}{$gene_end+1}{$upright}++;
		}
		my $downleft=$gene_start-$downstream_size;
		unless ($downleft>0) {
			$downleft=1;
		}
		if ($downleft<$gene_start) {
			$downstream{$ref}{$downleft}{$gene_start-1}++;
		}
	}
}

open (UTR5, ">", $ARGV[2].".UTR5.bed") || die "Error: can not write UTR5\n";
foreach my $ref1 (sort keys %utr5) {
	foreach my $left1 (sort {$a<=>$b} keys %{$utr5{$ref1}}) {
		foreach my $right1 (sort {$a<=>$b} keys %{$utr5{$ref1}{$left1}}) {
			print UTR5 $ref1, "\t", $left1-1 , "\t", $right1, "\n";
		}
	}
}
close UTR5;
open (UTR3, ">", $ARGV[2].".UTR3.bed") || die "Error: can not write UTR3\n";
foreach my $ref1 (sort keys %utr3) {
	foreach my $left1 (sort {$a<=>$b} keys %{$utr3{$ref1}}) {
		foreach my $right1 (sort {$a<=>$b} keys %{$utr3{$ref1}{$left1}}) {
			print UTR3 $ref1, "\t", $left1-1 , "\t", $right1, "\n";
		}
	}
}
close UTR3;
open (DOWNSTREAM, ">", $ARGV[2].".downstream$downstream_size.bed") || die "Error: can not write DOWNSTRAM\n";
foreach my $ref1 (sort keys %downstream) {
	foreach my $left1 (sort {$a<=>$b} keys %{$downstream{$ref1}}) {
		foreach my $right1 (sort {$a<=>$b} keys %{$downstream{$ref1}{$left1}}) {
			print DOWNSTREAM $ref1, "\t", $left1-1 , "\t", $right1, "\n";
		}
	}
}
close DOWNSTREAM;
open (PROMOTER, ">", $ARGV[2].".Promoter$promoter_size.bed") || die "Error: can not write PROMOTER\n";
foreach my $ref1 (sort keys %promoter) {
	foreach my $left1 (sort {$a<=>$b} keys %{$promoter{$ref1}}) {
		foreach my $right1 (sort {$a<=>$b} keys %{$promoter{$ref1}{$left1}}) {
			print PROMOTER $ref1, "\t", $left1-1 , "\t", $right1, "\n";
		}
	}
}
close PROMOTER;
open (GENEBODY, ">", $ARGV[2].".GENEBODY.bed") || die "Error: can not write GENEBODY\n";
foreach my $ref1 (sort keys %cds) {
	foreach my $left1 (sort {$a<=>$b} keys %{$cds{$ref1}}) {
		foreach my $right1 (sort {$a<=>$b} keys %{$cds{$ref1}{$left1}}) {
			print GENEBODY $ref1, "\t", $left1-1 , "\t", $right1, "\n";
		}
	}
}
close GENEBODY;
