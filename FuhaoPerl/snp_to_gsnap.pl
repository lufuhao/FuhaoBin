#!/usr/bin/env perl

use strict;
use warnings;
use Bio::DB::Fasta;
use constant USAGE =><<END;

USAGE: $0 <SNP-file> <reference-FASTA> | iit_store -o <prefix>.iit

Descriptions:
    Generates a SNP index formatted for use by GSNAP's iit_store

SNP-file
#chr    pos    Allele1    Allele2    ...

END
die USAGE unless (scalar(@ARGV)==2);
my $snpFile = $ARGV[0];
my $refFile = $ARGV[1];

my $FASTA = Bio::DB::Fasta->new($refFile);
my $lnum=0;

open (SNP, $snpFile) || die "Error: can not open snp file: $snpFile\n";
while (my $line=<SNP>) {
	chomp $line;
	$lnum++;
	next if ($line =~ /^#/);
	my @arr=split (/\t/, $line);
	if (scalar(@arr)<3) {
		print STDERR "Warnings: invalid columns at line($lnum): $line\n";
		next;
	}
	my $chr=shift @arr;
	my $pos=shift @arr;
	if ($pos =~ /,/ or $pos !~ /^\d+$/) {
		print STDERR "Warnings: invalid position at line($lnum): $line\n";
		next;
	}
	my $ref = $FASTA->seq( $chr, $pos => $pos );
	my %hash=();
	foreach my $allele (@arr) {
		$hash{$allele}=1;
	}
	@arr=();
	my @unique=keys(%hash);
	%hash=();
	foreach my $allele (@unique) {
		$allele='N' if (length($allele)>1);
		unless (length ("$ref$allele") == 2) {
			die "Error: invalid allele at line($lnum): $line\n     ref: $ref Allele: $allele\n";
		}
		if ($allele ne $ref) {
			print ">snp\t$chr:$pos\t$ref$allele\n";
		}
	}
}
close (SNP);



exit 0;
