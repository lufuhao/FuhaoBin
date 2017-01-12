#! /usr/bin/perl -w
#
# File: solexa_parallelise_using_ssaha.pl
# Time-stamp: <19-Feb-2009 14:43:44 jit>
# $Id: $
#
# Copyright (C) 2008 by Pathogene Group, Sanger Center
#
# Author: JIT
# Description: a parallelised script to split the Illumina reads to subsets
#ftp://ftp.sanger.ac.uk/pub/project/pathogens/Uruguay_2013/Module_7_Annotation/backup/bam2augustusHints.pl
#modified by Fu-Hao Lu
#20150515
use strict;

my $PI = `echo $$` ;

if (@ARGV != 1 ) {
    print "$0 bam\n" ;
	exit;
}

my $cutoff = 10;
my $bam = shift;
my %intron_bam = () ; 
print STDERR "start parsing...\n" ;
open (IN, "samtools view $bam |") or die "oooops\n" ;
while (<IN>) {
    chomp;
    my @r=split (/\s+/, $_);
    next unless $r[5] =~ /N/;
    next unless $r[4] >= $cutoff;
	#print "\n$r[5]\n";
	my @bits = $r[5] =~ /(\d+\w)/g;
	#print "@bits\n";
	my $pos = $r[3];
	for(0..@bits-2) {
		if($_ == 0) {
			if ($bits[$_] =~ /(\d+)M/) {
				$pos += $1;
			}
			else {
				last;
			}
		}
		elsif($bits[$_] =~ /(\d+)M/) {
			my $end = $pos + $1 - 1;
			#print "$r[2]\tb2h\texon\t$pos\t0\t.\t.\tmult=
			$intron_bam{$r[2]}->{$pos."\t".$end}->{exon}++;
			$pos += $1;
		}
		elsif($bits[$_] =~ /(\d+)N/) {
			my $end = $pos + $1 - 1;
			$intron_bam{$r[2]}->{$pos."\t".$end}->{"intron"}++;
			$pos += $1;
		}
	}
}
close(IN) ; 
print STDERR "all done! start outputting...\n" ;

for my $chr (sort keys %intron_bam) {
    for my $pos ( sort keys %{ $intron_bam{$chr} } ) {
		for my $type (sort keys  %{ $intron_bam{$chr}->{$pos}}) {
			print "$chr\tb2h\t$type\t$pos\t0\t.\t.\tmult=$intron_bam{$chr}{$pos}->{$type}\;src=E\n" ;
		}
    }
}
# print "$bam.introns.hints.gff produced\n" ;
