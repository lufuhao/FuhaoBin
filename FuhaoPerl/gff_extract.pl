#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 in.gff3 id_list 0/1 out.gff3

v20150716

id_list: list of geneID (ID=***;)
0/1: 0=keep list; 1=exclude list

EOH
die USAGE if (scalar(@ARGV) !=4 or [0] eq '-h' or [0] eq '--help');

my $inputgff=$ARGV[0];
my $idlist=$ARGV[1];
my $reversecode=$ARGV[2];
my $outputgff=$ARGV[3];
die "Error: invalid input GFF3\n" unless (defined $inputgff and -s $inputgff);
die "Error: invalid gene ID list file\n" unless (defined $idlist and -s $idlist);
die "Error: invalid output GFF3\n" unless (defined $outputgff);
unlink $outputgff if (-e $outputgff);



### Read Id list into hash
my $linenum=0;
my %geneid;
open (IDLIST, "<$idlist") || die "Error: can not open ID list file: $idlist\n";
while (my $line1=<IDLIST>) {
	$linenum++;
	chomp $line1;
	if ($line1=~/^(\S+)/) {
		if (exists $geneid{$1}){
			print STDERR "Warnings: repeated ID: $1\n";
		}
		else {
			$geneid{$1}++;
		}
	}
	else {
		print STDERR "Warnings: invalid gene ID $line1 at line: $linenum\n of $idlist";
	}
}
close IDLIST;
print "### Summary ###\n";
print "\tTotal ID lines: $linenum\n";
print "\tValid ID lines: ".scalar(keys %geneid)."\n";



### open In and out GFF
if ($inputgff=~/\.gz$/) {
	print "Info: input GFF in GFF.gz format\n";
	open (GFFIN, "zcat $inputgff | ") || die "Error: can not open compressed input GFF3\n";
}
elsif ($inputgff=~/\.gff\d*$/) {
	print "Info: input GFF in plain GFF format\n";
	open (GFFIN, "< $inputgff") || die "Error: can not open compressed input GFF3\n";
}
else {
	die "Error: can not guess input GFF format\n";
}
open (GFFOUT, ">$outputgff") || die "Error: can not write GFF output: $outputgff\n";



### Parse input GFF
my $testprint=0;
my %mrnaid=();
while (my $line2=<GFFIN>) {
	if ($line2=~/^#/) {
		print GFFOUT $line2;
		next;
	}
	
	chomp $line2;
	my @arr=split(/\t/, $line2);
	if ($arr[2] =~/gene/i) {
		if ($arr[8]=~/ID=([^; ]+)/) {
			my $geneid=$1;
			if (exists $geneid{$geneid}) {
				$testprint= ($reversecode==0) ? 1 : 0;
			}
			else {
				$testprint=($reversecode==1) ? 1 : 0;
			}
		}
	}
=old
	elsif ($arr[2] =~/mrna/i) {
		my $mrnaid1;
		if ($arr[8]=~/ID=([^; ]+)/) {
			$mrnaid1=$1;
		}
		if ($arr[8]=~/Parent=([^; ]+)/) {
			if (exists $geneid{$1}) {
				if ($reversecode==0) {
					$testprint=1;
				}
				else {
					$mrnaid{$mrnaid1}++;
					$testprint=0;
				}
			}
			else {
				if ($reversecode==1) {
					$testprint=1;
				}
				else {
					$mrnaid{$mrnaid1}++;
					$testprint=0;
				}
			}
		}
	}
	elsif ($arr[2] =~/(CDS)|(exon)/i) {
		my $parentid;
		if ($arr[8]=~/Parent=([^; ]+)/) {
			$parentid=$1;
			if (exists $geneid{$parentid}) {
				$testprint= ($reversecode==0) ? 1 : 0;
			}
			else {
				$testprint=($reversecode==1) ? 1 : 0;
			}
			$testprint=0 if (exists $mrnaid{$parentid});
		}
	}
=cut
	print GFFOUT $line2."\n" if ($testprint);
}
close GFFIN;
close GFFOUT;
