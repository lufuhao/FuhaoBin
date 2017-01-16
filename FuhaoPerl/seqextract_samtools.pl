#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 all.fasta ids outfile

Descriptions:
  Extract fasta using samtools

Requirements:
  Samtools

v20170116

#ID file
seq0001
seq0002[spaceORtab]start[spaceORtab]end
seq0003:start-end

#Note: start and end are 1-based

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

my $fasta=$ARGV[0];
die "Error: invalid fasta file input\n" unless (defined $fasta and -s $fasta);
my $idfile=$ARGV[1];
die "Error: invalid ID/region file input\n" unless (defined $idfile and -s $idfile);
my $faout=$ARGV[2];
if (-e $faout) {
	print STDERR "Warnings: $0 existing output $faout deleted\n";
	unlink $faout;
}

my $path2samtools='samtools';

unless (-e "$fasta.fai") {
	if (system("$path2samtools faidx $fasta")) {
		die "Error: fasta input index failed\n";
	}
}


my $fanum=0;
open(POSITIONS,"< $idfile") || die "Error: can not open fasta file\n";
open (FAOUTPUT, "> $faout") || die "Error: can not write fasta file\n";
while(my $line=<POSITIONS>){
	$fanum++;
	chomp $line;
	my ($seqName,$begin,$end)=('', 0, 0);
	if ($line=~/^(\S+)\s+(\d+)\s+(\d+)$/) {
		$seqName=$1;
		$begin=$2;
		$end=$3;
	}
	elsif ($line=~/^(\S+):(\d+)-(\d+)$/) {
		$seqName=$1;
		$begin=$2;
		$end=$3;
	}
	elsif ($line=~/^(\S+)$/) {
		$seqName=$1;
	}
	else {
		print STDERR "Error: invalid ID line($fanum): $line in ID file\n";
		next;
	}
	
	my $region='';
	if ($begin==0 and $end==0) {
		$region=$seqName;
	}
	elsif ($begin=~/^\d+$/ and $end=~/^\d+$/ and $begin>0 and $end>0 and $begin<=$end) {
		$region="$seqName:$begin-$end";
	}
	else {
		print STDERR "Error: unknown ID line($.): $line in ID file\n";
		next;
	}
	open(SAMTOOLSOUT,"$path2samtools faidx $fasta '$region' 2> /dev/null |");
	my $outnum=0;
	my $outline='';
	while(my $out = <SAMTOOLSOUT>){
		$outline.=$out;
		$outnum++;
		if ($outnum>=2) {
			print FAOUTPUT $outline;
			$outline='';
		}
	}
	close SAMTOOLSOUT;
	
	if ($outnum<2) {
		print STDERR "Error: failed to extract line($fanum): $line\n";
		next;
	}
	
}
close POSITIONS;
close FAOUTPUT;
