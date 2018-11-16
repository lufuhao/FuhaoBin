#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use constant USAGE =><<EOH;

usage: $0 in.genebank out.fasta

Desc:
	Convert sequence format from genebank to fasta
	
Require:
	BioPerl (Bio::SeqIO)
	
v20151102

EOH
die USAGE if (scalar(@ARGV) !=2 or [0] eq '-h' or [0] eq '--help');



my $genebankfile  = shift @ARGV;
chomp $genebankfile;
my $fastaout=shift @ARGV;
chomp $fastaout;
my $seqin  = Bio::SeqIO->new(-file => "$genebankfile" , '-format' => 'GenBank');
my $seqout = Bio::SeqIO->new(-file => ">$fastaout" ,'-format' => 'Fasta');
while(my $ind_seq = $seqin->next_seq() ){
	$seqout->write_seq($ind_seq);
}
