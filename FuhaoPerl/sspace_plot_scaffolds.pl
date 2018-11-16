#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::FastaKit qw/IndexFasta NumSeq/;
use FuhaoPerl5Lib::CmdKit;
use constant USAGE =><<EOH;

usage: $0 seqid.group original.fasta scaffolded.fasta reference.fasta outputprefix

Requirements
    seqextract_samtools.pl
    mum.stat

View SSPACE scaffolding joining using MUMmerplot
    seqid.group
        scaffold1	contig1	contig2	contig3
    sspace.input.fasta
    final.fasta
    reference.fasta
    outputprefix


v20160927

EOH
die USAGE if (scalar(@ARGV) !=5 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



die "Error: invalid evidence input\n" unless (defined "$ARGV[0]" and -s "$ARGV[0]");
my $inputids=$ARGV[0];
die "Error: invalid origin fasta\n" unless (defined "$ARGV[1]" and -s "$ARGV[1]");
my $originfasta=$ARGV[1];
die "Error: invalid scaffolded fasta\n" unless (defined "$ARGV[2]" and -s "$ARGV[2]");
my $scaffoldedfasta=$ARGV[2];
die "Error: invalid scaffolded fasta\n" unless (defined "$ARGV[3]" and -s "$ARGV[3]");
my $referenceseq=$ARGV[3];
die "Error: invalid output\n" unless (defined $ARGV[4]);
my $outputprefix=$ARGV[4];
my $linenum=0;
my $seqnum=0;
my %idhash=();



unless (IndexFasta($originfasta)) {
	die "Error: create fasta index failed using CDBfasta: $originfasta\n";
}
unless (IndexFasta($scaffoldedfasta)) {
	die "Error: create fasta index failed using CDBfasta: $scaffoldedfasta\n";
}
unless (IndexFasta($referenceseq)) {
	die "Error: create fasta index failed using CDBfasta: $scaffoldedfasta\n";
}


my $scaffoldfile="scaffolds.id";
my $contigsfile="contigs.id";
my $scaffoldnum=0;
open (IDGROUP, "< $inputids") || die "Error: can not open ID input: $inputids\n";
while (my $line=<IDGROUP>) {
	chomp $line;
	$linenum++;
	my @arr=split(/\t/, $line);
	die "Error: number elements <2\n" unless (scalar(@arr)>1);
	my $scaffoldid=shift @arr;
	
	unlink $scaffoldfile if (-e $scaffoldfile);
	close SCAFFOLDS if (defined fileno(SCAFFOLDS));
	open (SCAFFOLDS, "> $scaffoldfile") || die "Error: can not write $scaffoldfile\n";
	print SCAFFOLDS $scaffoldid, "\n";
	close SCAFFOLDS;
	unless (exec_cmd_return("seqextract_samtools.pl $scaffoldedfasta $scaffoldfile $scaffoldfile.fa")) {
		die "Error: scaffolds extraction running failed: $scaffoldfile.fa\n";
	}
	unless (-s "$scaffoldfile.fa") {
		die "Error: scaffolds extraction output failed: $scaffoldfile.fa\n";
	}
	my $number_of_seq=NumSeq("$scaffoldfile.fa");
	unless ($number_of_seq == 1) {
		die "Error: invalid scaffolds extraction number: $scaffoldfile.fa\n";
	}
	unlink $scaffoldfile;
	unless (exec_cmd_return("mum.stat -r $referenceseq -q $scaffoldfile.fa -o $scaffoldfile.fa.$scaffoldnum.after")) {
		die "Error: mum.stat failed: $scaffoldfile.fa\n";
	}
	
	unlink $contigsfile unless (-e $contigsfile);
	close CONTIGS if (defined fileno(CONTIGS));
	open (CONTIGS, "> $contigsfile") || die "Error: can not write $contigsfile\n";
	foreach my $id (@arr) {
		print CONTIGS $id, "\n";
	}
	close CONTIGS;
	unless (exec_cmd_return("seqextract_samtools.pl $originfasta $contigsfile $contigsfile.fa")) {
		die "Error: scaffolds extraction failed: $contigsfile.fa\n";
	}
	$number_of_seq=NumSeq("$contigsfile.fa");
	unless ($number_of_seq == scalar(@arr)) {
		die "Error: invalid scaffolds extraction number: $contigsfile.fa\n";
	}
	unlink $contigsfile;
	unless (exec_cmd_return("mum.stat -r $referenceseq -q $contigsfile.fa -o $contigsfile.fa.$scaffoldnum.before")) {
		die "Error: mum.stat failed: $contigsfile.fa\n";
	}
	
	
	$scaffoldnum++;
}
close IDGROUP;



