#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::FastaKit qw/IndexFasta ExtractFastaSamtoolsID RunEmbossStretcher AnalysisEmbossStretcherOutput/;
use constant USAGE =><<EOH;

usage: $0 fasta id output

Requirements:
	samtools
	EMBOSS stretcher


v20160809

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


my $fasta=$ARGV[0];
my $idfile=$ARGV[1];
my $output=$ARGV[2];
my $path_samtools='samtools';
my $seq1='TEMP_seq1.fa';
my $seq2='TEMP_seq2.fa';
my $stretcherout='TEMP_seq1.vs.TEMP_seq2.stretcherout';
my $linenum=0;



die "Error: index fasta failed\n" unless (IndexFasta($fasta, $path_samtools)),

open (IDFILES, "<$idfile") || die "Error: can not open fasta Id file\n";
open (OUTPUT, "> $output") || die "Error: can not write output\n";
while (my $line=<IDFILES>) {
	chomp $line;
	$linenum++;
	&CleanTemp();
	my ($clusterno, $seq1name, $seq2name)=split(/\t/, $line);
	unless (ExtractFastaSamtoolsID($fasta, $seq1, $seq1name)) {
		print STDERR "Warnings: Line($linenum): Cluster($clusterno): failed to extract $seq1name\n";
		print OUTPUT "$clusterno\t$seq1name\t$seq2name\t#\t#\t#\t#\tError:seq1\n";
		next;
	}
	unless (ExtractFastaSamtoolsID($fasta, $seq2, $seq2name)) {
		print STDERR "Warnings: Line($linenum): Cluster($clusterno): failed to extract $seq2name\n";
		print OUTPUT "$clusterno\t$seq1name\t$seq2name\t#\t#\t#\t#\tError:seq2\n";
		next;
	}
	unless (RunEmbossStretcher($seq1, $seq2, 'nucl', $stretcherout)) {
		print STDERR "Warnings: Line($linenum): Cluster($clusterno): failed to run stretcher: $seq1name, $seq2name\n";
		print OUTPUT "$clusterno\t$seq1name\t$seq2name\t#\t#\t#\t#\tstretcher\n";
		next;
	}
	my ($test, $identity, $similarity, $gaps, $scores)=AnalysisEmbossStretcherOutput($stretcherout);
	unless ($test) {
		print STDERR "Warnings: Line($linenum): Cluster($clusterno): failed to get identity: $seq1name, $seq2name\n";
		print OUTPUT "$clusterno\t$seq1name\t$seq2name\t#\t#\t#\t#\tstretcher_out\n";
		next;
	}
	print OUTPUT "$clusterno\t$seq1name\t$seq2name\t$identity\t$similarity\t$gaps\t$scores\tsuccess\n";
}
&CleanTemp();
close IDFILES;
close OUTPUT;


sub CleanTemp {
	unlink $seq1 if (-e $seq1);
	unlink $seq2 if (-e $seq2);
	unlink $stretcherout if (-e $stretcherout);
}
