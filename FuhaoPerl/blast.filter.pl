#!/usr/bin/perl
###usage: script.pl input minlength output
our $input=$ARGV[0];
chomp $input;
our $minlength=$ARGV[1];
chomp $minlength;
our $output=$ARGV[2];
chomp $output;

open (INPUT, "$input") ||  die "Can not input\n";
open (OUTPUT, ">>$output") || die "Can not output\n";

while (our $blast_line=<INPUT>) {
	chomp $blast_line;
	our @blast_comp=split(/\t/, $blast_line);
#0	1	2	3	4		5	6	7	8	9	10	11		12	13
#qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore	qcovs	qcovhsp
	if ($blast_comp[3]>=$minlength) {
		print OUTPUT $blast_line."\n";
	}
	else {
	next;
	}
}
close OUTPUT;
close INPUT;
