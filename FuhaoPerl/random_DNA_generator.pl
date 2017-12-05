#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;
use FuhaoPerl5Lib::FastaKit qw/RandomDNAgenerator/;
use constant USAGE=><<EOH;

SYNOPSIS:
    perl $0 -l length -n num

V20171205

Descriptions:
    Generate random DNA sequences for test

Options:
    -h  Print this help/usage;
    -l  length of random sequence [1000]
    -g  GC content of random sequence [0.5]
    -n  number of random sequence requested [1]
    -o  output full [PATH/]name [Random_sequence.fasta]


Example:
	perl $0 -l 1000 -n 1

Author:
    Fu-Hao Lu
    Post-Doctoral Scientist in Micheal Bevan laboratory
    Cell and Developmental Department, John Innes Centre
    Norwich NR4 7UH, United Kingdom
    E-mail: Fu-Hao.Lu\@jic.ac.uk

EOH
die USAGE unless (scalar(@ARGV)>0);

our($opt_g,$opt_l,$seq_length,$GC_value,$seq_number,$opt_n,$opt_v,$opt_h,$opt_o,$output);
getopts('hg:l:n:o:');
die USAGE if (defined $opt_h);


# set GC value and sequence length
$seq_length=((defined $opt_l) ? $opt_l : 1000);
$GC_value=((defined $opt_g) ? $opt_g : 0.5);
$seq_number=((defined $opt_n) ? $opt_n : 1);
$output=((defined $opt_o) ? $opt_o : "Random_sequence.fasta");

unless (RandomDNAgenerator($output, $seq_number, $seq_length, $GC_value)) {
	die "Error: RandomDNAgenerator running error\n";
}
