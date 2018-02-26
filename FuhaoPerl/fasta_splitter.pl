#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FindBin qw($Bin);
use FuhaoPerl5Lib::FastaKit qw/SplitFastaByNumber SplitFastaByLength ReadFastaLength/;
use Data::Dumper qw/Dumper/;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: LUFUHAO20180226

Requirements:
    Programs:
    Modiles: Scalar::Util, Cwd, Getopt::Long, FindBin

Descriptions:
    Determine the insert size given pairs of seqing data by
    mapping them to a reference.

Options:
    --help|-h
        Print this help/usage;
    --input|-i
        Fasta file ro be splited
    --num|-n
        Num of seq each file
    --length|-l
        Number of bases each file
    --prefix|-p
        Output file prefix
    --debug
        Debug mode
    --verbose
        Detailed output for trouble-shooting;
    --version|v!
        Print current SCRIPT version;

Example:
    perl $0 

Author:
    Fu-Hao Lu
    Post-Doctoral Scientist in Micheal Bevan laboratory
    Cell and Developmental Department, John Innes Centre
    Norwich NR4 7UH, United Kingdom
    E-mail: Fu-Hao.Lu\@jic.ac.uk

EOH
###HELP ends#########################################################
die USAGE unless @ARGV;



###Receving parameter################################################
my ($help, $verbose, $debug, $version);
my ($input, $prefix);
my $bynumber=0;
my $bylength=0;
GetOptions(
	"help|h!" => \$help,
	"input|i:s" => \$input,
	"prefix|p:s" => \$prefix,
	"num|n:i" => \$bynumber,
	"length|l:i" => \$bylength,
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$version) or die USAGE;
($help or $version) and die USAGE;



### Defaults ########################################################
$debug=0 unless (defined $debug);
$verbose=0 unless (defined $verbose);



### input and output ################################################
die "Error: invalid fasta input\n" unless (defined $input and -s $input);
if (defined $prefix and $prefix=~/^\S+$/) {
	print "### use defined prefix: $prefix\n";
}
else {
	$prefix=$input;
	$prefix=~s/^.*\///;
	print "### use default prefix: $prefix\n";
}



### Main ############################################################
if ((defined $bynumber and $bynumber=~/^\d+$/) or (defined $bylength and $bylength=~/^\d+$/)) {
	if ((defined $bynumber and $bynumber=~/^\d+$/ and $bynumber>0) and (defined $bylength and $bylength=~/^\d+$/ and $bylength>0)) {
		die "Error: can not specify both -n and -l\n";
	}
	elsif (defined $bynumber and $bynumber=~/^\d+$/ and $bynumber>0) {
		unless (SplitFastaByNumber($input, $bynumber, $prefix)) {
			die "Error: split failed\n";
		}
	}
	elsif (defined $bylength and $bylength=~/^\d+$/ and $bylength>0) {
		unless (SplitFastaByLength($input, $bylength, $prefix)) {
			die "Error: split failed\n";
		}
	}
}
else {
	die "Error: specify bynumber -n or bylength -l\n";
}
