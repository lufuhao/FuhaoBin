#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use FuhaoPerl5Lib::GffKit qw/Gff3Sort/;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 -i in.gff[.gz] -o out.gff3[.gz] -l chr.order.list
Version: v20210511

Requirements:
    Programs:
    Modiles: Getopt::Long

Descriptions:
    Sort GFF3 file by gene blocks
    Support gzipped input and output GFF3
    Need bgzip in PATH if out gzipped GFF3

Options:
    --help | -h
        Print this help/usage;
    --input | -i
        Input GFF3 file in flat or gzipped format
    --output | -o
        Input GFF3 file in flat or gzipped format
    --list | -l
        Chr list in desired order
    --verbose
        Detailed output for trouble-shooting;
    --version|v!
        Print current SCRIPT version;

Example:
    perl $0 

Author:
    Fu-Hao Lu
    Professor, PhD
    State Key Labortory of Crop Stress Adaptation and Improvement
    College of Life Science
    Jinming Campus, Henan University
    Kaifeng 475004, P.R.China
    E-mail: lufuhao\@henu.edu.cn
EOH
###HELP ends#########################################################
die USAGE unless @ARGV;



###Receving parameter################################################
my ($help, $verbose, $debug, $version);
my ($input, $output, $chr);

GetOptions(
	"help|h!" => \$help,
	"input|i:s" => \$input,
	"output|o:s" => \$output,
	"list|l:s" => \$chr,
#	!:s:i
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$version) or die USAGE;
($help or $version) and die USAGE;



### Defaults ########################################################
$debug=0 unless (defined $debug);
$verbose=0 unless (defined $verbose);



### input and output ################################################



### Main ############################################################
my $test=Gff3Sort($input, $output, $chr);
unless ($test) {
	print STDERR "Error: Gff3Sort failed\n";
}
else {
	print "Info: Gff3Sort succeed\n";
}

