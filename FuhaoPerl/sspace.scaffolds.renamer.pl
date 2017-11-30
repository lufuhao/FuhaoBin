#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FindBin qw($Bin);
use File::Basename;
use FuhaoPerl5Lib::FastaKit qw/SspaceOutRenamer/;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: LUFUHAO20171113

Requirements:
    Programs:
    Modiles: Scalar::Util, Cwd, Getopt::Long, FindBin

Descriptions:
    Rename SSPACE output fasta seqID to a certain pattern

Options:
    --help | -h
        Print this help/usage;
    --input | -h
        SSpace output fasta
    --source| -s
        
    --output | -o
        Renamed Fasta output [basename(\$input).rename.fa]
    --prefix | -p
        Prefix of new name [Default: MergedScaffold_]
    --suffix | -x
        Suffix of new name [Default: null]
    --base | -b
        Increment from this base [Default: 000000001]
    --verbose
   	Detailed output for trouble-shooting;
    --version | v!
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
my ($input, $scafold_source, $output, $name_prefix, $name_suffix, $name_base);

GetOptions(
	"help|h!" => \$help,
	"input|i:s" => \$input,
	"source|s:s" => \$scafold_source,
	"output|o:s" => \$output,
    "prefix|p:s" => \$name_prefix,
    "suffix|x:s" => \$name_suffix,
    "base|b:s" => \$name_base,
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$version) or die USAGE;
($help or $version) and die USAGE;



### Defaults ########################################################
$debug=0 unless (defined $debug);
$verbose=0 unless (defined $verbose);
$name_prefix="MergedScaffold_" unless (defined $name_prefix);
$name_suffix='' unless (defined $name_suffix);
$name_base='000000001' unless (defined $name_base);
$output=basename($input).".rename.fa" unless (defined $output);


### input and output ################################################
unless (defined $input and -s $input) {
	die "Error: invalid fasta input\n";
}
unless (defined $output) {
	die "Error: invalid fasta output\n";
}
unlink $output if (-e $output);



### Main ############################################################
unless (SspaceOutRenamer($input, $output, $name_base, $name_prefix, $name_suffix)) {
	die "Error: rename error\n";
}
