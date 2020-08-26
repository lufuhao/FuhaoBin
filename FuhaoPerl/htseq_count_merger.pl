#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FindBin qw($Bin);
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: v20200825

Descriptions:
    collect HTseq xxx.count files into one matrix file

Options:
    --help|-h
        Print this help/usage;
    --input|-i
        Input .count file list, sep by comma， file name ending with .count
    --dir|-d
        Input .count file path
    --output|-o
        Count matrix output
    --verbose
        Detailed output for trouble-shooting
    --version|v!
        Print current SCRIPT version

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
my $input="";
my $output="MyOuput.count.matrix";
my $indir="";

GetOptions(
	"help|h!" => \$help,
	"input|i:s" => \$input,
	"output|o:s" => \$output,
	"dir|d:s" >= \$indir,
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$version) or die USAGE;
($help or $version) and die USAGE;



### Defaults ########################################################
$debug=0 unless (defined $debug);
$verbose=0 unless (defined $verbose);



### input and output ################################################
my @count_file_list=();
my @count_file_list1=();
my @count_file_list2=();
unless ($input eq "") {
	@count_file_list1=split(/,/, $input);
}
unless ($indir eq "") {
	@count_file_list2=glob "$indir/*";
}
if (scalar(@count_file_list1)>0 and scalar(@count_file_list1)>0) {
	@count_file_list=(@count_file_list1, @count_file_list2);
}
elsif (scalar(@count_file_list1)>0) {
	@count_file_list=(@count_file_list1);
}
elsif (scalar(@count_file_list2)>0) {
	@count_file_list=(@count_file_list2);
}
else {
	die "Error: empty file list\n";
}



### Main ############################################################
my $header;
my @sample;
my %hash;

foreach my $file (@count_file_list) {
	if ($file =~ /^\w+.*\.count/) {
		push @sample, $file;
		$header.= "\t$file";
		open (INPUT, "<", $file) or die "Error: can not open file: $file\n";
		while (my $line=<INPUT>) {
			chomp $line;
			next if ($line =~ /^\W+/);
			my @array = split (/\t/, $line);
			$hash{$array[0]} -> {$file} = $array[1];
		}
		close INPUT;
	}
}

open (OUTPUT, ">", $output) || die "Error: can not ";
print OUTPUT "$header\n";
foreach my $gene (sort keys %hash) {
    print OUTPUT "$gene";
    foreach my $file (@sample) {
        print OUTPUT "\t".$hash{$gene} -> {$file};
    }
    print OUTPUT "\n";
}
close OUTPUT;



