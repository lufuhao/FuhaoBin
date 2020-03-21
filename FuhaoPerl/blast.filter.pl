#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use constant USAGE=><<EOH;
SYNOPSIS:

perl $0 --input my.fa [Options]

Version: 20200220

Requirements:
    Programs:
        Perl: Getopt::Long

Descriptions:
    Determine the insert size given pairs of seqing data by
    mapping them to a reference.
Options:
    --help|-h
        Print this help/usage;
    --input|-i <FILE>
        Input file in BLAST outfmt 6 tabular format
    --identity | -p <Identity_percentage>
        Minimum identity percentage; default: 0
    --length | -l <length>
        Minimum alignment length; default: 0
    --evalue | -e <evalue>
        Maximum E-value; default:1
    --number | -n
        Maximum number of subjects for each query, default: 99999999
    --output | -o <outFile>
    --verbose
        Detailed output for trouble-shooting;
    --version|v!
        Print current SCRIPT version;

Example:
    perl $0 

Author:
    Fu-Hao Lu
    State Key Lab of Crop Stress Adaptation and Improvement
    Henan University
    Kaifeng 475004, P. R. China
    E-mail: lufuhao\@henu.edu.cn

EOH

###HELP ends#########################################################
die USAGE unless @ARGV;



###Receving parameter################################################
my ($input, $output);
my $help=0;
my $minlength=0;
my $minIdentity=0;
my $evalue=1;
my $maxnumber=99999999;
my $debug=0;
my $version=0;
my $verbose=0;

GetOptions(
	"help|h!" => \$help,
	"input|i:s" => \$input,
	"output|o:s" => \$output,
	"length|l:s" => \$minlength,
	"identity|p:f" => \$minIdentity,
	"evalue|e:s" => \$evalue,
	"number|n:i" => \$maxnumber,
	"output|o:s" => \$output,
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$version) or die USAGE;

($help or $version) and die USAGE;



### Defaults ########################################################
my %idnumber=();
$evalue=sprintf("%.10f", $evalue) if ($evalue=~/e/i);
unless (defined $output) {
	$output=$input."BlastFilter";
}



### input and output ################################################
print "#############  Configure #############\n";
print "## Input:       : $input\n";
print "## Output:      : $output\n";
print "## Min length   : $minlength\n";
print "## Min Identity : $minIdentity\n";
print "## Max E-value  : $evalue\n";
print "## Max Number   : $maxnumber\n";



### Main ############################################################
open (INPUT, "< $input") ||  die "Can not open input: $input\n";
open (OUTPUT, " > $output") || die "Can not write output: $output\n";

while (our $blast_line=<INPUT>) {
	if ($blast_line=~/^#/) {
		print OUTPUT $blast_line."\n";
		next;
	}

	chomp $blast_line;
	my @blast_comp=split(/\t/, $blast_line);
#0	1	2	3	4		5	6	7	8	9	10	11		12	13
#qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore	qcovs	qcovhsp
	my $test_print=1;
#control identity
	unless ($blast_comp[2]>=$minIdentity) {
		$test_print=0;
	}
#control length
	unless ($blast_comp[3]>=$minlength) {
		$test_print=0;
	}
#control evalue
	unless ($blast_comp[10]<=$evalue) {
		$test_print=0;
	}
#
	if ($test_print) {
		unless (exists $idnumber{$blast_comp[0]}) {
			%idnumber=();
		}
		$idnumber{$blast_comp[0]}{$blast_comp[1]}++;
		if (scalar (keys %{$idnumber{$blast_comp[0]}})<=$maxnumber) {
			print OUTPUT $blast_line."\n";
		}
	}
}
close OUTPUT;
close INPUT;




#####################################################################
###                         sub functions                         ###
#####################################################################
### ReadSam
###&ReadSam(sam,ref, 1/2/3)
###Global:
###Dependency:
###Note:



