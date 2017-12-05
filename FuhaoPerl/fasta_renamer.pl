#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;


###HELP################################
use constant USAGE =><<END;

SYNOPSIS:

 $0 sequences.file [Options]

Version 20171205

Sequences.file	Youe sequences file in a format of fasta,
    genbank, scf, pir, embl, raw, gcg, ace, 
    bsml, swissprot, fastq and phd/phred;

    Note: Be careful to your seq file extension name
      It will try to guess the file format using file 
      extension name if you do not specify it using 
      "--format" or "-fmt";

Options:
    --help|-h
        Prints help/usage;
    --input|-i [file]
        fasta file input fa[.gz]
    --prefix|-p <str>
        The prefix of ID you want [''];
    --num_digit|-n <int>
        Number of digits behind the prefix
        eg: 4 = PREFIX0000, PREFIX0001, PREFIX0002...
        PREFIX+originalID if NOT defined
    --output|-o <str>
        output file name, default "input_basename.Final.fa";
    --version|-v
        Print current SCRIPT version;
    --verbose
        Print detailed output on STDOUT;

Example:
    perl fasta_stat.pl contigs.fa --min_length 500 -fmt fasta

Author:
    Fu-Hao Lu
    Post-Doctoral Scientist in Micheal Bevan laboratory
    Cell and Developmental Department, John Innes Centre
    Norwich NR4 7UH, United Kingdom
    E-mail: Fu-Hao.Lu\@jic.ac.uk
END
###HELP ends##########################
die USAGE unless (@ARGV);


###Receving parameter#################
my ($help, $input, $output, $prefix, $num_digit, $format, $log, $verbose, $version);

GetOptions(
	"help|h!" => \$help,
	"input|i=s" => \$input,
	"prefix|p:s" => \$prefix,
	"num_digit|n:i" => \$num_digit,
	"output|o:s" => \$output,
	"log:s" => \$log,
	"verbose!" => \$verbose,
	"version|v!" => \$version,
) or die USAGE;
($help or $version) and die USAGE;
###Receving parameter ends############



###Defaults###########################
$prefix='' unless (defined $prefix);



###INPUT and DeFault##################
die "Can not fould the input file specified: $input\n" unless (defined $input and -s $input);
my ($input_basename, $input_base_no_ext)=('','');
($input_basename=$input)=~ s/.*\///s;
($input_base_no_ext=$input_basename)=~s/^(.*)\.\w+$/$1/;
$output=$input_base_no_ext.".Final.fa" unless (defined $output);
my ($output_dir, $output_basename, $output_base_no_ext)=('', '','');
unless (($output_dir=$output)=~ s/(.*)\/.*$/$1/s) {
	$output_dir='.';
}
($output_basename=$output)=~ s/.*\///s;
($output_base_no_ext=$output_basename)=~s/^(.*)\.\w+$/$1/;
$log=$output_dir.'/'.$output_base_no_ext.".log" unless (defined $log);
print "\n\n\n##### Summary #####\n";
print "Input: $input\n";
print "Output: $output\n";
print "LOG: $log\n\n\n";



if (-e $output) {
	if (defined $verbose and $verbose >0) {
		print "The output file name already exists, Need to overwite [y] or re-specify [n]? default[y]:\n";
		my $output_exist=<STDIN>;
		chomp $output_exist;
		$output_exist='y' unless ($output_exist eq '');
		if (lc($output_exist) eq ('y'||'yes')) {
			unlink ($output);
		}
		else {
			die "Please input new output file name.\n";
		}
	}
	else {
#		rename("$output", &backup($output_basename));
		unlink ($output);
		print "###warning:\nThe output file exists and is renamed\n###\n";
	}
}
#$num_digit=8 unless (defined $num_digit);



###INPUT and Default ends############
my $count=0;
chomp $input;
#print "$input\n".length($input)."\n";
if ($input=~/(\.fa$)|(\.fas$)|(\.fasta$)/i) {
	print STDERR "InFo: input in fasta format\n";
	open (INPUT, "$input") || die "Can not find the input\n";
}
elsif ($input=~/(\.gz$)|(\.gzip$)/i) {
	print STDERR "Info: input in gz format\n";
	open (INPUT, "zcat $input |") || die "Can not find the input\n";
}
else {
	die "Can not read input format\n";
}
my $input_line_number=0;
my $total_sequences=0;
open (OUTPUT, ">$output") || die "Error: can not output\n";
open (LOG, ">$log") || die "Error: can not write log\n";
while (my $input_line=<INPUT>) {
	chomp $input_line;
	$input_line_number++;
	my $original_ID='';
	my $newID='';
	if ($input_line=~m/^>(\S+)\s*/ ) {
		$original_ID=$1;
		$total_sequences++;
		die "Error: invalid fasta ID at line $input_line_number: $input_line\n" if ($input_line eq '' or $original_ID eq '');
		$count++;
		if (defined $num_digit) {
			my $new_8d=&Full8D($count);
			$newID=$prefix.$new_8d;
		}
		else {
			$newID=$prefix.$original_ID;
		}
		die "Error: invalid new ID at line $input_line_number: $input_line\n" if ($input_line eq '' or $original_ID eq '');
		$input_line=~s/^>\S+/>$newID/;
		print LOG $newID."\t".$original_ID."\t".$input_line."\n";
	}
	print OUTPUT $input_line."\n";
}
close OUTPUT;
close INPUT;
close LOG;
print "Number of sequences: $total_sequences\n";

####################################
###  sub functions##################
####################################
sub backup {
	my $BUori_name=shift @_;
	my $BU_new_name='';
	$BU_new_name=$BUori_name.".bak.".int(rand(10000));
	if (-e "$output_dir/$BU_new_name") {
		&backup($BUori_name);
	}
	else {
		return ("$output_dir/$BU_new_name");
	}
}



###
###Global: $num_digit
sub Full8D {
	my $FDnum=shift @_;
	my $FD_8D='';
	die "Wrong num_digit detected\n" unless (defined $num_digit and $num_digit>0);
	if (length($FDnum)<=$num_digit) {
		$FD_8D='0'x ($num_digit-length($FDnum)).$FDnum;
	}
	else {
		die "There are more than 1E".($num_digit+1)." sequences exists\n";
	}
	return $FD_8D;
}
