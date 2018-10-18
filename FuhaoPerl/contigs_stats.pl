#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;


###HELP################################
use constant USAGE =><<END;

SYNOPSIS:

perl fasta_stat.pl sequences.file [Options]
Version 20150519

Sequences.file	Sequences file in a format of fasta,
		genbank, scf, pir, embl, raw, gcg, ace, 
		bsml, swissprot, fastq and phd/phred;
		
		Note: Be careful to your seq file extension name
		It will try to guess the file format using file 
		extension name if you do not specify it using 
		"--format" or "-f";
		
Options:
	--help/-h
		Prints help/usage;
	--input/-i <File>
		[Msg] Sequence files in fasta.
	--min_length/-l <INT>
		The minimum sequence length you want to 
		count in the statistics (default:0);
	--format/-f <STR>
		fasta, genbank, scf, pir, embl, raw, gcg,
		ace, bsml, swissprot, fastq and phd/phred;
	--verbose/-v
		Detailed output for trouble-shooting;
	--version
		Print current SCRIPT version;

Example:
	perl fasta_stat.pl -i contigs.fa -l 500 -f fasta

Author:
	Fu-Hao Lu
	Post-Doctoral Scientist in Micheal Bevan laboratory
	Cell and Developmental Department, John Innes Centre
	Norwich NR4 7UH, United Kingdom
	E-mail: Fu-Hao.Lu\@jic.ac.uk
END
###HELP ends##########################
die USAGE unless @ARGV;



###Receving parameter#################
my ($help, $input, $min_length, $format, $verbose, $version);

GetOptions(
	"help|h!" => \$help,
	"input|i:s" => \$input,
	"min_length|l:i" => \$min_length,
	"format|f:s" => \$format,
	"verbose|v!" => \$verbose,
	"version!" => \$version) or die USAGE;

($help or $version) and die USAGE;
###Receving parameter ends############



###INPUT and DeFault##################
die "Error: invalid sequence file: $input\n" unless (-s $input);
if (defined $min_length) {
	chomp $min_length;
}
else {
	$min_length=0;
}
if (-e $input) {
	if (defined $verbose and $verbose>0) {
		print "The file to be assessed: $input\n" if (-e $input) || die "Can not input file?\n";
	}
}
else {die "Can not find input file\n";}
if (defined $format) {
	chomp $format;
}
###INPUT and Default ends############


###Guess Format######################
sub guess_format {
	my $s = shift @_;
#	print $s."\n";     #For TEST
	$s =~ s/^.*\.(.+)$/$1/;
#	print $s."\n";   #For TEST
	my $failed = 0;
	my $seq_format;
 	SW: {
	if ($s =~ /(^fasta)|(^fast)|(^fst)|(^fsa)|(^ft)|(^fs)|(^fa)|(^fas)/i) {$seq_format = 'Fasta'; last SW};
	if ($s =~ /(^fastq)|(^fq)/i) {$seq_format = 'Fastq'; last SW};
	if ($s =~ /(lfasta)|(lfast)|(lfst)|(lfsa)|(lft)|(lfs)/i) {$seq_format = 'LabeledFasta'; last SW};
	if ($s =~ /(embl)|(emb)|(em)|(eml)/i) {$seq_format = 'EMBL'; last SW};
	if ($s =~ /(genebank)|(genbank)|(genb)|(geneb)|(gbank)|(gb)/i) {$seq_format = 'GenBank'; last SW};
	if ($s =~ /(swissprot)|(sprt)|(swissp)|(sprot)|(sp)|(spr)/i) {$seq_format = 'Swissprot'; last SW};
	if ($s =~ /pir/i) {$seq_format = 'PIR'; last SW};
	if ($s =~ /gcg/i) {$seq_format = 'GCG'; last SW};
	if ($s =~ /scf/i) {$seq_format = 'SCF'; last SW};
	if ($s =~ /ace/i) {$seq_format = 'Ace'; last SW};
	if ($s =~ /phd/i) {$seq_format = 'phd'; last SW};
	if ($s =~ /phred|phd/i) {$seq_format = 'phred'; last SW};
	if ($s =~ /raw/i) {$seq_format = 'raw'; last SW};
	if ($s =~ /bsml/i) {$seq_format = 'bsml'; last SW};
	$failed++;
	}
	return eval{$failed ? 0 : $seq_format};
}


our ($file_basename, $file_dirname);
our $file_fullname="";
if ($input=~ m/\//g) {
	($file_basename=$input)=~ s/.*\///s;
	($file_dirname=$input)=~ s/(.*)\/.*$/$1/s;
	#$input=~ /.*\/(.*)$/; our $file_basename=$1;
	#$input=~ /(.*)\/.*$/; our $file_dirname=$1;
}
else {
	$file_basename=$input;
	$file_dirname=$ENV{'PWD'};
}
$file_fullname=$file_dirname."/".$file_basename;
#print $input."\n".$file_dirname."\n".$file_basename."\n".$file_fullname."\n";    #For test
if (-e $file_fullname) {
	print "Processing file: $file_fullname\n" if (defined $verbose and $verbose>0);
}
else{
	die "\nFile name conversion failed.\nPlease check the file input file path and name.\n";
}
if (! defined $format) {
#	print & guess_format($file_basename) || die "\nCan not guess file format\nPlease use --help for help\n";
	$format= & guess_format($file_basename) || die "\nCan not guess file format\nPlease use --help for help\n";
	print "The input format is supposed to be $format\n" if (defined $verbose and $verbose>0);
}
else {
	print "The input format is specified to be $format\n" if (defined $verbose and $verbose>0);
}

###Guess Format ends#################



###Get seq length array##############
our @seq_length=();
our @sorted_seq_length=();
our $longest_contig_id="";
our $longest_contig_length=0;
our $length_total=0;
#our $shortest_contig_id="";
#our $shortest_contig_length=0;

sub n_cal {
	my $n_cal_quartile = shift @_;
	my $n_cal_count=0;
	foreach my $n_cal_seq_length_ind (@sorted_seq_length) {
		$n_cal_count+=$n_cal_seq_length_ind;
		if ($n_cal_count >= $length_total*$n_cal_quartile/100){
			return ($n_cal_seq_length_ind);
			last;
		}
	}
}

sub num_contig {
	my $num_contig_quartile=shift @_;
	my $num_contig_count=0;
	foreach my $num_contig_seq_length_ind (@sorted_seq_length) {
		if ($num_contig_seq_length_ind >=$num_contig_quartile){
			$num_contig_count++;
		}
	}
	return ($num_contig_count);
}


our $seqio_obj=Bio::SeqIO->new(-file=>"$input", -format=>"$format");
while (our $seq_obj=$seqio_obj->next_seq) {
	if ($seq_obj->length > $min_length) {
		unshift @seq_length, $seq_obj->length;
		$length_total+=$seq_obj->length;
		if ($seq_obj->length > $longest_contig_length) {
			$longest_contig_id=$seq_obj->id;
			$longest_contig_length=$seq_obj->length;
		}
#		if ($seq_obj->length < $shortest_contig_length){
#			$shortest_contig_length=$seq_obj->length;
#			$shortest_contig_id=$seq_obj->id;
#		}
	}
}
@sorted_seq_length=sort {$b <=> $a} @seq_length;
print "N\tn:0\tn:200\tn:500\tn:n50\tmin\tn10\tn20\tn30\tn40\tn50\tn60\tn70\tn80\tn90\tmax\tmax_id\tfile\n";
print $length_total."\t".&num_contig(0)."\t".&num_contig(200)."\t".&num_contig(500)."\t".&num_contig(&n_cal(50))."\t".$sorted_seq_length[scalar(@sorted_seq_length)-1]."\t".&n_cal(10)."\t".&n_cal(20)."\t".&n_cal(30)."\t".&n_cal(40)."\t".&n_cal(50)."\t".&n_cal(60)."\t".&n_cal(70)."\t".&n_cal(80)."\t".&n_cal(90)."\t".$longest_contig_length."\t".$longest_contig_id."\t".$file_fullname."\n";
