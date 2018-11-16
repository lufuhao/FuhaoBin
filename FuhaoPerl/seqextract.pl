#!/usr/bin/env perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long;



###HELP################################
use constant USAGE=><<END;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: LUFUHAO20140807

		
Options:
	--help/-h
		Print this help/usage;
	--input|-i <FILE>
		[necessary] Input fasta file to extract from 
	--min_length/-l <INT>
		[Optional] The minimum sequence length you want to 
		count in the statistics (default:0);
	--seq_ids/-s <STR>
		[Optional] fasta IDs, comma-separated
	--seqid_file/-f <FILE>
		[Optional] File containing all sequence IDs, each 
		seq_id each line
	--reverse|-r
		Reverse selction against '--seq_ids' and/or 
		'--seqid_file', but not agaigst 'min_length'
	--output/-o <FILE>
		[Optional] New file for the extracted sequences
	--verbose
		[Optional] Detailed output for trouble-shooting;
	--version/-v
		[Optional] Print current SCRIPT version;

Example:
	perl $0 --input my.fasta --min_length 200 \
		-s id01,id02,id03 \
		--seqid_file myID.file
		--output Extracted.my.fa
		

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
our ($help, $input, $min_length, $seqid, @seqid_arr, $seqid_file, $reverse, $output, $verbose, $ver);

GetOptions(
	"help|h!" => \$help,
	"input|i=s" => \$input,
	"min_length|l:i" => \$min_length,
	"seq_ids|s:s" => \$seqid,
	"seqid_file|f:s" => \$seqid_file,
	"reverse|R|r!" => \$reverse,
	"output|o:s" => \$output,
	"verbose!" => \$verbose,
	"version|v!" => \$ver) or die USAGE;

($help or $ver) and die USAGE;

###DEFAULT
our ($input_basename, $input_base_no_ext)=('','');
($input_basename=$input)=~ s/.*\///s;
($input_base_no_ext=$input_basename)=~s/^(.*)\.\w+$/$1/;
$reverse=0 unless (defined $reverse);
#our $input_dirname=$input=~ s/(.*)\/.*$/$1/s;
$min_length=0 unless (defined $min_length);
@seqid_arr=();
unless (defined $output) {
	if ($reverse >0) {
		if ((defined $seqid) || (defined $seqid_file)) {
			$output=$input_base_no_ext.".ReverseSelected.fa";
		}
		else {
			die "Please specify the Sequence ID list by '--seq_ids' and/or '--seqid_file'\n";
		}
	}
	elsif ((defined $seqid) || (defined $seqid_file)) {
		$output=$input_base_no_ext.".extracted.fa";
	}
	else {
		$output=$input_base_no_ext.".filtered.fa"
	}
}
our ($output_basename, $output_base_no_ext, $sum_output);
($output_basename=$output)=~ s/.*\///s;
($output_base_no_ext=$output_basename)=~s/^(.*)\.\w+$/$1/;
$sum_output=$output_base_no_ext.'.extract.log';
print "INPUT: $input\nINPUT_basename: $input_basename\nINPUT_base_no_ext: $input_base_no_ext\nmin_length: $min_length\nReverse: $reverse\nOUTPUT: $output\nSum_output:$sum_output\n" if (defined $verbose);
###Receving parameter ends############



###Test###############################
if (defined $seqid) {
	if ($seqid=~m/,/) {
		@seqid_arr=split(/,/, $seqid);
	}
	else {
		push (@seqid_arr, $seqid);
	}
}
if (defined $seqid_file) {
	if (-e $seqid_file) {
		open (IDFILE, $seqid_file) || die "Can not open file specified by '--seqid_file/-sf'\n";
		while (<IDFILE>) {
			chomp $_;
#			print $_."\n";         ###test
			push (@seqid_arr, $_);
		}
		close IDFILE;
	}
	else {
		die "Can not open file specified by '--seqid_file/-sf'\n";
	}
}
#print "@seqid_arr\n";   ###test
if (-e $output) {
	if (defined $verbose) {
		print "The output $output specified already exists\nNeed to replease [y] or reinput [n], default[y]\n";
		my $output_replace_test=<STDIN>;
		chomp $output_replace_test;
		print $output_replace_test."\n";
		$output_replace_test='y' if ($output_replace_test eq '');
		if (lc $output_replace_test eq ('y'||'yes')) {
			unlink($output);
			`touch $output`;
		}
		else {
			die "Please specify another name for '--output'\n";
		}
	}
	else {
		unlink($output);
		`touch $output`;
	}
}
if (-e $sum_output) {
	if (defined $verbose) {
		print "The output $sum_output specified already exists\nNeed to replease [y] or reinput [n], default [y]:\n";
		my $sum_output_replace_test=<STDIN>;
		chomp $sum_output_replace_test;
		$sum_output_replace_test='y' if ($sum_output_replace_test eq '');
		if (lc $sum_output_replace_test eq ('y'||'yes')) {
			unlink($sum_output);
		}
		else {
			die "Please backup your $sum_output file then\n";
		}
	}
	else {
		unlink($sum_output);
	}
}



###ExtractSeq
our %found=();
our %failed=();
our %ignored=();
open (SUM, ">>$sum_output") || die "Can not write to NOTFOUNd file\n";
if ($reverse>0) {
	&ReverseExtractSelectedSeq();
	&PrintMessage();
}
elsif ((defined $seqid) || (defined $seqid_file)) {
	&ExtractSelectedSeq();
	&PrintMessage();
}
elsif (defined $min_length) {
	&FilterSeq();
	&PrintMessage();
}
else {
	die "Not input any filter ID or sequence length\n"
}
close SUM;





#######################################################################
###   sub functions   #################################################
#######################################################################

### Extract seq and write id not found to list
sub ExtractSelectedSeq {
	my $ESseqio_obj=Bio::SeqIO->new(-file=>$input, -format=>'fasta');
	my ($ESid, $ESseq_obj)=('','');
	my $ESoutput_obj=Bio::SeqIO->new(-file=>">$output", -format=>'fasta');
	while ($ESseq_obj=$ESseqio_obj->next_seq){
		foreach $ESid (@seqid_arr) {
			if ($ESseq_obj->id eq $ESid){
				if ($ESseq_obj->length >= $min_length) {
					$ESoutput_obj->write_seq($ESseq_obj);
					$found{$ESid}++;
				}
				else {
					$failed{$ESid}=$ESseq_obj->length;
				}
			}
		}
	}
}

###Simply filtre sequences with length
sub FilterSeq {
	my $FSseqio_obj=Bio::SeqIO->new(-file=>$input, -format=>'fasta');
	my $FSseq_obj='';
	my $FSoutput_obj=Bio::SeqIO->new(-file=>">$output", -format=>'fasta');
	while ($FSseq_obj=$FSseqio_obj->next_seq){
		if ($FSseq_obj->length >= $min_length) {
			$FSoutput_obj->write_seq($FSseq_obj);
			$found{$FSseq_obj->id}++;
		}
		else {
			$failed{$FSseq_obj->id}=$FSseq_obj->length;
		}
	}
}


###ReverseSelection
sub ReverseExtractSelectedSeq {
	my $REseqio_obj=Bio::SeqIO->new(-file=>$input, -format=>'fasta');
	my ($REid, $REseq_obj)=('','');
	my $REoutput_obj=Bio::SeqIO->new(-file=>">$output", -format=>'fasta');
	while ($REseq_obj=$REseqio_obj->next_seq){
		foreach $REid (@seqid_arr) {
			if ($REseq_obj->id eq $REid){
				$ignored{$REseq_obj->id}=$REseq_obj->length;
			}
		}
		if (not (exists $ignored{$REseq_obj->id})) {
			if ($REseq_obj->length >= $min_length) {
				$REoutput_obj->write_seq($REseq_obj);
				$found{$REseq_obj->id}++;
			}
			else {
				$failed{$REseq_obj->id}=$REseq_obj->length;
			}
		}
	}
}



sub PrintMessage {
	&PrintNotFound();
	&PrintFailed();
	&PrintIgnored();
	&PrintFound();
}


sub PrintNotFound {
	print SUM '#'x '30'."\nList of contigs not found:\n";
	print '#'x '30'."\nList of contigs not found:\n" if (defined $verbose);
	my $PNid='';
	foreach $PNid (@seqid_arr) {
		unless (exists $found{$PNid} or exists $failed{$PNid} or $ignored{$PNid} ) {
			print "$PNid ---------not found\n" if (defined $verbose);
			print  SUM $PNid."\n";
		}
	}
}



sub PrintFailed {
	print SUM "\n\n\n".'#' x '30'."\nList of failed contigs because the length less than $min_length:\n";
	if (%failed) {
		print "\n\n\n".'#' x '30'."\nList of failed contigs because the length less than $min_length:\n" if (defined $verbose);
		my @PFfailed_id_list=();
		@PFfailed_id_list=keys %failed;
		foreach (@PFfailed_id_list) {
			print "$_ --------------------Failed for <length\n" if (defined $verbose);
			print SUM $_."\t". $failed{$_}."\n";
		}
	}
}



sub PrintIgnored {
	print SUM "\n\n\n".'#' x '30'; print SUM "\nList of ignored contigs for reverse selection:\n";
	if (%ignored) {
		print "\n\n\n".'#' x '30'."\nList of ignored contigs for reverse selection:\n" if (defined $verbose);
		my @PIignored_id_list=();
		@PIignored_id_list=keys %ignored;
		foreach (@PIignored_id_list) {
			print "$_ --------------------ignored as reverse selection\n" if (defined $verbose);
			print SUM $_."\t". $ignored{$_}."\n";
		}
	}
}



sub PrintFound {
	print SUM "\n\n\n".'#' x '30'."\nList of found contigs (output to $output file):\n";
	if (%found) {
		print "\n\n\n".'#' x '30'."\nList of found contigs (output to $output file):\n" if (defined $verbose);
		my @Ffound_id_list=();
		@Ffound_id_list=keys %found;
		foreach (@Ffound_id_list) {
			print "$_ -------------------------------Found\n" if (defined $verbose);
			print SUM $_."\n";
		}
	}
}
