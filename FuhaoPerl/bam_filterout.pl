#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 [Options]
Version: LUFUHAO20150519

Requirements:
	Programs: samtools

Descriptions:
	Filter out some references in BAM

Options:
	--help|-h
		Print this help/usage;
	--input|-i	<File>
		[Msg] BAM/SAM input
	--idlist|-s	<list>
		[Opt] IDS to filter out, comma-delimited
	--idfile|-f <File>
		[Opt] Files contsining IDs, 1 ID/line
	--output|-o	<File>
		[Opt] BAM output, default: input_base.filterout.bam
	--filhead|-h
		[Opt] Also clean header
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
my ($input, $output, $idlist, $idfile, $filheader);

GetOptions(
	"help|h!" => \$help,
	"input|i:s" => \$input,
	"idlist|s:s" => \$idlist,
	"idfile|f:s" => \$idfile,
	"output|o:s" => \$output,
	"filhead|d!" => \$filheader,
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$version) or die USAGE;
($help or $version) and die USAGE;



### Defaults ########################################################
$debug=0 unless (defined $debug);
$verbose=0 unless (defined $verbose);
$output = &RetrvNoExt($input).'.filterout.bam' unless (defined $output);


### input and output ################################################
die "Error: invalid BAM input: $input\n" unless (defined $input and -s $input);
die "Error: existing output: $output\n" if (defined $output and -e $output);
die "Error: invalid seq ID file: $idfile\n" if (defined $idfile and ! -s $idfile);



### Main ############################################################
###Check Id list
my %idarr=();
my @ids=();
if (defined $idlist) {
	@ids=split(/,/, $idlist);
	map {$idarr{$_}++} @ids;
}
if (defined $idfile) {
	open (IDFILE, "<", $idfile) || die "Error: can not open ID file specified by --idfile/-f $idfile\n";
	while (my $line1=<IDFILE>) {
		chomp $line1;
		(my $ind_id=$line1)=~s/^(\S+)\s*.*$/$1/;
		$idarr{$ind_id}++;
	}
	close IDFILE;
}
print "Info: total numbers of IDs to be ecluded: ", scalar(keys %idarr), "\n" if ($verbose or $debug);



### Parse BAM and filter out
my %excluded_header=();
my %excluded_alignments=();
if ($input =~ /\.sam/i) {
	open (BAMIN, "samtools view -h -S $input | ") || die "Error: can not open SAM input: $input\n";
}
elsif ($input =~ /\.bam/i) {
	open (BAMIN, "samtools view -h $input | ") || die "Error: can not open BAM input: $input\n";
}
else {
	die "Error: can not guess input formar SAM/BAM: $input\n";
}
open (BAMOUT, " | samtools view -S -h -b - > $output") || die "Error: can not write BAM output: $output\n";
while (my $line2=<BAMIN>) {
	if ($line2=~/\@/) {
		if ($filheader) {
			if ($line2=~/\@SQ.*SN:(\S+)\s*.*$/) {
				if (exists $idarr{$1}) {
					$excluded_header{$1}++;
					next;
				}
			}
		}
		print BAMOUT $line2;
	}
	else {
		my @arr=split(/\t/, $line2);
		if (exists $idarr{$arr[2]}) {
			$excluded_alignments{$arr[2]}++;
			next;
		}
		else {
			print BAMOUT $line2;
		}
	}
}
close BAMIN;
close BAMOUT;
if ($verbose or $debug) {
	print "Info: Total excluded headers IDs: ", scalar(keys %excluded_header), "\n";
	print "Info: Total excluded alignments IDs: ", scalar(keys %excluded_alignments), "\n";
}
exit 0;


#####################################################################
###                         sub functions                         ###
#####################################################################
### Retrieve filebasename without extension
### &RetrvNoExt(file)
### Global: 
sub RetrvNoExt {
	my $RNE_ori=shift @_;
	chomp $RNE_ori;
	my $RNE_new='';
	my $RNE_base='';
	($RNE_base=$RNE_ori)=~s/.*\///s;
	($RNE_new=$RNE_base)=~s/^(\S+)\.\w+$/$1/;
	return $RNE_new;
}
