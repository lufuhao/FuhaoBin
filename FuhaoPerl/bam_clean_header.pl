#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input in.bam --output out.bam
Version: LUFUHAO20150430

Requirements:
	Programs:
	Modiles: Scalar::Util, Cwd, Getopt::Long, FindBin

Descriptions:
	Determine the insert size given pairs of seqing data by
	mapping them to a reference.

Options:
	--help|-h
		Print this help/usage;
	--input|-i	<File>
		[Msg] Bam input
	--output|-o	<path/to/File.bam>
		[Opt] Bam output
	--out1rg|-n <path/to/File.bam>
		[Opt] Bam output without RG field
	--debug
		[Opt] Output detailed info
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
our ($help, $verbose, $debug, $ver);
our ($input, $output, $outnogroup);

GetOptions(
	"help|h!" => \$help,
	"input|i:s" => \$input,
	"output|o:s" => \$output,
	"out1rg|n:s" => \$outnogroup,
	"debug|" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$ver) or die USAGE;
($help or $ver) and die USAGE;



### Defaults ########################################################
$debug=0 unless (defined $debug);



### input and output ################################################
die "Error: invalid input\n" unless (-s $input);
unlink $output if (-e $output);
unlink $outnogroup if (defined $outnogroup and -e $outnogroup);


### Main ############################################################
&SamCleanHeader($input, $output, $outnogroup);



#####################################################################
###                         sub functions                         ###
#####################################################################
### Clean those sequences in header but not exist in alignments
###&SamCleanHeader(bam_input, bam_output)
###Global: $debug;
###Dependency: 
###Note: 
sub SamCleanHeader {
	my ($SCHbam_input, $SCHbam_output, $SCHoutnogroup)=@_;
	die "SUB(SamCleanHeader)Error: BAM input file not found\n" unless (defined $SCHbam_input and -s $SCHbam_input);
	my $SCHnew_readgroup=0;
	my %SCHreahgroup=();
	unless (defined $SCHbam_output and $SCHbam_output ne '' and $SCHbam_output !~ m/\s+/) {
		my $SCHbam_prefix=&RetrvNoExt($SCHbam_input);
		$SCHbam_output=$SCHbam_prefix.".CleanHeader.bam";
	}
	unlink ($SCHbam_output) if (-e $SCHbam_output);
##COMMENT: read bam/sam
	if ($SCHbam_input=~/sam$/i) {
		open (SCHBAMINPUT1, "samtools view -S $SCHbam_input|") || die "SUB(SamCleanHeader)Error: sam open error\n";
	}elsif ($SCHbam_input=~/bam$/i) {
		open (SCHBAMINPUT1, "samtools view $SCHbam_input|") || die "SUB(SamCleanHeader)Error: bam open error\n";
	}
	my %SCHseqID=();
	while (my $SCHline1=<SCHBAMINPUT1>) {
		my @SCHarr1=split(/\t/, $SCHline1);
		$SCHseqID{$SCHarr1[2]}++;
		if ($SCHline1=~/RG:Z:(\S+)/) {
#			print $1."\n";###For test###
			if (defined $1 and $1 ne '') {
				$SCHnew_readgroup=1;
				$SCHreahgroup{$1}++;
			}
		}
	}
	close SCHBAMINPUT1;
	print STDERR "SUB(SamCleanHeader)Info: total of ".scalar(keys %SCHseqID)." REF sequences detected\n" if ($debug);
	print STDERR "SUB(SamCleanHeader)Info: total of ".scalar(keys %SCHreahgroup)." readgroups detected\n" if ($debug);
##COMMENT: write output
	if ($SCHbam_input=~/sam$/i) {
		open (SCHBAMINPUT2, "samtools view -h -S $SCHbam_input|") || die "SUB(SamCleanHeader)Error: sam open error\n";
	}elsif ($SCHbam_input=~/bam$/i) {
		open (SCHBAMINPUT2, "samtools view -h $SCHbam_input|") || die "SUB(SamCleanHeader)Error: bam open error\n";
	}
	open (SCHBAMOUTPUT, "|samtools view -bS - > $SCHbam_output") || die "SUB(SamCleanHeader)Error: bam write error\n";
	if (defined $SCHoutnogroup) {
		open (SCHBAMOUTNOGROUP, "|samtools view -bS - > $SCHoutnogroup") || die "SUB(SamCleanHeader)Error: bam2 write error\n";
	}
	my $SCHtest_norg= (defined $SCHoutnogroup) ? 1 : 0;
	while (my $SCHline=<SCHBAMINPUT2>) {
		if ($SCHline=~m/^\@SQ\s+SN:(\S+)\s+LN:\d+/) {
			unless (exists $SCHseqID{$1}) {
				print STDERR "SUB(SamCleanHeader)Info: ignore $SCHline\n" if ($debug);
				next;
			}
		}
		if ($SCHline=~m/^\@RG/) {
			if (scalar(keys %SCHreahgroup) >0) {
				if ($SCHnew_readgroup==1) {
					foreach (keys %SCHreahgroup) {
						print SCHBAMOUTPUT '@RG'."\tID:$_\tSM:$_\tPL:Illumina\tLB:$_\n";
					}
					$SCHnew_readgroup=0;
				}
				else {
					next;
				}
			}
			else {
				print SCHBAMOUTPUT $SCHline;
				next;
			}
			if ($SCHtest_norg==0) {
				print SCHBAMOUTNOGROUP '@RG'."\tID:nogroup\tSM:nogroup\tPL:Illumina\tLB:nogroup\n" if (defined $SCHoutnogroup);
				$SCHtest_norg=1;
			}
		}
		else {
			$SCHtest_norg=0;
		}
		print SCHBAMOUTPUT $SCHline;
		$SCHline=~s/RG:Z:\S+/RG:Z:nogroup/g;
		print SCHBAMOUTNOGROUP $SCHline if ($SCHtest_norg==0 and defined $SCHoutnogroup);
	}
	close SCHBAMOUTPUT;
	close SCHBAMINPUT2;
	close SCHBAMOUTNOGROUP if (defined $SCHoutnogroup);
}

###Retrieve filebasename without extension
###& RetrvNoExt(file)
###Global:

sub RetrvNoExt {
	my $RNE_ori=shift @_;
	chomp $RNE_ori;
	my $RNE_new='';
	my $RNE_base='';
	($RNE_new=$RNE_base)=~s/^(\S+)\.\w+$/$1/;
	return $RNE_new;
}

