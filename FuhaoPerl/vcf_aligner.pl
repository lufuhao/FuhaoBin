#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FindBin qw($Bin);
use FuhaoPerl5Lib::VcfKit;
use FuhaoPerl5Lib::CmdKit;
use FuhaoPerl5Lib::MiscKit;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: LUFUHAO20150618

Requirements:
	Programs:
	Modiles: Scalar::Util, Cwd, Getopt::Long, FindBin

Descriptions:
	Determine the insert size given pairs of seqing data by
	mapping them to a reference.

Options:
	--help|-h
		Print this help/usage;
	--idlist           <fasta ID list>
	--reference|-f     <Reference.fa>
	--vcf1             AABBDD.vcf.gz
	--vcf2             AABB.vcf.gz
	--vcf3             AA.vcf.gz
	--vcf4             DD.vcf.gz
	--script           Script to extract BAM
	--output           Output.vcf
	--faillog          failed log
	--finallog         Final log
	--failregion       Failed region
	--problemvcf       Problematic.vcf
	--debug            Debug mode
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
my ($idlist, $reference, $vcf_aabbdd, $vcf_aabb, $vcf_aa, $vcf_dd);
my ($output, $faillog, $finallog, $failedregion, $problemvcf);
my ($bam_extract_sh);

GetOptions(
	"help|h!" => \$help,
	"idlist:s" => \$idlist,
	"reference|f:s" => \$reference,
	"vcf1:s" => \$vcf_aabbdd,
	"vcf2:s" => \$vcf_aabb, 
	"vcf3:s" => \$vcf_aa, 
	"vcf4:s" => \$vcf_dd,
	"script:s" => \$bam_extract_sh,
	"output|o:s" => \$output,
	"faillog:s" => \$faillog, 
	"finallog:s" => \$finallog, 
	"failregion:2" => \$failedregion,
	"problemvcf:" => \$problemvcf,
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$version) or die USAGE;
($help or $version) and die USAGE;



### Defaults ########################################################
$debug=0 unless (defined $debug);
$verbose=0 unless (defined $verbose);
my $path_tabix='tabix';
my $path_bgzip='bgzip';
my $freebayes_min_coverage=3;
my $freebayes_min_alternative_count=3;
my $min_mapq=0;
my @seq_ids;
my $path_freebayes='freebayes';
my $path_cdbfasta='cdbfasta';
my $vcfhead_AABBDD='AABBDD';
my $vcfhead_AABB='AABB';
my $vcfhead_AA='AA';
my $vcfhead_DD='DD';

### input and output ################################################
die "(IO)Error: invalid fasta ID file\n" unless (defined $idlist and -s $idlist);
die "(IO)Error: invalid AABBDD VCF file\n" unless (defined $vcf_aabbdd and -s $vcf_aabbdd);
die "(IO)Error: invalid AABB VCF file\n" unless (defined $vcf_aabb and -s $vcf_aabb);
die "(IO)Error: invalid AA VCF file\n" unless (defined $vcf_aa and -s $vcf_aa);
die "(IO)Error: invalid DD VCF file\n" unless (defined $vcf_dd and -s $vcf_dd);



###Prepare: 
my ($test_cdbfasta, $cdbfasta_index)=&CdbFasta($reference, $path_cdbfasta);
if (! $test_cdbfasta) {
	die "(IO)Error: index fasta failed: $test_cdbfasta\n";
}

### Main ############################################################
open (IDLIST, "< $idlist") || die "(IO)Error: can not open ID file: $idlist\n";
open (FAILLOG, "> $faillog") || die "(IO)Error: can not write fail log: $faillog\n";
open (FINALLOG, "> $finallog") || die "(IO)Error: can not write final log: $finallog\n";
open (REGIONFILE, ">$failedregion") || die "(IO)Error: can not write final log: $failedregion\n";
open (OUTVCF, ">$output") || die "(IO)Error: can not write VCF output: $output\n";
open (PROBLEMVCF, ">$problemvcf") || die "(IO)Error: can not write problematic VCF: $problemvcf\n";
my $num_lines=0;

while (my $line=<IDLIST>) {
	chomp $line;
	$num_lines++;
	@seq_ids=split(/\s+/, $line);
	my ($vcfobj1, $vcfobj2, $vcfobj3, $vcfobj4);
	
	if (! &ExtractVcf($vcf_aabbdd, $line, "clust$num_lines.AABBDD.vcf.gz", $path_tabix, $path_bgzip)) {
		print STDERR "Error:\tExtractAABBDD\t$num_lines\tfailed\n";
		next;
	}
	$vcfobj1=&LoadVcf("clust$num_lines.AABBDD.vcf.gz", 20, 1);
	
	if (! &ExtractVcf($vcf_aabb, $line, "clust$num_lines.AABB.vcf.gz", $path_tabix, $path_bgzip)) {
		print STDERR "Error:\tExtractAABB\t$num_lines\tfailed\n";
		next;
	}
	$vcfobj2=&LoadVcf("clust$num_lines.AABB.vcf.gz", 20, 1);
	
	if (! &ExtractVcf($vcf_aa, $line, "clust$num_lines.AA.vcf.gz", $path_tabix, $path_bgzip)) {
		print STDERR "Error:\tExtractAA\t$num_lines\tfailed\n";
		next;
	}
	$vcfobj3=&LoadVcf("clust$num_lines.AA.vcf.gz", 20, 1);
	if (! &ExtractVcf($vcf_dd, $line, "clust$num_lines.DD.vcf.gz", $path_tabix, $path_bgzip)) {
		print STDERR "Error:\tExtractDD\t$num_lines\tfailed\n";
		next;
	}
	$vcfobj4=&LoadVcf("clust$num_lines.DD.vcf.gz", 20, 1);
	
	close REGIONFILE if (defined fileno(REGIONFILE));

	
	my %posit2keep=();
	my @posit2rerun=();
	my @region=();
	my $rerun_freebayes=0;
	my %posit2run=();
	my %sum_alleles=();
	foreach my $chrom (keys %{$vcfobj1}) {
		my @region_startend=();
		foreach my $pos (keys %{${$vcfobj1}{$chrom}}) {
			my $rerun_aabbdd=0;
			my $hexa_allele='';
			my $rerun_aabb=0;
			my $tetra_allele='';
			my $rerun_aa=0;
			my $ura_allele='';
			my $rerun_dd=0;
			my $tau_allele='';
			my @aabbdd_geno=();
			my @aabbdd_alleles=();
			@region_startend=($pos, $pos);
			my %aabbdd_allele2num=();
			
			if (exists ${${${${$vcfobj1}{$chrom}}{$pos}}{'col10'}}{'GT'}) {
				$hexa_allele=${${${${$vcfobj1}{$chrom}}{$pos}}{'col10'}}{'GT'};
				if ($hexa_allele=~/\//) {
					@aabbdd_geno=&SplitGenotypes($hexa_allele, 2);
					if (exists ${${${$vcfobj1}{$chrom}}{$pos}}{'allele2num'}) {
						%aabbdd_allele2num=%{${${${$vcfobj1}{$chrom}}{$pos}}{'allele2num'}};
						@aabbdd_alleles=keys %aabbdd_allele2num;
						my $varmaxlength=&MaxLength(@aabbdd_alleles);
						$region_startend[1]+=$varmaxlength-1;
						if (scalar(@aabbdd_alleles)==1) {
							${$posit2keep{$chrom}}{$pos}=$region_startend[1];
							next;
						}
					}
					else {
						print FAILLOG "Failed\t$chrom\t$pos\tNoAllele2num\n";
						push (@posit2rerun, "$chrom:$region_startend[0]-$region_startend[1]");
						$rerun_aabbdd=1;
						next;
					}
				}
				elsif ($hexa_allele eq '.') {
					print FINALLOG "Failed\t$chrom\t$pos\tJust.\n";
					next;
				}
				else {
					print FAILLOG "Failed\t$chrom\t$pos\tNo/\n";
					push (@posit2rerun, "$chrom:$region_startend[0]-$region_startend[1]");
					$rerun_aabbdd=1;
					next;
				}
			}
			else {
				print FAILLOG "Failed\t$chrom\t$pos\tNoGT\n";
				push (@posit2rerun, "$chrom:$region_startend[0]-$region_startend[1]");
				$rerun_aabbdd=1;
				next;
			}
			
			
			AABB: {
				if (exists ${${${${$vcfobj2}{$chrom}}{$pos}}{'col10'}}{'GT'}) {
					if (${${${${$vcfobj2}{$chrom}}{$pos}}{'col10'}}{'GT'}=~/\//) {
						my @aabb_geno=&SplitGenotypes(${${${${$vcfobj2}{$chrom}}{$pos}}{'col10'}}{'GT'}, 2);
						my @aabb_alleles=();
						foreach (@aabb_geno) {
							push (@aabb_alleles, ${${${${$vcfobj2}{$chrom}}{$pos}}{'num2allele'}}{$_});
						}
						my @aabb_correct=();
						foreach (@aabb_alleles) {
							if (exists $aabbdd_allele2num{$_}) {
								push (@aabb_correct, $aabbdd_allele2num{$_});
							}
							else {
								push (@aabb_correct, $_);
								$rerun_aabb=1;
							}
						}
						$tetra_allele=join('/', @aabb_correct);
					}
					elsif (${${${${$vcfobj2}{$chrom}}{$pos}}{'col10'}}{'GT'} eq '.') {
						$tetra_allele='.';
						last AABB;
					}
					else {
						$rerun_aabb=1;
						last AABB;
					}
				}
				else {
					$rerun_aabb=1;
					last AABB;
				}
			}###AABB
			


			AA: {
				if (exists ${${${${$vcfobj3}{$chrom}}{$pos}}{'col10'}}{'GT'}) {
					if (${${${${$vcfobj3}{$chrom}}{$pos}}{'col10'}}{'GT'}=~/\//) {
						my @aa_geno=&SplitGenotypes(${${${${$vcfobj3}{$chrom}}{$pos}}{'col10'}}{'GT'}, 2);
						my @aa_alleles=();
						foreach (@aa_geno) {
							push (@aa_alleles, ${${${${$vcfobj3}{$chrom}}{$pos}}{'num2allele'}}{$_});
						}
						my @aa_correct=();
						foreach (@aa_alleles) {
							if (exists $aabbdd_allele2num{$_}) {
								push (@aa_correct, $aabbdd_allele2num{$_});
							}
							else {
								push (@aa_correct, $_);
								$rerun_aa=1;
							}
						}
						$ura_allele=join ('/', @aa_correct);
					}
					elsif (${${${${$vcfobj3}{$chrom}}{$pos}}{'col10'}}{'GT'} eq '.') {
						$ura_allele='.';
						last AA;
					}
					else {
						$rerun_aa=1;
						last AA;
					}
				}
				else {
					$rerun_aa=1;
					last AA;
				}
			}###AA
			
			

			DD: {
				if (exists ${${${${$vcfobj4}{$chrom}}{$pos}}{'col10'}}{'GT'}) {
					if (${${${${$vcfobj4}{$chrom}}{$pos}}{'col10'}}{'GT'}=~/\//) {
						my @dd_geno=&SplitGenotypes(${${${${$vcfobj4}{$chrom}}{$pos}}{'col10'}}{'GT'}, 2);
						my @dd_alleles=();
						foreach (@dd_geno) {
							push (@dd_alleles, ${${${${$vcfobj4}{$chrom}}{$pos}}{'num2allele'}}{$_});
						}
						my @dd_correct=();
						foreach (@dd_alleles) {
							if (exists $aabbdd_allele2num{$_}) {
								push (@dd_correct,$aabbdd_allele2num{$_});
							}
							else {
								push (@dd_correct, $_);
								$rerun_dd=1;
							}
						}
						$tau_allele=join('/', @dd_correct);
					}
					elsif (${${${${$vcfobj4}{$chrom}}{$pos}}{'col10'}}{'GT'} eq '.') {
						$tau_allele='.';
						last DD;
					}
					else {
						$rerun_dd=1;
						last DD;
					}
				}
				else {
					$rerun_dd=1;
					last DD;
				}
			}###DD
			if ($rerun_aabbdd==1 or $rerun_aabb==1 or $rerun_aa==1 or $rerun_dd==1) {
				$rerun_freebayes=1;
				my @sort_startend=sort {$a <=>$b} @region_startend;
				push (@region, "$chrom:$sort_startend[0]-$sort_startend[1]");
				if (exists ${$posit2run{$chrom}}{$sort_startend[0]}) {
					print STDERR "Warnings: rerun1: existing hash chrom:pos $chrom:$sort_startend[0]\n";
				}
				else {
					${$posit2run{$chrom}}{$sort_startend[0]}=$sort_startend[1];
				}
			}
			else {
				${$posit2keep{$chrom}}{$pos}=$region_startend[1];
				@{${$sum_alleles{$chrom}}{$pos}}=($hexa_allele, $tetra_allele, $ura_allele, $tau_allele);
			}
		}
	}
	
	if (scalar(@region) <1) {
		if ($rerun_freebayes==1) {
			print FAILLOG "Warnings\t$line\t???\tRunfreebayesWithoutRegion\n";
			$rerun_freebayes=0;
		}
	}
	
	if ($rerun_freebayes==1) {
		my $bamregion=join (' ', @region);
		unless (defined $bamregion and $bamregion=~/\S+/) {
			print FAILLOG "Warnings\t$line\t???\tEmptyBamRegion\n";
		}
		###Extract, merge, index bam
		if (! &exec_cmd_return("$bam_extract_sh -i $bamregion -o Clust.$num_lines.bam")) {
			print STDERR "Error:\tExtactBam\t$num_lines\tfailed\n";
			print FAILLOG "Warnings\t$line\t???\tExtractBam\n";
			print REGIONFILE "$bamregion\n";
		}
		###extract fasta
		if (! &CdbYank($cdbfasta_index, "Clust.$num_lines.fa", \@seq_ids, 'cdbyank')) {
			print STDERR "Error:\tcdbyank\t$num_lines\tfailed\n";
			print FAILLOG "Warnings\t$line\t???\tCdbyank\n";
			print REGIONFILE "$bamregion\n";
		}
		my $freebayes_addcmd=" --min-coverage $freebayes_min_coverage --min-alternate-count $freebayes_min_alternative_count --min-mapping-quality $min_mapq --pooled-discrete ";
		if (! &RunFreebayes("Clust.$num_lines.fa", @seq_ids, "Clust.$num_lines.bam", 3, "Clust.$num_lines.rerun.vcf.gz", 0, $freebayes_addcmd)) {
			print STDERR "Error:\tRunfreebayes\t$num_lines\tfailed\n";
			print FAILLOG "Warnings\t$line\t???\tRunfreebayes\n";
			print REGIONFILE "$bamregion\n";
		}
	}
	
	
	my %problemsites=();
	my %posit2keep_rerun=();
	my $vcfmerge=&LoadVcf("clust$num_lines.AABBDD.vcf.gz", 20, 1);
	foreach my $chrom (keys %{$vcfmerge}) {
		foreach my $pos (keys %{${$vcfmerge}{$chrom}}) {
			my $hexa_allele='';
			my $tetra_allele='';
			my $ura_allele='';
			my $tau_allele='';
			if (exists ${${${$vcfmerge}{$chrom}}{$pos}}{'allele2num'}) {
				my @alleles=keys %{${${${$vcfmerge}{$pos}}{$chrom}}{'allele2num'}};
				my $maxlength=&MaxLength(@alleles);
				if (exists $posit2keep{$chrom}) {
					my $test_overlap=0;
					foreach my $pos2 (keys %{$posit2keep{$chrom}}) {
						my $ovlapreturn=&TestIntersect($pos, $pos+$maxlength-1, $pos2, ${$posit2keep{$chrom}}{$pos2});
						if ($ovlapreturn eq '1') {
							$test_overlap=1;
							last;
						}
						elsif ($ovlapreturn eq '?') {
							print STDERR "WrongOverlap: $chrom ($pos, $pos+$maxlength-1, $pos2, ${$posit2keep{$chrom}}{$pos2})\n";
						}
					}
					next if ($test_overlap==1);
				}
				if (exists $posit2run{$chrom}) {
					my $test_overlap=0;
					foreach my $pos3 (keys %{$posit2run{$chrom}}) {
						my $ovlapreturn=&TestIntersect($pos, $pos+$maxlength-1, $pos3, ${$posit2keep{$chrom}}{$pos3});
						if ($ovlapreturn eq '1') {
							$test_overlap=1;
							last;
						}
						elsif ($ovlapreturn==1) {
							print STDERR "WrongOverlap2: $chrom ($pos, $pos+$maxlength-1, $pos3, ${$posit2keep{$chrom}}{$pos3})\n";
						}
					}
					next if ($test_overlap==0);
				}
				my $this_pos_is_problematic=0;
				PROBLEMATIC: { 
				if (exists ${${${$vcfmerge}{$chrom}}{$pos}}{$vcfhead_AABBDD}) {
					if (exists ${${${${$vcfmerge}{$chrom}}{$pos}}{$vcfhead_AABBDD}}{'GT'}) {
						$hexa_allele=exists ${${${${$vcfmerge}{$chrom}}{$pos}}{$vcfhead_AABBDD}}{'GT'};
						if (exists ${${${${$vcfmerge}{$chrom}}{$pos}}{$vcfhead_AABB}}{'GT'}) {
							$tetra_allele=${${${${$vcfmerge}{$chrom}}{$pos}}{$vcfhead_AABB}}{'GT'};
						}
						else {
							$this_pos_is_problematic=1;
							last PROBLEMATIC;
						}
						if (exists ${${${${$vcfmerge}{$chrom}}{$pos}}{$vcfhead_AA}}{'GT'}) {
							$ura_allele=${${${${$vcfmerge}{$chrom}}{$pos}}{$vcfhead_AA}}{'GT'};
						}
						else {
							$this_pos_is_problematic=1;
							last PROBLEMATIC;
						}
						if (exists ${${${${$vcfmerge}{$chrom}}{$pos}}{$vcfhead_DD}}{'GT'}) {
							$tau_allele=${${${${$vcfmerge}{$chrom}}{$pos}}{$vcfhead_DD}}{'GT'};
						}
						else {
							$this_pos_is_problematic=1;
							last PROBLEMATIC;
						}
						
					}
					else {
						$this_pos_is_problematic=1;
						last PROBLEMATIC;
					}
				}}###PROBLEMATIC
				if ($this_pos_is_problematic==1) {
					${$problemsites{$chrom}}{$pos}++;
				}
				else {
					if (exists ${$sum_alleles{$chrom}}{$pos}) {
						print STDERR "Error: existing sum_allele hash: chrom:pos $chrom:$pos\n";
					}
					else {
						@{${$sum_alleles{$chrom}}{$pos}}=($hexa_allele, $tetra_allele, $ura_allele, $tau_allele);
						${$posit2keep_rerun{$chrom}}{$pos}++;
					}
				}
			}
			else {
				print REGIONFILE "$chrom:$pos\tNoAllele2Num\n";
			}
		}
	}
	close ORIGINVCF if (defined fileno(ORIGINVCF));
	unless (open("clust$num_lines.AABBDD.vcf.gz") ) {
		print STDERR "Error: open Original VCF file: clust$num_lines.AABBDD.vcf.gz\n";
		next;
	}
	close RERUNVCF if (defined fileno(RERUNVCF));
	unless (open(RERUNVCF, "zcat Clust.$num_lines.rerun.vcf.gz |") ) {
		print STDERR "Error: open Original VCF file: Clust.$num_lines.rerun.vcf.gz\n";
		next;
	}
	my $num_original=0;
	my $num_rerun=0;
	
	
	while (my $line2=<ORIGINVCF>) {
		chomp $line2;
		next if ($line2=~/^#/);
		my @arr2=split(/\t/, $line2);
		if (exists ${$posit2keep{$arr2[0]}}{$arr2[1]}) {
			print OUTVCF $line2."\n";
			$num_original++;
		}
	}
	close ORIGINVCF;
	while (my $line3=<RERUNVCF>) {
		chomp $line3;
		next if ($line3=~/^#/);
		my @arr3=split(/\t/, $line3);
		my @new_arr3=();
		if (exists ${$posit2keep_rerun{$arr3[0]}}{$arr3[1]}) {
			for (my $col=0; $col<9; $col++) {
				push (@new_arr3, $arr3[$col]);
			}
		}
		if (exists ${${${${$vcfmerge}{$arr3[0]}}{$arr3[1]}}{$vcfhead_AABBDD}}{'genotypes'}) {
			push (@new_arr3, ${${${${$vcfmerge}{$arr3[0]}}{$arr3[1]}}{$vcfhead_AABBDD}}{'genotypes'});
		}
		else {
			print STDERR "Error: rerun AABBDD genotypes not found: chrom:pos $arr3[0]:$arr3[1]\n";
			next;
		}
		my $new_line3=join("\t", @new_arr3);
		if (exists ${$problemsites{$arr3[0]}}{$arr3[1]}) {
			print PROBLEMVCF $new_line3."\n";
			next;
		}
		else {
			print OUTVCF $new_line3."\n";
			$num_rerun++;
		}
	}
}

close REGIONFILE;
close FAILLOG;
close FINALLOG;
close IDLIST;
close OUTVCF;
close PROBLEMVCF;
#####################################################################
###                         sub functions                         ###
#####################################################################
#my $arr=&SplitGenotypes('1/0/1', 2);
#print "Arr: @{$arr}\n";
### ReadSam
### &ReadSam(genotype(1/0/1), ret_code(1/2/3))
### Global:
### Dependency:
### Note: 
### ret_code: 1=original array: 1 0 1; 2= non_redundant_array: 0 1
sub SplitGenotypes {
	my ($SGgeno, $SGret_code)=@_;
	my @SGarr=split(/\/|\|/, $SGgeno);
	if ($SGret_code==1) {
		return @SGarr;
	}
	my %SGgenohash=();
	foreach (@SGarr) {
		$SGgenohash{$_}++;
	}
	my @SGarr2=keys %SGgenohash;
	if ($SGret_code==2) {
		return @SGarr2;
	}
}
